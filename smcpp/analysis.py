import numpy as np
import functools
import json
import sys
import os.path
import scipy.optimize
import ad
import multiprocessing
import os
import concurrent.futures as futures

try:
    cpu_count = os.cpu_count()
except:
    cpu_count = multiprocessing.cpu_count()

def thread_pool():
    return futures.ThreadPoolExecutor(cpu_count)

from . import estimation_tools, _smcpp, util, logging, optimizer, jcsfs, spline
from .contig import Contig
from .model import SMCModel, SMCTwoPopulationModel

logger = logging.getLogger(__name__)

_model_cls_d = {cls.__name__: cls for cls in (SMCModel, SMCTwoPopulationModel)}

class BaseAnalysis:
    "Base class for analysis of population genetic data."
    def __init__(self, files, args):
        # Misc. parameter initialiations
        self._N0 = args.N0
        self._penalty = args.regularization_penalty
        self._niter = args.em_iterations
        args.solver_args = {}
        if args.factr:
            args.solver_args['factr'] = args.factr

        # Data-related stuff
        self._load_data(files)
        self._validate_data()
        self._recode_nonseg(args.nonseg_cutoff)


    ## PRIVATE INIT FUNCTIONS
    def _load_data(self, files):
        ## Parse each data set into an array of observations
        logger.info("Loading data...")
        self._files = files
        self._contigs = estimation_tools.load_data(files)
        pops = set(c.pid for c in self._contigs)
        unique_pops = list({x for p in pops for x in p})
        assert len(unique_pops) <= 2, (
                "Only one or two populations are supported, but the "
                "following were found in the data: %r" % unique_pops)
        assert len(unique_pops) <= 2
        self._populations = tuple(unique_pops)
        
        for c in self._contigs:
            assert len(c.n) == len(c.a)
            assert c.a.max() <= 2
            assert c.a.min() >= 0
            assert c.a.sum() == 2
            assert c.data.shape[1] == 1 + 3 * len(c.n)
            logger.debug("Contig(pid=%r, fn=%r, n=%r, a=%r)", c.pid, c.fn, c.n, c.a)
        logger.info("%d population%s", self.npop, "" if self.npop == 1 else "s")

    def _validate_data(self):
        for c in self._contigs:
            assert c.data.flags.c_contiguous
            bad = (np.all(c.data[:, 1::3] == c.a[None, :], axis=1) &
                   np.all(c.data[:, 2::3] == c.n[None, :], axis=1))
            if np.any(bad):
                logger.error("In file %s, observations %s:", c.fn, np.where(bad)[0])
                raise RuntimeError("Error: data set contains sites where every "
                        "individual is homozygous recessive. Please encode / "
                        "fold these as non-segregating (homozygous dominant).")

    def _recode_nonseg(self, cutoff):
        self._contigs = estimation_tools.recode_nonseg(self._contigs, cutoff)

    def _perform_thinning(self, thinning):
        # thin each dataset
        ns = np.array([sum(c.n) for c in self._contigs])
        if isinstance(thinning, int):
            thinning = np.array([thinning] * len(self._contigs))
        if thinning is None:
            thinning = 500 * ns
        if np.any(thinning > 1):
            logger.info("Thinning...")
            new_data = estimation_tools.thin_dataset(self._data, thinning)
            self._contigs = [Contig(data=d, pid=c.pid, fn=c.fn, n=c.n, a=c.a) 
                             for c, d in zip(self._contigs, new_data)]
        elif np.any(ns > 0):
            logger.warn("Not thinning yet undistinguished lineages are present")

    def _normalize_data(self, length_cutoff, no_filter):
        ## break up long spans
        self._contigs, attrs = estimation_tools.break_long_spans(self._contigs, length_cutoff)
        w, het = np.array([a[2:] for k in attrs for a in attrs[k]]).T
        self._het = avg = np.average(het, weights=w)
        if self._het == 0:
            logger.error("Data contain *no* mutations. Inference is impossible.")
            sys.exit(1)
        n = len(het)
        if n == 1:
            avg = 0.
            sd = np.inf
        else:
            var = np.average((het - avg) ** 2, weights=w) * (n / (n - 1.))
            sd = np.sqrt(var)
            logger.debug("Average/sd het:%f(%f)", avg, sd)
            if not no_filter:
                logger.debug("Keeping contigs within +-3 s.d. of mean")
        logger.debug("Average heterozygosity (derived / total bases) by data set (* = dropped)")
        ci = 0
        tpl = "%15d%15d%15d%12g"
        new_contigs = []
        for fn, key in zip(self._files, attrs):
            logger.debug(fn + ":")
            for attr in attrs[key]:
                het = attr[-1]
                mytpl = tpl
                if no_filter or abs(het - avg) <= 3 * sd:
                    new_contigs.append(self._contigs[ci])
                else:
                    mytpl += " *"
                logger.debug(mytpl % attr)
                ci += 1
        self._contigs = new_contigs

    def _init_hidden_states(self, prior_model, M):
        if prior_model is not None:
            d = json.load(open(prior_model, "rt"))
            model = _model_cls_d[d['model']['class']].from_dict(d['model'])
        else:
            model = self._model
        ## choose hidden states based on prior model
        dm = model.distinguished_model
        hs = estimation_tools.balance_hidden_states(dm, M)
        self._hidden_states = np.sort(
                np.unique(np.concatenate([self._model.distinguished_model._knots, hs]))
            )
        logger.debug("%d hidden states:\n%s" % (len(self._hidden_states), str(self._hidden_states)))

    def _init_inference_manager(self, fold):
        ## Create inference object which will be used for all further calculations.
        logger.debug("Creating inference manager...")
        d = {}
        self._ims = {}
        for c in self._contigs:
            k = (c.pid, tuple(c.n), tuple(c.a))
            d.setdefault(k, []).append(c)
        for pid, n, a in d:
            k = (pid, n, a)
            data = [contig.data for contig in d[k]]
            if len(pid) == 1:
                im = _smcpp.PyOnePopInferenceManager(n[0], data, 
                        self._hidden_states, k, fold)
            else:
                assert len(pid) == 2
                im = _smcpp.PyTwoPopInferenceManager(n[0], n[1], a[0], a[1], 
                        data, self._hidden_states, k, fold)
            im.model = self._model
            im.theta = self._theta
            im.rho = self._rho
            self._ims[k] = im
        self._model.randomize()

    def _init_bounds(self, Nmin):
        ## Construct bounds
        # P(seg) is at most theta * 2 * N_max / H_n << 1
        # For now, just base bounds off of population 1.
        sample_size = 2 + max(sum(c.n) for c in self._contigs)
        Hn = np.log(sample_size)
        Nmax = 1 / (2 * self._theta * Hn) / 5.
        logger.debug("Nmax calculated to be %g" % Nmax)
        self._bounds = (Nmin, Nmax)

    @property
    def _data(self):
        return [c.data for c in self._contigs]

    def run(self):
        'Perform the analysis.'
        self._optimizer.run(self._niter)

    def Q(self, k=None):
        'Value of Q() function in M-step.'
        # q1, q2, q3 = self._im.Q(True)
        qq = 0.
        with thread_pool() as executor:
            fs = []
            for na in self._ims:
                fs.append(executor.submit(self._ims[na].Q))
            for x in futures.as_completed(fs):
                qq += x.result()
        qr = -self._penalty * self.model.regularizer()
        logger.debug(("Q", float(qq), [qq.d(x) for x in self.model.dlist]))
        logger.debug(("reg", float(qr), [qr.d(x) for x in self.model.dlist]))
        return qq + qr

    def E_step(self):
        'Perform E-step.'
        logger.info('Running E-step')
        with thread_pool() as executor:
            fs = []
            for na in self._ims:
                fs.append(executor.submit(self._ims[na].E_step))
            futures.wait(fs)
        logger.info('E-step completed')

    def loglik(self):
        'Log-likelihood of data after most recent E-step.'
        ll = 0
        with thread_pool() as executor:
            fs = []
            for na in self._ims:
                fs.append(executor.submit(self._ims[na].loglik))
            for x in futures.as_completed(fs):
                ll += x.result()
        return ll - self._penalty * float(self.model.regularizer())

    @property
    def model(self):
        return self._model

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, r):
        self._rho = r
        for im in self._ims.values():
            im.rho = r

    @property
    def N0(self):
        return self._N0

    @property
    def npop(self):
        'The number of populations contained in this analysis.'
        return len(self._populations)

    def dump(self, filename):
        'Dump result of this analysis to :filename:.'
        d = {'N0': self._N0, 'theta': self._theta, 'rho': self._rho}
        d['model'] = self.model.to_dict()
        json.dump(d, open(filename + ".json", "wt"), sort_keys=True, indent=4)


class Analysis(BaseAnalysis):
    '''A dataset, model and inference manager to be used for estimation.'''
    def __init__(self, files, args):
        BaseAnalysis.__init__(self, files, args)

        # Perform initial filtering for weird contigs
        self._normalize_data(args.length_cutoff, args.no_filter)

        # Initialize members
        self._init_parameters(args.theta, args.rho)
        self._init_bounds(args.Nmin)
        self._init_model(args.pieces, args.N0, args.t1,
                args.tK, args.offset, args.knots, args.spline)

        if not args.no_initialize:
            self._hidden_states = np.array([0., np.inf])
            self._init_inference_manager(args.fold)
            self._init_optimizer(args, files, args.outdir, args.block_size,
                    args.algorithm, args.tolerance, learn_rho=False)
            self._optimizer.run(1)

        # Thin the data
        self._perform_thinning(args.thinning)

        # Continue initializing
        self._init_hidden_states(args.prior_model, args.M)
        self._init_inference_manager(args.fold)
        self._init_optimizer(args, files, args.outdir, args.block_size,
                args.algorithm, args.tolerance, learn_rho=True)

    def _init_parameters(self, theta=None, rho=None):
        ## Set theta and rho to their default parameters
        self._L = sum([d[:,0].sum() for d in self._data])
        logger.info("%.2f Gb of data", self._L * 1e-9)
        if theta is not None:
            self._theta = theta
        else:
            # Compute watterson's estimator while not really accounting
            # for any sort of population structure or missing data.
            # TODO This could be improved.
            watt = []
            for c in self._contigs:
                ma = (c.a[None, :] * (c.data[:, 1::3] != -1)).sum(axis=1)
                sample_size = ma + c.data[:, 3::3].sum(axis=1)
                nseg = np.maximum(0, c.data[:, 1::3]).sum(axis=1) + c.data[:, 2::3].sum(axis=1) 
                ss0 = sample_size > 0
                span = c.data[ss0, 0]
                nseg0 = nseg[ss0] > 0
                watt.append((np.average(nseg0 / np.log(sample_size[ss0]), weights=span), len(ss0)))
                assert(np.isfinite(watt[-1][0]))
            x, w = np.array(watt).T
            self._theta = np.average(x, weights=w)
        logger.info("theta: %f", self._theta)
        self._rho = rho or self._theta
        assert np.all(np.isfinite([self._rho, self._theta]))
        logger.info("rho: %f", self._rho)

    def _init_model(self, pieces, N0, t1, tK, offset, knots, spline_class):
        ## Initialize model
        # FIXME currently disabled.
        # exponential_pieces = args.exponential_pieces or []
        pieces = estimation_tools.extract_pieces(pieces)
        fac = 2. * N0
        t1 /= fac
        tK /= fac
        time_points = estimation_tools.construct_time_points(t1, tK, pieces, 0.)
        logger.debug("time points in coalescent scaling:\n%s", str(time_points))
        try:
            num_knots = int(knots)
            knot_spans = np.ones(num_knots, dtype=int)
            knots = np.cumsum(estimation_tools.construct_time_points(t1, tK, knot_spans, offset))
        except ValueError:
            knots = [float(x) for x in knots.split(",")]
        logger.debug("knots in coalescent scaling:\n%s", str(knots))
        spline_class = {"cubic": spline.CubicSpline,
                        "bspline" : spline.BSpline,
                        "akima": spline.AkimaSpline, 
                        "pchip": spline.PChipSpline}[spline_class]
        y0 = self._het / (2. * self._theta)
        logger.debug("Long term avg. effective population size: %f", y0)
        y0 = np.log(y0)
        if self.npop == 1:
            self._model = SMCModel(time_points, knots, spline_class, self._populations[0])
            self._model[-1] = y0
        else:
            split = tK - t1  # just pick the midpoint as a starting value.
            split /= 2. * N0
            mods = []
            for pid in self._populations:
                mods.append(SMCModel(time_points, knots, spline_class, pid))
                mods[-1][-1] = y0
            self._model = SMCTwoPopulationModel(mods[0], mods[1], split)

    def _init_optimizer(self, args, files, outdir, block_size, algorithm, tolerance, learn_rho):
        if self.npop == 1:
            self._optimizer = optimizer.SMCPPOptimizer(
                self, algorithm, tolerance, block_size, args.solver_args)
            # Also optimize knots in 1 pop case. Not yet implemented
            # for two pop case.
            # self._optimizer.register(optimizer.KnotOptimizer())
        elif self.npop == 2:
            self._optimizer = optimizer.TwoPopulationOptimizer(
                self, algorithm, tolerance, block_size, args.solver_args)
            smax = np.sum(self._model.distinguished_model.s)
            self._optimizer.register(
                optimizer.ParameterOptimizer("split", (0., smax), "model"))
        self._optimizer.register(optimizer.AnalysisSaver(outdir))
        if learn_rho:
            self._optimizer.register(
                optimizer.ParameterOptimizer("rho", (1e-6, 1e-2)))

    ## END OF PRIVATE FUNCTIONS
    @property
    def bounds(self):
        return self._bounds


class SplitAnalysis(BaseAnalysis):
    def __init__(self, files, args):
        BaseAnalysis.__init__(self, files, args)
        assert self.npop == 2
        self._init_model(args.pop1, args.pop2)
        self._init_bounds(.001)

        self._hidden_states = np.array([0., np.inf])
        self._init_inference_manager(False)
        self._init_optimizer(args, files, args.outdir, args.algorithm, args.tolerance, args.block_size, False)
        # Hack to only estimate split time.
        self._optimizer.run(1)

        # After inferring initial split time, thin
        self._perform_thinning(args.thinning)
        self._normalize_data(args.length_cutoff, args.no_filter)

        self._init_hidden_states(args.pop1, args.M)
        self._init_inference_manager(False)
        self._init_optimizer(args, files, args.outdir, args.algorithm, args.tolerance, args.block_size)

    def _init_optimizer(self, args, files, outdir, algorithm, tolerance, block_size, save=True):
        self._optimizer = optimizer.TwoPopulationOptimizer(self, algorithm, tolerance, block_size, args.solver_args)
        smax = np.sum(self._model.distinguished_model.s)
        self._optimizer.register(optimizer.ParameterOptimizer("split", (0., smax), "model"))
        if save:
            self._optimizer.register(optimizer.AnalysisSaver(outdir))

    def _init_model(self, pop1, pop2):
        d = json.load(open(pop1, "rt"))
        self._theta = d['theta']
        self._rho = d['rho']
        m1 = _model_cls_d[d['model']['class']].from_dict(d['model'])
        d = json.load(open(pop2, "rt"))
        m2 = _model_cls_d[d['model']['class']].from_dict(d['model'])
        self._model = SMCTwoPopulationModel(m1, m2, np.sum(m1.s) * 0.5)
