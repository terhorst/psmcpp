import multiprocessing
from logging import getLogger
import signal
import collections
import wrapt
import ad
logger = getLogger(__name__)
# Package imports
from . import _smcpp
from .population import Population


@wrapt.decorator
def _fix_derivatives(wrapped, instance, args, kwargs):
    # Here we must patch up derivatives. Because they are being
    # serialized the objects will not be equal so calls to
    # x.d(v) will not produce the correct answers.
    models = args[0]
    tds = []
    for m in models:
        tds.append({})
        for x in m.x.flat:
            if isinstance(x, ad.ADF):
                for d in x.d():
                    if hasattr(d, 'tag') and d.tag is not None:
                        tds[-1][d.tag] = d
    ret = wrapped(*args, **kwargs)
    for r, td in zip(ret, tds):
        if isinstance(r, ad.ADF):
            r = [r]
        for rr in r:
            dr = rr.d()
            keys = list(dr)
            for k in keys:
                new_k = td.get(k.tag, dr[k])
                val = dr[k]
                del dr[k]
                dr[new_k] = val
    return ret


class Worker(multiprocessing.Process):
    def __init__(self, pipe, population):
        multiprocessing.Process.__init__(self)
        self._pipe = pipe
        self._population = population

    def run(self):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        self._population = Population(*self._population)
        while True:
            task, args = self._pipe.recv()
            # logger.debug((task, args))
            if task == "exit":
                # logger.debug("exiting")
                self._pipe.send(True)
                self._pipe.close()
                break
            f = getattr(self._population, task)
            self._pipe.send(f(*args))

class InferenceService(object):
    def __init__(self, populations):
        '''Initialize the inference service with a sequence of populations. 
        Each population consists of group of data sets.'''
        # Initialize workers
        self._parent_pipes, self._child_pipes = list(zip(*[multiprocessing.Pipe() for _ in populations]))
        self._workers = [Worker(pipe, pop) for pipe, pop in zip(self._child_pipes, populations)]
        self._npop = len(populations)
        logger.debug("starting workers")
        for worker in self._workers:
            worker.start()
        logger.debug("finished initializing workers")

    def _send_message(self, message, args=None):
        try:
            if args is None:
                args = [[]] * self._npop
            for p, a in zip(self._parent_pipes, args):
                p.send((message, a))
            return [p.recv() for p in self._parent_pipes]
        except KeyboardInterrupt:
            self.close()
            raise

    def __del__(self):
        self.close()

    def close(self):
        for p in self._parent_pipes:
            p.send(("exit", None))
        self._parent_pipes = []

    @_fix_derivatives
    def Q(self, models, coords):
        self.set_params(models, coords)
        return self._send_message("Q")

    def E_step(self):
        return self._send_message("E_step")

    @_fix_derivatives
    def penalize(self, models):
        return self._send_message("penalize", [[m] for m in models])

    def set_params(self, models, coords):
        coords = [coords] * len(models)
        return self._send_message("set_params", list(zip(models, coords)))

    def loglik(self):
        return self._send_message("loglik")

    def dump(self, file):
        self._send_message("dump", file)

    @property
    def coords(self):
        return self._send_message("coords")

    @property
    def sfs(self):
        return self._send_message("sfs")

    @property
    def precond(self):
        return self._send_message("precond")

    @property
    def theta(self):
        return self._send_message("theta")

    @property
    def model(self):
        return self._send_message("model")

# Used for debugging, does not fork()
class DumbInferenceService(InferenceService):
    def __init__(self, populations):
        '''Initialize the inference service with a sequence of populations. 
        Each population consists of group of data sets.'''
        # Initialize workers
        self._populations = [Population(*pop) for pop in populations]
        self._npop = len(populations)

    def _send_message(self, message, args=None):
        if args is None:
            args = [[]] * self._npop
        return [getattr(p, message)(*a) for p, a in zip(self._populations, args)]

    def close(self):
        pass