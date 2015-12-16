#/usr/bin/env python2.7
# Complete example showing how to use the package for inference
from __future__ import division
import numpy as np
import scipy.optimize
import pprint
import multiprocessing
import sys
import itertools
from collections import Counter
import sys
import time
import argparse
import psmcpp._pypsmcpp
from psmcpp.lib.util import config2dict
import ConfigParser as configparser
import cPickle as pickle

np.set_printoptions(linewidth=120, precision=6, suppress=True)

def exp_quantiles(M, h_M):
    hs = -np.log(1. - np.linspace(0, h_M, M, False) / h_M)
    hs = np.append(hs, h_M)
    hs[0] = 0
    return hs

parser = argparse.ArgumentParser("smc++")
parser.add_argument("--debug", action="store_true", default=False, help="display a lot of debugging info")
parser.add_argument("--comment", default=None, type=str)
parser.add_argument('config', type=argparse.FileType('r'), help="config file")
parser.add_argument('data', type=argparse.FileType('rb'), help="data file in smcpp format")
args = parser.parse_args()

psmcpp._pypsmcpp.do_progress(args.debug)
smcpp_data = pickle.load(args.data)
obs_list = smcpp_data['obs']
n = smcpp_data['n']
config = configparser.SafeConfigParser()
config.readfp(args.config)
print(config2dict(config))

## Calculate observed SFS for use later
osfs = []
for ol0 in obs_list:
    obsfs = np.zeros([3, n - 1])
    for r, c1, c2 in ol0[ol0[:, 1:].sum(axis=1) > 0]:
        obsfs[c1, c2] += r
    obsfs /= ol0[:, 0].sum()
    obsfs[0, 0] = 1. - obsfs.sum()
    osfs.append(obsfs)
obsfs = np.mean(osfs, axis=0)
print(" - Observed sfs:")
print(obsfs)

# Emission mask
em = np.arange(3 * (n - 1), dtype=int).reshape([3, n - 1])
# em[0, 3:] = 3
# em[1] = 4
# em[2] = 5

# Model parameters
mbs = config.getfloat('advanced', 'minimum block size')

try:
    ts = np.array(eval(config.get('model parameters', 'ts')))
    s = ts[1:] - ts[:-1]
    K = len(s)
except configparser.NoOptionError:
    t_1 = config.getfloat('model parameters', 't_1')
    t_K = config.getfloat('model parameters', 't_K')
    Ks = config.get('model parameters', 'K').split("+")
    pieces = []
    for piece in Ks:
        try:
            num, span = list(map(int, piece.split("*")))
        except ValueError:
            span = int(piece)
            num = 1
        pieces += [span] * num
    s = np.logspace(np.log10(t_1), np.log10(t_K), sum(pieces))
    s = np.concatenate(([t_1], s[1:] - s[:-1]))
    sp = np.zeros(len(pieces))
    count = 0
    for i, p in enumerate(pieces):
        sp[i] = s[count:(count+p)].sum()
        count += p
    s = sp
print("time points", s)
print(np.cumsum(s))

## Compute hidden states
try:
    hs = np.array(eval(config.get('hidden states', 'hidden states')))
except configparser.NoOptionError:
    M = config.getint('hidden states', 'M')
    h_M = config.getfloat('hidden states', 'h_M')
    hs = exp_quantiles(M, h_M)
if hs[0] != 0:
    raise Exception("First hidden state interval must begin at 0")
print("hidden states", hs)

# Load additional params
N0 = config.getfloat('parameters', 'N0')
mu = config.getfloat('parameters', 'mu')
rho = config.getfloat('parameters', 'rho')
block_size = config.getint('advanced', 'block size')

t_start = time.time()
try:
    thinning = config.getint('advanced', 'thinning')
except configparser.NoOptionError:
    thinning = n

try:
    lambda_penalty = config.getfloat("advanced", "lambda penalty")
except configparser.NoOptionError:
    lambda_penalty = 0.0

im = psmcpp._pypsmcpp.PyInferenceManager(n - 2, obs_list, hs,
        4.0 * N0 * mu, 4.0 * N0 * rho,
        block_size, thinning, em)

try:
    im.hj = config.getboolean('advanced', 'use hj')
except configparser.NoOptionError:
    im.hj = True
print("using hj: {hj}".format(hj=im.hj))

K = len(s)
x0 = np.ones([2, K])
a, b = x0
flat = not config.getboolean('model parameters', 'piecewise exponential')
print("Using piecewise exponential: {notflat}".format(notflat=not flat))
if flat:
    b = a
else:
    b += 0.1

im.setParams((a,b,s),False)
im.Estep()

llold = -np.inf
bounds = np.array([[0.1, 20.0]] * K + [[0.15, 19.9]] * K).reshape([2, K, 2])

# Optimization part
precond = 1. / s
precond[-1] = 1. / (15.0 - np.sum(s))

def optimize_fullgrad(iter, coords, x0, factr=1e9):
    print("Optimizing factr {factr}".format(factr=factr))
    def fprime(x):
        x0c = x0.copy()
        # Preconditioner (sort of)
        for xx, cc in zip(x, coords):
            x0c[cc] = xx * precond[cc[1]]
        global s
        aa, bb = x0c
        if flat:
            bb = aa
        print(aa)
        print(bb)
        print(s)
        im.setParams((aa, bb, s), coords)
        print("done")
        res = im.Q(lambda_penalty)
        lls = np.array([ll for ll, jac in res])
        jacs = np.array([jac for ll, jac in res])
        ret = [-np.mean(lls, axis=0), -np.mean(jacs, axis=0)]
        for i, cc in enumerate(coords):
            ret[1][i] *= precond[cc[1]]
        if flat:
            print(ret[1])
        else:
            print(ret[1].reshape([2,ret[1].shape[0] // 2]))
        print(ret[0])
        # reg = im.regularizer()
        # print("regularizer: ", LAMBDA_PENALTY * reg)
        # add penalty
        # esfs, jac = psmcpp._pypsmcpp.sfs(n, (aa, bb, s), 0.0, hs1[-1], 4 * N0 * theta, coords)
        # diff = esfs[0, 0] - obsfs[0, 0]
        # penalty = ALPHA_PENALTY * diff**2
        # print("penalty", penalty)
        # ret[0] += penalty
        # ret[1] += 2 * ALPHA_PENALTY * diff * jac[0, 0]
        # print(x)
        # print(ret[0])
        # print(ret[1])
        # print("------")
        return ret
    # print("gradient check")
    # xx0 = np.array([x0[cc] / precond[cc[1]] for cc in coords])
    # f0, fp = fprime(xx0)
    # for i, cc in enumerate(coords):
    #     x0c = xx0.copy()
    #     x0c[i] += 1e-8
    #     f1, _ = fprime(x0c)
    #     print(i, cc, f1, f0, (f1 - f0) / 1e-8, fp[i])
    res = scipy.optimize.fmin_l_bfgs_b(fprime, [x0[cc] / precond[cc[1]] for cc in coords], 
            None, bounds=[tuple(bounds[cc] / precond[cc[1]]) for cc in coords], disp=False, factr=factr)
    # print(res)
    return np.array([x * precond[cc[1]] for x, cc in zip(res[0], coords)])

break_loop = False
import signal, sys
def print_state():
    global a, b, s
    d = {'a': a, 'b': b, 's': s, 'argv': sys.argv, 't_start': t_start, 't_now': time.time(), 'config': config2dict(config), 'comment': args.comment}
    pprint.pprint(repr(d))
    return d

def signal_handler(signal, frame):
    print("State...")
    print_state()
# signal.signal(signal.SIGINT, signal_handler)

def run_iteration(i, coords, factr):
    global x0
    global llold
    global im
    # for j in range(K * di // 3, K * (di + 1) // 3)]
    ret = optimize_fullgrad(i, coords, x0, factr)
    for xx, cc in zip(ret, coords):
        x0[cc] = xx
    print("************** ITERATION %d ***************" % i)
    print(a)
    print(b)
    im.setParams((a, b, s), False)
    im.Estep()
    ll = np.sum(im.loglik(0.0))
    print(" - New loglik:" + str(ll))
    print(" - Old loglik:" + str(llold))
    if ll < llold:
        print("*** Log-likelihood decreased")
    if llold == -np.inf:
        ret = 1.
    else:
        ret = (llold - ll) / llold
    llold = ll
    return ret

i = 0
ca = [0]
if not flat:
    ca.append(1)
while i < 20:
    coords = [(aa, j) for aa in ca for j in range(K)]
    # run_iteration(i, coords, (10**(12 - i / 10.)))
    run_iteration(i, coords, 1e9)
    esfs = psmcpp._pypsmcpp.sfs(n, (a,b,s), 0.0, hs[-1], 4 * N0 * mu, False)
    print("calculated sfs")
    print(esfs)
    print("observed sfs")
    print(obsfs)
    i += 1
d = print_state()
