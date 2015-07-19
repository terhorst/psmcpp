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
import bitarray

import matplotlib
matplotlib.use('Agg')
import psmcpp.scrm, psmcpp.bfgs, psmcpp._pypsmcpp, psmcpp.util, psmcpp.plotting

def norm_counter(c, nn): 
    s = sum(c.values()); 
    d = {k:np.array([1.*v/s]) for k,v in c.items()}
    ret = np.zeros([3, nn - 1])
    for k in d:
        ret[k] = d[k]
    return ret

def p_emit(obs, sfs):
    return (obs[None, ...] * np.log(sfs)).sum(axis=(1, 2))

def unpack(iterable):
    for span, x in iterable:
        for i in range(span):
            yield x

def pack(seq):
    iterable = iter(seq)
    x = next(iterable)
    i = 1
    for xp in iterable:
        if xp == x:
            i += 1
        else:
            yield (i, x)
            x = xp
            i = 1
    yield (i, x)

def subset_data(d, start, end):
    pnew = np.logical_and(start <= d[1], d[1] < end)
    inew = np.where(pnew)[0]
    return (int(end - start), d[1][pnew] - start,
            [bitarray.bitarray([h[i] for i in inew]) for h in d[2]],
            list(pack(itertools.islice(unpack(d[3]), start, end))))

num_threads = 2
block_size = 25
num_samples = 100
np.set_printoptions(linewidth=120, precision=6, suppress=True)
N0 = 10000
rho = 1e-8
theta = 2.5e-8
L = int(float(sys.argv[1]))
a = np.array([10., 1., .1, 2.])
b = np.array([1., 1., .1, 2.])
s = np.array([5000.0, 50000.0, 10000., 10000.]) / 25.0 / (2 * N0)
true_parameters = (a, b, s)
width = 2000
M = 10

nns = [2, 3, 5, 10, 25]
n = max(nns)
demography = psmcpp.scrm.demography_from_params((a * 2.0, b * 2.0, s))
print(" ".join(map(str, demography)))
full_data = psmcpp.scrm.simulate(n, N0, theta, rho, L, demography, include_trees=True, seed=int(sys.argv[2]))

start = 0
end = L
mid = (end - start) // 2
span = mid / 2
plot_start = int((mid - span / 2) // block_size)
plot_end = int((mid + span / 2) // block_size)
ma_window = int((plot_end - plot_start) / width)
print(start, end, mid, plot_start, plot_end)
data = subset_data(full_data, start, end)
ct = [(c1, psmcpp.scrm.newick_to_dists(c2, "123")) for c1, c2 in data[3]]

pd = {}
gm = {}
al = {}
bt = {}
bs = {}
true_pos = {}
ims = {}
oo = {}
oo_orig = {}
for nn in nns:
    # subset data. some sites might be non-segregating in the subsample.
    seg = [i for i, pos in enumerate(data[1]) if any(not h[i] for h in data[2][:nn])]
    segdata = [bitarray.bitarray(ba[i] for i in seg) for ba in data[2][:nn]]
    dsub = (data[0], data[1][seg], segdata, data[3])
    obs = psmcpp.scrm.hmm_data_format(dsub, (0, 1))
    oo_orig[nn] = obs.copy()
    ol = obs
    ol[ol[:, 1] == 1, 2] = 0
    ol[np.logical_and(ol[:, 1] == 0, ol[:, 2] > 1), 2] = 0
    ol[ol[:, 1] == 2, 1:] = [0, 0]
    oo[nn] = ol
    hidden_states = np.array([0., np.inf])
    im = psmcpp._pypsmcpp.PyInferenceManager(nn - 2, [obs], hidden_states,
            4.0 * N0 * theta, 4.0 * N0 * rho,
            block_size, num_threads, num_samples, 10, [0], None)
    hidden_states = im.balance_hidden_states((a, b, s), M)
    print("balanced hidden states", hidden_states)
    em = np.arange(3 *  (nn - 1), dtype=int).reshape([3, nn - 1])
    em[0] = em[2] = 0
    em[1] = 1
    print(em)
    # em = np.zeros([3, nn - 1], dtype=int)
    # em[0] = 0
    # em[1] = 1
    # em[2] = 2
    # em[0, 1:2] = 3
    # em[1, 1:2] = 4
    im = psmcpp._pypsmcpp.PyInferenceManager(nn - 2, [obs], hidden_states,
            4.0 * N0 * theta, 4.0 * N0 * rho,
            block_size, num_threads, num_samples, 5, [0], em)
    im.seed = np.random.randint(0, sys.maxint)
    im.setParams((a, b, s), False)
    im.Estep()
    ims[nn] = im
    gamma = im.gammas()[0][:, plot_start:plot_end]
    alpha = im.alphas()[0][:, plot_start:plot_end]
    beta = im.betas()[0][:, plot_start:plot_end]
    Bs = im.Bs()[0][:, plot_start:plot_end]
    beta /= beta.sum(axis=0)
    gm[nn] = np.zeros([M, width])
    al[nn] = np.zeros([M, width])
    bt[nn] = np.zeros([M, width])
    bs[nn] = np.zeros([M, width])
    for i in range(width):
        for dct, ary in (gm, gamma), (al, alpha), (bt, beta), (bs, Bs):
            dct[nn][:, i] = np.mean(ary[:, (i * ma_window):((i + 1) * ma_window)], axis=1)
    # pd[nn] = psmcpp.inference.posterior_decode_score(0, 1, block_size, hidden_states, gamma, data[3])
    # Only care about first since they are all the same

import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=1 * len(nns) + 1, sharex=True, sharey=True, figsize=(30,15))
true_pos = np.zeros([2, width])
for i, ll in enumerate(((0, 1), (1, 2))[:2]):
    coal_times = np.searchsorted(hidden_states, 
            [x[frozenset(ll)] 
                for x in itertools.islice(unpack(ct), plot_start * block_size, plot_end * block_size)]) - 1
    for j in range(width):
        true_pos[i, j] = np.mean(coal_times[(j * block_size * ma_window):((j + 1) * block_size * ma_window)])
    axes[-1].step(range(width), true_pos[i])
plt.set_cmap("cubehelix")
ai = 0
for dct in (gm,):
    for nn in sorted(dct):
        ax = axes[ai]
        ai += 1
        im = ax.imshow(dct[nn][::-1], extent=[0,width,-0.5,M-0.5],aspect='auto', vmin=0.0, vmax=1.0)
        ax.step(range(width), true_pos[0], color=(0, 1., 1.))
        ax.set_ylabel("n=%d" % nn)
        ax.set_ylim([-0.5, M - 0.5])
axes[-1].set_ylabel("True hid. st.")
axes[-1].set_ylim([-0.5, M - 0.5])
axes[-1].xaxis.set_ticks(np.arange(0, width + 1, 100))
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.92, 0.29, 0.02, 0.55])
fig.colorbar(im, cax=cbar_ax)
psmcpp.plotting.save_pdf(fig, "posterior_decoding_heatmap_%d_%d.pdf" % (start, end))
plt.close(fig)

aoeu

fig, ax = psmcpp.plotting.pretty_plot()
for nn, color in zip(pd, psmcpp.plotting.palette[1:]):
    print(n, color)
    pd = np.mean([moving_average(test_posterior(n)) for _ in range(5)], axis=0)
    ax.plot(pd[::100], color=color, label=n)
psmcpp.plotting.save_pdf(fig, "posterior_decoding.pdf")

