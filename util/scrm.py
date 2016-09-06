
import bitarray
import numpy as np
import math
import os
import sh
from subprocess import check_output
import subprocess
import itertools
import logging
import sys
from collections import Counter
import re
import six
from six.moves import cStringIO as StringIO

import smcpp.util

logger = logging.getLogger(__name__)
scrm = sh.Command(os.environ['SCRM_PATH'])

# Stuff related to posterior decoding
# Format trees
def tree_obs_iter(l1, l2, trees):
    fs = frozenset([l1, l2])
    for sec, d, _ in trees:
        for i in range(sec):
            yield d[fs]

def true_hidden_states(trees, distinguished_lineages):
    tb = []
    M = len(hs) - 1
    for block in smcpp.util.grouper(tree_obs_iter(trees), block_size):
        a = np.zeros([M, 1])
        c = Counter(block)
        s = sum(c.values())
        for k in c:
            ip = np.searchsorted(hs, k) - 1
            a[ip] = 1. * c[k] / s
        tb.append(a)
    return np.array(tb).T

def splitter(_iter, key="//"):
    def f(line):
        if line.strip() == key:
            f.i += 1
        return f.i
    f.i = 0
    for (k, subiter) in itertools.groupby(_iter, f):
        yield (k, (line for line in subiter if line.strip() != key))

def demography_from_params(params):
    demography = []
    ct = 0.0
    z = list(zip(*params))
    for ai, bi, si in z[:-1]:
        beta = (np.log(ai) - np.log(bi)) / si
        demography += ['-eN', ct, ai]
        demography += ['-eG', ct, beta]
        ct += si
    demography += ['-eN', ct, z[-1][0]]
    demography += ['-eG', ct, 0.0]
    return demography

def print_data_stats(positions, haps):
    gaps = positions[1:] - positions[:-1]
    print(("Minimum gap: %d" % np.min(gaps)))
    print(("Average gap: %d" % np.mean(gaps)))
    print(("# seg. sites: %d" % gaps.shape[0]))
    i = np.argmin(gaps)
    print((positions[(i-3):(i+3)]))

def parse_scrm(n, L, output, include_trees):
    coal_times = []
    ts = 0
    for line in output:
        if line.startswith("segsites"):
            break
        if not include_trees:
            continue
        l = line.strip()
        k = l.index("(") - 1
        span = int(l[1:k])
        ts += span
        coal_times.append((span, l[(k+1):]))
    positions = next(output).strip()
    if positions:
        positions = (np.fromstring(positions[11:], sep=" ").astype('float32') * (L - 1)).astype('int')
        # ignore trailing newline
        haps = np.zeros([n, len(positions)], dtype=np.int8)
        # haps = []
        i = 0
        scrm_sfs = None
        for line in output:
            ls = line.strip()
            if not ls:
                continue
            if ls.startswith("SFS:"):
                scrm_sfs = [int(x) for x in ls[5:].split()]
            else:
                haps[i] = np.fromstring(str(line.strip()), np.uint8) - 48
                i += 1
                # .append(bitarray.bitarray(str(line).strip()))
                # print(len(haps))
        # haps = [bitarray.bitarray(str(line).strip()) for line in output if line.strip()] 
        uniqpos = np.concatenate(([True], positions[1:] != positions[:-1]))
        ret = (L, positions[uniqpos], haps[:, uniqpos])
        if scrm_sfs is not None:
            c = Counter(haps.sum(axis=0))
            psfs = [0] * (n - 1)
            for s in c:
                assert 0 < s < n
                psfs[s - 1] = c[s]
            if not np.array_equal(psfs, scrm_sfs):
                raise RuntimeError("sfs mismatch in parser")
        if include_trees:
            ret += (coal_times,)
        return ret
    return None


def simulate(n, N0, theta, rho, L, include_trees=False, scrm_args=[]):
    # scrm will emit positions in [0, L] (inclusive).
    seeds = np.random.randint(0, six.MAXSIZE, size=3)
    r = 4 * N0 * rho * (L - 1)
    t = 4 * N0 * theta * L
    args = [n, 1, '-p', int(math.log10(L)) + 2, '-t', t, '-r', r, L, '-seeds'] + list(seeds) + scrm_args
    if include_trees:
        args.append("-T")
    print(args)
    output = scrm(*args, _iter=True)
    cmd_line, seed, _, _ = [line.strip() for line in itertools.islice(output, 4)]
    ret = parse_scrm(n, L, output, include_trees)
    haps = ret[2]
    c = Counter(haps.sum(axis=0))
    psfs = [0] * (n - 1)
    for s in c:
        assert 0 < s < n
        psfs[s - 1] = c[s]
    return ret


def distinguished_sfs(n, M, N0, theta, demography, t0=0.0, t1=np.inf):
    seeds = np.random.randint(0, six.MAXSIZE, size=3)
    t = 4 * N0 * theta
    args = [n + 2, M, '-t', t, '-seeds'] + list(seeds) + demography
    if t0 > 0.0 or t1 < np.inf:
        args.append("-T")
    cmd = os.environ['SCRM_PATH'] + " " + " ".join(map(str, args))
    output = scrm(*args, _iter=True)
    avgsfs = np.zeros([3, n + 1], dtype=float)
    fs = frozenset([0, 1])
    m = 0
    for k, lines in splitter(output):
        if k == 0:
            continue
        sfs = np.zeros([3, n + 1], dtype=float)
        if t0 > 0.0 or t1 < np.inf:
            import smcpp._newick as _newick
            newick = next(lines)
            d12 = _newick.tmrca(newick, "1", "2")
            if not t0 < d12 < t1:
                continue
        m += 1
        ss = next(lines)
        segsites = int(ss.strip().split(" ")[1])
        if segsites == 0:
            avgsfs[0, 0] += 1
            continue
        next(lines) # positions
        bits = np.array([np.fromstring(str(line.strip()), np.uint8) - 48 for line in lines if line.strip()])
        assert (n + 2, segsites) == bits.shape
        for col in range(segsites):
            sfs[bits[:2, col].sum(), bits[2:, col].sum()] += 1
        avgsfs += sfs
    print((m, M))
    return avgsfs / m


def sfs(n, M, N0, theta, demography):
    t = 4 * N0 * theta
    args = [n + 2, M, '-t', t, '-oSFS'] + demography
    cmd = os.environ['SCRM_PATH'] + " " + " ".join(map(str, args))
    output = check_output(
            """%s | grep SFS | tail -n+2 | cut -f2- -d' ' | Rscript -e 'cat(colSums(read.table(file("stdin"))))'""" % cmd, shell=True)
    ret = np.array([float(x) for x in output.split()])
    assert ret.shape == (n + 1,)
    ret = np.append([M - ret.sum(),], ret)
    assert ret.sum() == M
    return ret / M

            
def empirical_transition(M, N0, rho, demography, hidden_states):
    r = 4 * N0 * rho
    args = [2, M, '-r', r, 2, '-T'] + demography
    output = scrm(*args).splitlines()
    ctre = re.compile(r'^\[.*\]\(\d:([^,]+),')
    H = hidden_states.shape[0] - 1
    M = np.zeros([H, H])
    h1 = []
    h2 = []
    for k, lines in splitter(output):
        if k == 0:
            continue
        l1 = next(lines)
        # Return in units of 4 * N0
        ct1 = float(ctre.match(l1).group(1)) * 4 * N0
        l2 = None
        try:
            l2 = next(lines)
            ct2 = float(ctre.match(l2).group(1)) * 4 * N0
        except:
            ct2 = ct1
        h1.append(float(ct1))
        h2.append(float(ct2))
    ch1, ch2 = np.searchsorted(hidden_states, [h1, h2])
    c = Counter(list(zip(ch1, ch2)))
    for a, b in c:
        M[a - 1, b - 1] = c[(a, b)]
    return M

def _pet_helper(args):
    return empirical_transition(*args)

def parallel_empirical_transition(*args):
    import multiprocessing as mp
    p = mp.Pool(32)
    res = list(p.imap(_pet_helper, [args for _ in range(16)]))
    p.terminate()
    p.close()
    del p
    return np.sum(res, axis=0)

if __name__ == "__main__":
    L = 10
    data = simulate(5, 2.0, 1e-8, 1e-8, L, ['-n', 1, 1])
    hmmd = hmm_data_format(data, (0, 1))
