from __future__ import print_function, division
import numpy as np
import smcpp._smcpp, smcpp.model, smcpp.spline, smcpp.util, smcpp.estimation_tools
from smcpp.logging import setup_logging, getLogger
import scipy.optimize
import sys
import os.path

logger = getLogger(__name__)

# Jonathan's test for random data
def test_stitch():
    hs = np.concatenate([[0.0], np.logspace(-3, 1, 30), [np.inf]])
    s = np.diff(hs)[:10]
    K = 10
    model = smcpp.model.SMCModel(s, np.logspace(-2, 0, 5))
    model.randomize()
    n = 30
    fakeobs = [[1, -1, 0, 0], [1, 1, 0, 0], [10, 0, 0, 0], [10, -1, 0, 0],
               [200000, 0, 0, n - 2], [1, 1, n - 4, n - 2]]
    fakeobs *= 20
    obs = np.array(fakeobs, dtype=np.int32)
    stitchpoints = np.concatenate([[0.], np.cumsum(obs[:, 0])[::10], [obs[:, 0].sum()]]).astype(np.int32)
    im = smcpp._smcpp.PyOnePopInferenceManager(n - 2, [obs], hs, 0, stitchpoints, None)
    im.model = model
    im.theta = 0.0025000000000000001
    im.rho = 0.0025

    for _ in range(3):
        im.rho_vals[:] = 10 ** np.random.uniform(-2, -4, size=len(im.rho_vals))
        im.rho = im.rho  # this is needed to trigger the dirty flag and recompute transitions.
        im.E_step()
        print(im.rho_vals)
        print(im.loglik())

# Check match to the rjmcmc code
def test_likelihood():
    setup_logging(2, os.path.join(".", ".debug.txt"))

    #hs = np.concatenate([[0.0], np.logspace(-3, 1, 30), [np.inf]])
    hs = np.concatenate([[0.0], np.logspace(-3, 1, 3), [np.inf]])

    ## Initialize model
    demo = smcpp.util.human
    a = demo['a']
    b = demo['b']
    s = demo['s']
    model = smcpp.model.OldStyleModel(a,b,s)

    n = 30
    fakeobs = [[1, -1, 0, 0], [1, 1, 0, 0], [10, 0, 0, 0], [10, -1, 0, 0],
               [200000, 0, 0, n - 2], [1, 1, n - 4, n - 2]]
    # fakeobs *= 20
    fakeobs *= 3
    obs = np.array(fakeobs, dtype=np.int32)

    ## Init rho map -- vector of rho values at each SNP (length = # of snps)
    # obs = [obs[0][1:]]
    # ob = obs[0]
    # snp_blocks, snp_pos = find_SNPs(ob)
    # background = 0.0044 * ((snp_pos < 47500) | (snp_pos >= 52500))
    # hotspot = 0.044 * ((snp_pos >= 47500) & (snp_pos < 52500))
    # rho_map = background + hotspot

    stitchpoints = np.concatenate([[0.], np.cumsum(obs[:, 0])[::10], [obs[:, 0].sum()]]).astype(np.int32)
    im = smcpp._smcpp.PyOnePopInferenceManager(n - 2, [obs], hs, 0, stitchpoints, None)
    im.model = model
    im.theta = 0.0025000000000000001
    im.rho = 0.0025

    im.rho_vals[:] = [0.0044] * len(im.rho_vals)
    im.rho = im.rho  # this is needed to trigger the dirty flag and recompute transitions.
    im.E_step()
    print(im.rho_vals)
    print(im.loglik())


class StitchOptimizer():
    def __init__(self, stitchpoints, files):
        # Hidden states, number = 32
        hs = np.concatenate([[0.0], np.logspace(-3, 1, 30), [np.inf]])

        # Initialize model
        demo = smcpp.util.human
        a = demo['a']
        b = demo['b']
        s = demo['s']
        model = smcpp.model.OldStyleModel(a, b, s)

        #Init obs
        contigs = smcpp.estimation_tools.load_data(files)
        obs = []
        for c in contigs:
            obs.append(np.array(c.data,dtype=np.int32))

        # Initialize IM
        n = 50
        self._im = smcpp._smcpp.PyOnePopInferenceManager(n - 2, obs , hs, 0, stitchpoints, None)
        self._im.model = model
        self._im.theta = 0.0025
        self._im.rho = 0.0025
        self._im.rho_vals[:] = [0.0025]*(len(stitchpoints) - 1)
        #self._im.rho_vals[:] = [0.0025]


    def optimize(self, x0):
        #bnds = [(0.0,0.5)]
        bnds = [(0.001,0.2)]*len(x0)
        return scipy.optimize.minimize(self.helper, x0, bounds=bnds)

    def ll_surface(self, rhos, mode_left, mode_right):
        compute = lambda r: self.helper([mode_left, r, mode_right])
        ll = [-1.0*compute(r) for r in rhos]
        return ll

    def helper(self, rhovals):
        self._im.rho_vals[:] = rhovals
        self._im.rho = self._im.rho
        self._im.E_step()
        print(self._im.rho_vals)
        print(self._im.loglik())

        return -1.0*self._im.loglik()

def optimizer(stitchpoints, filename, num_composite = 25):
    files = []
    pi = np.random.permutation(25)
    for i in pi:
        file = "stitched_data/" + filename + ".smc." + str(i) + ".txt"
        print(file)
        files.append(file)
    so = StitchOptimizer(stitchpoints, files)
    x = so.optimize([0.008] * (len(stitchpoints)- 1))
    return x

# Tests flat human recomb with no stitchpoints
def test_no_hotspot(filename, num_composite):
    stitchpoints = [0,100000]
    optimal = optimizer(stitchpoints, filename , num_composite)
    return optimal

def test_hotspot(filename, hotspot_kb=2, num_composite=25):
    left = 50000 - hotspot_kb*1000/2
    right = 50000 + hotspot_kb*1000/2
    stitchpoints = [0,left,right,100000]
    print("Stitchpoints: ",stitchpoints)
    optimal = optimizer(stitchpoints, filename, num_composite)
    return optimal

# Tests flat human recomb with stitchpoints
# Width of middle stitchpoints
# num_composite - The number of samples for the composite likelihood
def test_flat(filename, out, num_composite=25):
    return print(filename, ": ", test_no_hotspot(filename,num_composite).x, file=out)

def test_tenx(length, out,  num_composite=25):
    print("Tenx_" + str(length) + "k:" ,test_hotspot("tenx_" + str(length) + "k", hotspot_kb=length).x, file=out)

def test_fiftyx(length, out, num_composite=25):
    print("Fiftyx_" + str(length) + "k:", test_hotspot("fiftyx_" + str(length) + "k", hotspot_kb=length).x, file=out)

def test_hundredx(length, out, num_composite=25):
    print("Hundredx_" + str(length) + "k:", test_hotspot("hundredx_" + str(length) + "k", hotspot_kb = length).x, file=out)


def save_likelihood_surface(surface_file, modes_file):
    rhos = np.linspace(0.0044, 0.1, 100)
    stitchpoints = [0,47500,52500,100000]
    modes = np.loadtxt(modes_file)
    replicates = len(modes[:,0])

    arr = np.zeros((25*replicates + 1, len(rhos)))
    #arr[0,:] = rhos ## STILL NEED TO SAVE STUFF! How do we want to store? Single file or multiple?
    for i in range(replicates):
        for j in range(25):
            files = ["stitched_data/tenx_5k_" + str(i) + ".smc." + str(j) + ".txt"]
            so = StitchOptimizer(stitchpoints, files)
            lls = so.ll_surface(rhos,modes[i,0],modes[i,2])
            arr[25*i+j,:] = lls

    arr[-1,:] = rhos
    np.savetxt(surface_file, arr)

def compute_modes(replicates):
    stitchpoints = [0,47500,52500,100000]
    print("Stitchpoints: ",stitchpoints)
    arr = np.zeros((replicates,3))
    for i in range(replicates):
        files = []
        for j in range(25):
            files.append("stitched_data/tenx_5k_" + str(i) + ".smc." + str(j) + ".txt")
        so = StitchOptimizer(stitchpoints, files)
        opt = so.optimize([0.008] * (len(stitchpoints)- 1))
        arr[i,:] = opt.x

    np.savetxt("tenx_5k_modes", arr)








if __name__ == "__main__":
    #save_likelihood_surface("tenx_5k_surfaces.txt", "tenx_5k_modes")
    test_likelihood()
    #compute_modes(20)
    # fname = "composite_ll_variance.txt"
    # if os.path.isfile(fname):
    #     raise RuntimeError("Don't overwrite the file!!")
    # else:
    #     out = open(fname, "wt")
    # rhos = np.linspace(0.0044, 0.1, 100)
    # stitchpoints = [0,47500,52500,100000]
    # #out = sys.stdout
    # #test_hotspot("flat_new",hotspot_kb=5)
    # # Tenx_5k: [ 0.00546596  0.04808182  0.00548466]


    # #print("Rhos: ",list(rhos), file=out)
    # replicates = 20
    # arr = np.zeros((replicates + 1, len(rhos)))
    # arr[0,:] = rhos
    # for i in range(replicates):
    #     print("Replicate: ", i, file=sys.stdout)
    #     # subsample
    #     pi = np.random.permutation(25)
    #     files = []
    #     for j in pi[:5]:
    #         files.append("stitched_data/tenx_5k.smc." + str(j) + ".txt")
    #     so = StitchOptimizer(stitchpoints, files)
    #     lls = so.ll_surface(rhos)
    #     arr[i+1,:] = lls

    # np.savetxt(fname, arr)

