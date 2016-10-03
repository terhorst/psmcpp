import numpy as np
import smcpp._smcpp, smcpp.model, smcpp.spline, smcpp.util, smcpp.estimation_tools
import scipy.optimize

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

        print obs[0]

        # Initialize IM
        n = 50
        self._im = smcpp._smcpp.PyOnePopInferenceManager(n - 2, obs , hs, 0, stitchpoints, None)
        self._im.model = model
        self._im.theta = 0.0025
        self._im.rho = 0.0025
        self._im.rho_vals[:] = [0.001, 0.01, 0.001]


    def optimize(self, x0):
        bnds = ((0.00001, 0.1), (0.0001, 0.1), (0.00001, 0.01))
        return scipy.optimize.minimize(self.helper, x0, bounds=bnds)

    def helper(self, rhovals):
        self._im.rho_vals[:] = rhovals
        self._im.rho = self._im.rho
        self._im.E_step()
        print(self._im.rho_vals)
        print(self._im.loglik())

        return -1.0*self._im.loglik()





if __name__ == "__main__":
    stitchpoints = [0,49000,51000,100000]
    files = []
    for i in range(5):
        file = "stitched_data/hundredx.smc." + str(i) + ".txt"
        files.append(file)
    so = StitchOptimizer(stitchpoints, files)
    x = so.optimize([0.0025,0.025,0.0025])
