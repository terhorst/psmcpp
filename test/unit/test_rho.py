import numpy as np
import smcpp._smcpp, smcpp.model, smcpp.spline, smcpp.util, smcpp.estimation_tools
import scipy.optimize


class RhoOptimizer():
    def __init__(self, files):
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
        self._im = smcpp._smcpp.PyOnePopInferenceManager(n - 2, obs , hs, 0, None)
        self._im.model = model
        self._im.theta = 0.0025
        self._im.rho = 90.0
        #self._im.rho_vals[:] = [0.001, 0.01, 0.001]


    def optimize(self, x0):
        #bnds = ((0.00001, 0.1), (0.0001, 0.1), (0.00001, 0.01))
        bnds = [(0.0, 0.5)]
        return scipy.optimize.minimize(self.helper, x0, bounds=bnds)

    def helper(self, rho):
        self._im.rho = rho[0]
        self._im.E_step()
        print(self._im.rho)
        print(self._im.loglik())

        return -1.0*self._im.loglik()



if __name__ == "__main__":
    #stitchpoints = [0,49000,51000,100000]
    files = []
    for i in range(25):
        file = "stitched_data/flat.smc." + str(i) + ".txt"
        files.append(file)
    ro = RhoOptimizer(files)
    x = ro.optimize([0.1])