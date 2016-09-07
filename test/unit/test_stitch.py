import numpy as np
import smcpp._smcpp, smcpp.model, smcpp.spline

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
