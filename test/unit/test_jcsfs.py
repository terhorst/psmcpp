import pytest
import numpy as np

import smcpp
from smcpp import util
from smcpp.jcsfs import JointCSFS
from smcpp.model import SMCModel

@pytest.fixture
def model1():
    s = np.diff(np.logspace(np.log10(.001), np.log10(4), 20))
    ret = SMCModel(s, np.logspace(np.log10(.001), np.log10(4), 5))
    ret[:] = np.log(np.arange(1, 6)[::-1])
    return ret

@pytest.fixture
def model2():
    s = np.diff(np.logspace(np.log10(.001), np.log10(4), 20))
    ret = SMCModel(s, np.logspace(np.log10(.001), np.log10(4), 5))
    ret[:] = np.log(np.arange(1, 6))
    return ret

@pytest.fixture
def jcsfs():
    return JointCSFS(2, 2, 2, 0, [0.0, 0.5, 1.0, np.inf])

np.set_printoptions(precision=3, linewidth=100)

def test_marginal_pop1(model1, model2):
    ts = [0.0, 0.5, 1.0, np.inf]
    n1 = 10 
    n2 = 8
    j = JointCSFS(n1, n2, 2, 0, ts, 100)
    for split in [0.1, 0.25, 0.5, 0.75, 1.0, 2.0]:
        jc = j.compute(model1, model2, split)
        for t1, t2, jjc in zip(ts[:-1], ts[1:], jc):
            A1 = smcpp._smcpp.raw_sfs(model1, n1, t1, t2)
            A2 = jjc.sum(axis=(-1, -2))
            A2[0, 0] = A2[-1, -1] = 0.
            # assert np.allclose(A1, A2, 1e-2, 0)

def test_marginal_pop2(model1, model2):
    ts = [0.0, 0.5, 1.0, np.inf]
    n1 = 10 
    n2 = 8
    j = JointCSFS(n1, n2, 2, 0, ts, 100)
    for split in [0.1, 0.25, 0.5, 0.75, 1.0, 2.0]:
        jc = j.compute(model1, model2, split)
        for t1, t2, jjc in zip(ts[:-1], ts[1:], jc):
            A1 = util.undistinguished_sfs(smcpp._smcpp.raw_sfs(model2, n2 - 2, 0., np.inf))
            A2 = jjc.sum(axis=(0, 1, 2))
            assert np.allclose(A1, A2, 1e-1, 0)

def test_jcsfs(jcsfs, model1, model2):
    jcsfs.compute(model1, model2, 0.25)
