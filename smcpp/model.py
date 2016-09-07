from __future__ import absolute_import, division, print_function
import numpy as np
import ad.admath

from . import spline, logging, util
from .observe import Observable, Observer, targets

logger = logging.getLogger(__name__)


# Dummy class used for JCSFS and a few other places
class PiecewiseModel(Observable):
    def __init__(self, a, s):
        assert len(a) == len(s)
        Observable.__init__(self)
        self.s = np.array(s)
        self.a = np.array(a)
        
    def stepwise_values(self):
        return self.a

    def regularizer(self):
        return (np.diff(self.a) ** 2).sum()

    def __getitem__(self, it):
        return self.a[it]

    def __setitem__(self, it, x):
        self.a[it] = x

    @property
    def dlist(self):
        ret = []
        for yy in self.a:
            try:
                ret += [d for d in yy.d() if d.tag is not None]
            except AttributeError:
                pass
        return list(set(ret))

class OldStyleModel(PiecewiseModel):
    def __init__(self, a, b, s):
        assert b[-1] == a[-1]
        ap = []
        sp = []
        for aa, bb, ss, cs in zip(a, b, s, util.cumsum0(s)[:-1]):
            print(aa,bb,ss,cs)
            if aa == bb:
                ap.append(aa)
                sp.append(ss)
            else:
                s0 = cs if cs > 0 else 1e-5
                s1 = s0 + ss
                # t = np.logspace(np.log10(s0), np.log10(s1), 40)
                t = np.linspace(s0, s1, 40)
                sp += np.diff(t).tolist()
                ap += (aa * (bb / aa) ** ((t[:-1] - s0) / (s1 - s0))).tolist()
        PiecewiseModel.__init__(self, ap, sp)


class SMCModel(Observable):
    def __init__(self, s, knots, spline_class=spline.PChipSpline):
        Observable.__init__(self)
        self._spline_class = spline_class
        self._s = np.array(s)
        self._cumsum_s = np.cumsum(s)
        self._knots = np.array(knots)
        self._trans = np.log
        # self._trans = lambda x: x
        self._spline = self._spline_class(self.transformed_knots)

    @property
    def s(self):
        return self._s

    @property
    def K(self):
        return len(self.knots)

    def randomize(self):
        self[:] += np.random.normal(0., .01, size=len(self[:]))

    @property
    def knots(self):
        return self._knots

    @property
    def transformed_knots(self):
        return self._trans(self._knots)

    def __setitem__(self, key, item):
        self._spline[key] = item
        self.update_observers('model update')

    def __getitem__(self, key):
        return self._spline[key]

    @property
    def dlist(self):
        ret = []
        for yy in self[:]:
            try:
                ret += [d for d in yy.d() if d.tag is not None]
            except AttributeError:
                pass
        return list(set(ret))

    def regularizer(self):
        ret = self._spline.roughness()
        ret += (np.diff(self[:]) ** 2).sum()
        if not isinstance(ret, ad.ADF):
            ret = ad.adnumber(ret)
        return ret

    def __call__(self, x):
        'Evaluate :self: at points x.'
        ret = np.array(
            ad.admath.exp(self._spline(self._trans(x)))
        )
        return ret

    def stepwise_values(self):
        return self(np.cumsum(self._s))

    def reset(self):
        self[:] = 0.

    def to_s(self, until=None):
        ret = []
        for ary in [self[:until], self.stepwise_values()]:
            ary = ary.astype('float')
            fmt = " ".join(["{:>5.2f}"] * len(ary))
            ret.append(fmt.format(*ary))
        return "\n" + "\n".join(ret)

    def to_dict(self):
        return {
                'class': self.__class__.__name__,
                's': list(self._s),
                'knots': list(self._knots),
                'spline_class': self._spline_class.__name__,
                'y': self[:].astype('float').tolist()
                }

    @classmethod
    def from_dict(cls, d):
        assert cls.__name__ == d['class']
        spc = getattr(spline, d['spline_class'])
        r = cls(d['s'], d['knots'], spc)
        r[:] = d['y']
        return r

    def distinguished_model(self, index):
        assert index == 0
        return self

    def copy(self):
        return SMCModel.from_dict(self.to_dict())


class SMCTwoPopulationModel(Observable, Observer):

    def __init__(self, model1, model2, split):
        Observable.__init__(self)
        self._models = [model1, model2]
        model1.register(self)
        model2.register(self)
        self._split = split

    # Propagate changes from submodels up
    @targets('model update')
    def update(self, message, *args, **kwargs):
        self.update_observers('model update')

    @property
    def split(self):
        return self._split

    @split.setter
    def split(self, x):
        self._split = x
        self.update_observers('model update')

    @property
    def split_ind(self):
        'Return k such that model2.t[k] <= split < model2.t[k + 1]'
        return np.searchsorted(self.model2.knots, self._split, side="right") - 1

    @property
    def s(self):
        return self.model1.s

    def splitted_models(self):
        # return [self.model1, ConcatenatedModel(self.model1, self.model2, self.split)]
        return [self.model1, _concat_models(self.model1, self.model2, self.split)]

    @property
    def model1(self):
        return self._models[0]

    @property
    def model2(self):
        return self._models[1]

    def distinguished_model(self, index=0):
        if index is not None:
            return self.splitted_models()[index]
        # The "apart" case
        s = self.model1.s
        a = self.model1.stepwise_values()
        cs = util.cumsum0(self.model1.s)
        cs[-1] = np.inf
        ip = np.searchsorted(cs, self._split)
        sp = np.diff(np.insert(cs, ip, self._split))
        sp[-1] = 1.
        s = sp[ip - 1:]
        s[0] = self.split
        a = np.insert(a[ip - 1:], 0, np.inf)
        return PiecewiseModel(a, s)

    @property
    def dlist(self):
        return list(set(self._models[0].dlist + self._models[1].dlist))

    def randomize(self):
        for m in self._models:
            m.randomize()

    def reset(self):
        for m in self._models:
            m.reset()

    def to_dict(self):
        return {'class': self.__class__.__name__,
                'model1': self._models[0].to_dict(),
                'model2': self._models[1].to_dict(),
                'split': float(self._split)}

    @classmethod
    def from_dict(cls, d):
        assert cls.__name__ == d['class']
        model1 = SMCModel.from_dict(d['model1'])
        model2 = SMCModel.from_dict(d['model2'])
        return cls(model1, model2, d['split'])

    def to_s(self):
        return "\nPop. 1:\n{}\nPop. 2:\n{}\nSplit: {:.3f}".format(
            self._models[0].to_s(), self._models[1].to_s(self.split_ind),
            self.split)

    # FIXME this counts the part before the split twice
    def regularizer(self):
        ret = sum([x.regularizer() for x in self.splitted_models()])
        if not isinstance(ret, ad.ADF):
            ret = ad.adnumber(ret)
        return ret

    def __getitem__(self, coords):
        a, cc = coords
        return self._models[a][cc]

    def __setitem__(self, coords, x):
        a, cc = coords
        # This will generate 'model updated' messages in the submodels.
        self._models[a][cc] = x


def _concat_models(m1, m2, t):
    # ip = np.searchsorted(m1._knots, t, side="right")
    ip = np.argmin(np.abs(m1._knots - t))
    nk = m1._knots.copy()
    nk[ip] = t
    ny = np.zeros_like(m2[:])
    ny[:ip] = m2[:ip]
    ny[ip] = ad.admath.log(m1(t).item())
    ny[ip + 1:] = m1[ip + 1:]
    ret = SMCModel(m1.s, nk, m1._spline_class)
    ret[:] = ny
    return ret
