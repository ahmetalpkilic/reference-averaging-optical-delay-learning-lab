import unittest
import numpy as np
from lab import average_lags,nested_references,sequences,lag_batch,N

class AveragingTests(unittest.TestCase):
    def test_average_of_lags_not_lag_of_average(self):
        x=np.array([[[1.,0.],[2.,1.],[4.,3.]],[[4.,2.],[2.,3.],[1.,0.]]])
        expected=sum(lag_batch(t[None])[0] for t in x)/2
        np.testing.assert_allclose(average_lags(x),expected)
        self.assertGreater(np.linalg.norm(expected-lag_batch(x.mean(axis=0)[None])[0]),.1)
    def test_nested_prefix_and_single_reference(self):
        refs=nested_references(np.random.default_rng(71),1)
        rng=np.random.default_rng(71)
        x=np.concatenate([sequences(rng,1,N,0) for _ in range(16)])
        for ell in (1,4,16):
            np.testing.assert_allclose(refs[ell][0],average_lags(x[:ell]),atol=1e-18)
    def test_invalid_reference_banks(self):
        for x in [np.empty((0,96,2)),np.zeros((2,1,2)),np.full((2,96,2),np.nan)]:
            with self.assertRaises(ValueError):average_lags(x)
    def test_independent_matrix_covariance_one_over_l(self):
        refs=nested_references(np.random.default_rng(73),5000)
        covs={ell:np.cov(a[:,[0,0,1],[0,1,1]],rowvar=False) for ell,a in refs.items()}
        for ell in (4,16):
            self.assertLess(np.linalg.norm(ell*covs[ell]-covs[1])/np.linalg.norm(covs[1]),.10)
    def test_shared_reference_dependence_and_correlation_floor(self):
        rng=np.random.default_rng(74);count=150000;ell=16
        independent=rng.normal(size=(count,ell));common=rng.normal(size=count)
        mean=independent.mean(axis=1)
        u,v=rng.normal(size=(2,count))
        self.assertAlmostEqual(np.cov(u-mean,v-mean)[0,1],1/ell,delta=.012)
        correlated=np.sqrt(.3)*common[:,None]+np.sqrt(.7)*independent
        self.assertAlmostEqual(correlated.mean(axis=1).var(),.3+.7/ell,delta=.006)

if __name__=='__main__':unittest.main()
