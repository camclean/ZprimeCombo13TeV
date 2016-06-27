# -*- coding: utf-8 -*-
from theta_auto.test_model import *

clean_workdir()

import unittest
import time

config.suppress_info = True
one_sigma = 0.6827

class TestMle(unittest.TestCase):
    def setUp(self):
        self.model_nobkg = simple_counting(s = 100, n_obs = 100)
        self.model_s2 = simple_counting(s = 1.0, n_obs = 11.0, b = 8.0, s2 = 1.5)
        self.model_bunc = simple_counting(s = 1.0, n_obs = 11.0, b = 8.0, b_uncertainty = 2.0)

    def test_mle(self):
        res = mle(self.model_s2, 'data', 1)
        self.assertEqual(len(res['s']['beta_signal']), 1)
        self.assertEqual(len(res['s2']['beta_signal']), 1)
        self.assertAlmostEqual(res['s']['beta_signal'][0][0], 3.0, places = 3)
        self.assertAlmostEqual(res['s2']['beta_signal'][0][0], 2.0, places = 3)
        res = mle(self.model_bunc, 'data', 1)
        self.assertEqual(len(res), 1)
        self.assertAlmostEqual(res['s']['bunc'][0][0], 0.0, places = 3)
        self.assertAlmostEqual(res['s']['beta_signal'][0][0], 3.0, places = 3)

    def test_pl(self):
        res = pl_interval(self.model_nobkg, 'data', n=1, cls = [one_sigma])
        # res is dict (sp) --> (cl) --> list of length n with (value, error) tuples
        self.assertEqual(len(res), 1)
        self.assertEqual(len(res['s'][one_sigma]), 1)
        # the interval should be about +- 10%:
        self.assertAlmostEqual(res['s'][one_sigma][0][0], 0.90, places = 2)
        self.assertAlmostEqual(res['s'][one_sigma][0][1], 1.10, places = 2)
        
    def test_chi2(self):
        res = mle(self.model_nobkg, 'toys:1.0', 1000, chi2 = True, signal_prior = 'fix:1.0')
        mean, width = get_mean_width(res['s']['__chi2'])
        self.assertTrue(abs(mean - 1.0) < 0.1)
        self.assertTrue(abs(width - math.sqrt(2)) < 0.2)
        m2 = multichannel_counting([10000., 10000.])
        res = mle(m2, 'toys:1.0', 3000, chi2 = True, signal_prior = 'fix:1.0')
        mean, width = get_mean_width(res['s']['__chi2'])
        self.assertTrue(abs(mean - 2.0) < 0.1)
        self.assertTrue(abs(width - math.sqrt(4)) < 0.2)
        
        # test absolute chi2 values:
        s = [10000., 10000.]
        n_obs = [10170., 9786.]
        model = multichannel_counting(s, n_obs = n_obs)
        res = mle(model, 'data', 1, chi2 = True, signal_prior = 'fix:1.0', nuisance_constraint = get_fixed_dist(model.distribution))
        expected_chi2 = sum(map(lambda (mu, n): (mu - n)**2 / mu, zip(s, n_obs)))
        self.assertAlmostEqual(res['s']['__chi2'][0], expected_chi2, places = 1)
        
        
    def test_signal_prior(self):
        res = mle(self.model_nobkg, 'toys:1.0', 10, signal_prior = 'flat:[0,1]')
        bs = [x[0] for x in res['s']['beta_signal']]
        self.assertTrue(all(map(lambda x: x<=1.0 and x>=0.0, bs)))
        
        res = mle(self.model_nobkg, 'toys:1.0', 10, signal_prior = 'fix:0.7')
        bs = [x[0] for x in res['s']['beta_signal']]
        self.assertTrue(all(map(lambda x: x<=0.7001 and x>=0.6999, bs)))
        
    def test_cov(self):
        res = mle(self.model_bunc, 'toys-asimov:1.0', 1, with_covariance = True)
        self.assertTrue('s' in res)
        res = res['s']
        self.assertTrue('__cov' in res)
        cov = res['__cov'][0]
        self.assertEqual(len(cov), 2)
        self.assertEqual(len(cov[0]), 2)
        self.assertEqual(len(cov[1]), 2)
        self.assertEqual(cov[0][1], cov[1][0])
        """
        pars = self.model_bunc.get_parameters(['s'])
        print cov
        print res
        for i, p in enumerate(pars):
            error_from_cov = math.sqrt(cov[i][i])
            error_from_fit = res[p][0][1]
            self.assertAlmostEqual(error_from_cov, error_from_fit)
        """
        
    def test_pl_termonly(self):
        model_termonly = Model()
        model_termonly.additional_nll_term = NLGauss(['beta_signal', 'b'], [1.0, 0.0], [(0.2**2, 0.0), (0.0, 0.1**2)])
        model_termonly.distribution.set_distribution('b', 'gauss', 0.0, inf, [-inf, inf])
        res = pl_interval(model_termonly, 'toys:0', n=1, cls = [one_sigma], signal_process_groups = {'': []})
        self.assertAlmostEqual(res[''][one_sigma][0][0], 0.80, places = 4)
        self.assertAlmostEqual(res[''][one_sigma][0][1], 1.20, places = 4)


class TestCls(unittest.TestCase):
    def setUp(self):
        self.model = simple_counting(s = 100, n_obs = 10100, b = 10000)
    
    def test_cls(self):
        exp, obs = cls_limits(self.model)
        
    def test_cls_freq(self):
        exp, obs = cls_limits(self.model, frequentist_bootstrapping = True)

    def test_asymptotic_cls(self):
        exp, obs = asymptotic_cls_limits(self.model)
        print exp, obs

class TestBB(unittest.TestCase):
    # test the Barlow-Beeston treatment of MC uncertainties: using a simple counting, no-background experiment, the relative uncertainties
    # should approximately add quadratically. Use very large signal to make sure we have an approx. Gauss case.
    def testbb(self):
        n = 10000.0
        for unc in (0.001, 0.01): # a relative uncertainty
            model_bb = simple_counting_bb(s = n, n_obs = n, s_uncertainty = unc * n)
            model_bb0 = simple_counting_bb(s = n, n_obs = n, s_uncertainty = 0.0)
            res = pl_interval(model_bb0, 'data', n = 1, cls = [one_sigma])
            lower, upper = res['s'][one_sigma][0]
            self.assertAlmostEqual(lower, 1 - math.sqrt(1.0 / n), places = 3)
            self.assertAlmostEqual(upper, 1 + math.sqrt(1.0 / n), places = 3)
    
            res = pl_interval(model_bb, 'data', n = 1, cls = [one_sigma])
            lower, upper = res['s'][one_sigma][0]
            self.assertAlmostEqual(lower, 1 - math.sqrt(1.0 / n + unc**2), places = 3)
            self.assertAlmostEqual(upper, 1 + math.sqrt(1.0 / n + unc**2), places = 3)
    
    
    def test_bb_twochannel(self):
        model_bb = template_counting_bb(s = [100.0], n_obs = [101.0], s_uncertainty = [10.0])
        res = pl_interval(model_bb, 'data', n = 1, cls = [one_sigma])
        oc_lower, oc_upper = res['s'][one_sigma][0]
        #print "oc", oc_lower, oc_upper
        
        model_bb = template_counting_bb(s = [100.0, 0.0], n_obs = [101.0, 2.0], s_uncertainty = [10.0, 1.0])
        res = pl_interval(model_bb, 'data', n = 1, cls = [one_sigma])
        lower, upper = res['s'][one_sigma][0]
        #print "tc", lower, upper
        
        # interval should have shifted upwards, but not change size much:
        self.assertTrue(upper > oc_upper)
        self.assertTrue(lower > oc_lower)
        relsize = (upper - lower) / (oc_upper - oc_lower)
        self.assertTrue(relsize > 0.97 and relsize < 1.03)

        
class TestBayes(unittest.TestCase):
    def setUp(self):
        self.model = simple_counting(s = 1.0, n_obs = 11.0, b = 8.0, b_uncertainty = 2.0)

    def test_bayesian_quantiles(self):
        t0 = time.time()
        res0 = bayesian_quantiles(self.model, 'data', 100, quantiles = [0.95])
        time0 = time.time() - t0
        quants0 = sorted(res0['s'][0.95])
        t0 = time.time()
        options = Options()
        options.set('main', 'n_threads', '2')
        res1 = bayesian_quantiles(self.model, 'data', 100, quantiles = [0.95], options = options)
        time1 = time.time() - t0
        quants1 = sorted(res1['s'][0.95])
        #print "real time elapsed: ", time0, time1
        #print "expected limits: ", quants0[len(quants0) / 2], quants1[len(quants1) / 2]
    
    def test_bayesian_quantiles2(self):
        model = simple_counting(s = 10000., n_obs = 11.0, b = 8.0, b_uncertainty = 2.0)
        q0 = bayesian_quantiles(model, 'data', 100)
        q0 = sorted(q0['s'][0.95])
        q1 = bayesian_quantiles(model, 'data', 100, hint_method = None)
        q1 = sorted(q1['s'][0.95])
        #print "without hint:", q1
        #print "with hint: ", q0
    
        
    def test_bayesian_posterior_model_prediction(self):
        res = bayesian_posterior_model_prediction(self.model, 'data', 10)
        self.assertTrue(len(res)==1)
        self.assertTrue(len(res['s']['obs'])==10)
        values = res['s']['obs'][0].get_values()
        uncs = res['s']['obs'][0].get_uncertainties()
        self.assertTrue(len(values)==1)
        self.assertTrue(len(uncs)==1)
        self.assertTrue(values[0] < 15.0 and values[0] > 10.0)
        self.assertTrue(uncs[0] < 5.0 and uncs[0] > 0.0)
        
class TestSqlite(unittest.TestCase):
    
    def test_many_columns(self):
        db = DbResult('../../test_many_cols.db') # that one should have been created by bin/test
        res = db.get_results('test_table', '*')
        for i in range(5000):
            colname = 'col%d' % i
            self.assertTrue(colname in res)
        res = db.get_results('test_table', ['col0', 'col1'])
        self.assertEqual(len(res), 2)
        self.assertTrue('col0' in res and 'col1' in res)
        res = db.get_results('test_table', ['col1000', 'col2000'])
        self.assertEqual(len(res), 2)
        self.assertTrue('col1000' in res and 'col2000' in res)
        
        
        
class TestModel(unittest.TestCase):
    def test_rebin(self):
        model = test_model.gaussoverflat(100, 1000)
        mle(model, 'toys-asimov:1.0', 1)
        model.rebin('obs', 2)
        mle(model, 'toys-asimov:1.0', 1)

class TestGetData(unittest.TestCase):

    def setUp(self):
        self.model = test_model.multichannel_counting([0., 0.], backgrounds = [1000., 2000.], b_uncertainty1 = [0.1, 0.1])

    def test_gd(self):
        res = make_data(self.model, 'toys:0.0', 1000, retval = 'data')
        # get info about total yields:
        mean0, unc0 = get_mean_width([h.get_value_sum() for h in res['s']['obs0']])
        mean1, unc1 = get_mean_width([h.get_value_sum() for h in res['s']['obs1']])
        self.assertTrue((mean0 - 1000.) < 20.)
        self.assertTrue((mean1 - 2000.) < 40.)
        # the uncertainty should be about sqrt(1000) + 10% (quadratically added) for
        # obs0 and sqrt(2000) + 10% for obs1:
        unc0_expected = math.sqrt(1000. + (100.)**2)
        unc1_expected = math.sqrt(2000. + (200.)**2)
        self.assertLess(abs(1 - unc0 / unc0_expected), 0.1)
        self.assertLess(abs(1 - unc1 / unc1_expected), 0.1)

        # this is how one would get the per-bin uncertainties from the histogrm list:
        hist0 = get_mean_width_hist(res['s']['obs0'])
        hist1 = get_mean_width_hist(res['s']['obs1'])
        self.assertAlmostEqual(hist0.get_values()[0], mean0, 2)
        self.assertAlmostEqual(hist0.get_uncertainties()[0], unc0, 2)
        self.assertAlmostEqual(hist1.get_values()[0], mean1, 2)
        self.assertAlmostEqual(hist1.get_uncertainties()[0], unc1, 2)

    def test_no_poisson(self):
        # now the prior uncertainty *without* statistical uncertainties, but including systematic ones.
        # Note that 'toys-asimov' as default does not vary systematics, so we have to switch that on again by explicitly setting nuisance_prior_toys ...
        res = make_data(self.model, input = 'toys-asimov:0.0', n = 1000, nuisance_prior_toys = self.model.distribution, retval = 'data')
        mean0, unc0 = get_mean_width([h.get_value_sum() for h in res['s']['obs0']])
        mean1, unc1 = get_mean_width([h.get_value_sum() for h in res['s']['obs1']])
        unc0_expected = 100.
        unc1_expected = 200.
        self.assertLess(abs(1 - unc0 / unc0_expected), 0.1)
        self.assertLess(abs(1 - unc1 / unc1_expected), 0.1)


    def test_asymmetric_uncertainty(self):
        # test the procedure for a very asymmetric uncertainty on the background -- like here: -1% / +10% systematic uncertainty
        b0 = 10000.
        m = test_model.simple_counting(s = 0.0, b = b0)
        m.add_asymmetric_lognormal_uncertainty('b_unc', 0.01, 0.1, 'b')
        # only use systematic uncertainties, see no_poisson:
        res = make_data(m, input = 'toys-asimov:0.0', n = 1000, nuisance_prior_toys = m.distribution, retval = 'data')
        # use the get_asymmetric_errors routine to calculate the uncertainty on the background.
        # note that we plugin in the "best estimate" for the background here (instead of using the mean
        # value in the toys) in order to compute the errors:
        mean, error_minus, error_plus = get_asymmetric_errors([h.get_value_sum() for h in res['s']['obs']], b0)
        print error_minus, error_plus
        # error_minus should be about 100, error_plus about 1000:
        error_plus_expected = b0 * (math.exp(0.1) - 1.0)
        error_minus_expected = b0 * (math.exp(0.01) - 1.0)
        self.assertLess(abs(error_minus / error_minus_expected - 1.), 0.1)
        self.assertLess(abs(error_plus / error_plus_expected - 1.), 0.15)
        

class TestRootModel(unittest.TestCase):
    
    @staticmethod
    def th1(name, xmin, xmax, data, uncertainties, root_dir):
        res = ROOT.TH1D(name, name, len(data), xmin, xmax)
        res.SetDirectory(root_dir)
        for i in range(len(data)):
            res.SetBinContent(i+1, data[i])
            res.SetBinError(i+1, uncertainties[i])
        return res
        
    
    def setUp(self):
        self.fn = 'test_root_model.root'
        f = ROOT.TFile(self.fn, "recreate")
        f.cd()
        hs = []
        # build a model with 2 channels (=two one-bin histos): "chan1" and "chan2" with one signal (100.0, 10.0),
        # one background ([100.0, 3.0]), and one background systematic uncertainty "unc" (+-10.0, +-0.3).
        # The MC uncertainties are all set to 10%. Data is (200.0, 13.0)
        h = TestRootModel.th1('chan1__s', 0.0, 1.0, [100.0], [10.0], f)
        hs.append(h)
        h = TestRootModel.th1('chan1__b', 0.0, 1.0, [100.0], [10.0], f)
        hs.append(h)
        h = TestRootModel.th1('chan2__s', 0.0, 1.0, [10.0], [1.0], f)
        hs.append(h)
        h = TestRootModel.th1('chan2__b', 0.0, 1.0, [3.0], [0.3], f)
        hs.append(h)
        h = TestRootModel.th1('chan2__b__unc__plus', 0.0, 1.0, [110.0], [0.0], f)
        hs.append(h)
        h = TestRootModel.th1('chan2__b__unc__minus', 0.0, 1.0, [90.0], [0.0], f)
        hs.append(h)
        h = TestRootModel.th1('chan1__DATA', 0.0, 1.0, [200.0], [0.0], f)
        hs.append(h)
        h = TestRootModel.th1('chan2__DATA', 0.0, 1.0, [13.0], [0.0], f)
        hs.append(h)
        f.Write()
        f.Close()
        del f
        
    
    def tearDown(self):
        try: os.unlink(self.fn)
        except RuntimeError: pass
    
    
    def test_basic(self):
        model = build_model_from_rootfile(self.fn)
        model.set_signal_processes('s')
        self.assertEqual(model.processes, set(['s', 'b']))
        self.assertEqual(set(model.get_parameters(['s'])), set(['unc', 'beta_signal']))
        self.assertEqual(set(model.observables.keys()), set(['chan1', 'chan2']))
        
        # filter out background uncertainty this time:
        model = build_model_from_rootfile(self.fn, histogram_filter = (lambda s: s.count('__')==1))
        model.set_signal_processes('s')
        self.assertEqual(model.processes, set(['s', 'b']))
        self.assertEqual(set(model.get_parameters(['s'])), set(['beta_signal']))
        self.assertEqual(set(model.observables.keys()), set(['chan1', 'chan2']))
        
    
    def test_mc_uncertainties(self):
        model = build_model_from_rootfile(self.fn)
        model.set_signal_processes('s')
        res = pl_interval(model, 'data', n = 1, cls = [one_sigma])
        nlower, nupper = res['s'][one_sigma][0]
        
        model = build_model_from_rootfile(self.fn, include_mc_uncertainties = True)
        model.set_signal_processes('s')
        res = pl_interval(model, 'data', n = 1, cls = [one_sigma])
        lower, upper = res['s'][one_sigma][0]
        
        #MC stat uncertainty should enlarge interval:
        self.assertTrue(lower < nlower and upper > nupper)
        
        
class MCMCHighdimtest(unittest.TestCase):        
    def test_100d(self):
        model = Model()
        ndim = 200
        parameters = ['p%d' % i for i in range(ndim)]
        mu = [0.0] * ndim
        cov = []
        for i in range(ndim):
            cov.append([1.0 if i==j else 0.0 for j in range(ndim)])
            model.distribution.set_distribution('p%d' % i, 'gauss', mean = 0.0, width = inf, range = (-inf, inf))
        model.additional_nll_term = NLGauss(parameters, mu, cov)
        options = Options()
        options.set('mcmc', 'diag', 'True')
        options.set('mcmc', 'strategy', 'asimov_widths')
        res = bayesian_quantiles(model, 'toys:0.0', 10, quantiles = [0.5], signal_process_groups = {'': []}, parameter = 'p0', iterations = 10000, options = options)
        m, w = get_mean_width(res[''][0.5])
        # less than 4sigma effect:
        self.assertTrue(abs(m) / w * math.sqrt(len(res[''])) < 4.0)
        options.set('mcmc', 'strategy', 'asimov_widths_1d')
        res = bayesian_quantiles(model, 'toys:0.0', 10, quantiles = [0.5], signal_process_groups = {'': []}, parameter = 'p0', iterations = 10000, options = options)
        m, w = get_mean_width(res[''][0.5])
        # less than 4sigma effect:
        self.assertTrue(abs(m) / w * math.sqrt(len(res[''])) < 4.0)
        options.set('mcmc', 'strategy', 'asimov_der_cov')
        res = bayesian_quantiles(model, 'toys:0.0', 10, quantiles = [0.5], signal_process_groups = {'': []}, parameter = 'p0', iterations = 10000, options = options)
        m, w = get_mean_width(res[''][0.5])
        # less than 4sigma effect:
        self.assertTrue(abs(m) / w * math.sqrt(len(res[''])) < 4.0)
        
mletests = unittest.TestLoader().loadTestsFromTestCase(TestMle)
suite2 = unittest.TestLoader().loadTestsFromTestCase(TestBB)
suite3 = unittest.TestLoader().loadTestsFromTestCase(TestRootModel)
bayes = unittest.TestLoader().loadTestsFromTestCase(TestBayes)
mcmc = unittest.TestLoader().loadTestsFromTestCase(MCMCHighdimtest)
cls = unittest.TestLoader().loadTestsFromTestCase(TestCls)
sqlite = unittest.TestLoader().loadTestsFromTestCase(TestSqlite)
model = unittest.TestLoader().loadTestsFromTestCase(TestModel)
gd = unittest.TestLoader().loadTestsFromTestCase(TestGetData)
#alltests = unittest.TestSuite([gd])
alltests = unittest.TestSuite([mletests, suite2, suite3, bayes, cls, sqlite, model, gd])

# verbose version:
res = unittest.TextTestRunner(verbosity=2).run(alltests)

# silent version:
#f = open('/dev/null', 'a')
#res = unittest.TextTestRunner(stream = f, descriptions = False, verbosity=0).run(alltests)

print "Failures=%d" % len(res.failures + res.errors)

