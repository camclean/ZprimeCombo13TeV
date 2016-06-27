# -*- coding: utf-8 -*-
# this file defines some simple models useful for testing theta-auto

from Model import *

# returns a model with a counting experiment in one bin with the given backgruond with a log-normal uncertainty.
# b_uncertainty is the absolute uncertainty on b.
#
# If s2 is not None, will return a model with two signal processes, "s" and "s2"
def simple_counting(s, n_obs=None, b=0.0, b_uncertainty=0.0, s2 = None):
    model = Model()
    if n_obs is not None:
        model.set_data_histogram('obs', Histogram(0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, [float(s)]))
    model.set_histogram_function('obs', 's', hf_s)
    if s2 is not None:
        hf_s2 = HistogramFunction()
        hf_s2.set_nominal_histo(Histogram(0.0, 1.0, [float(s2)]))
        model.set_histogram_function('obs', 's2', hf_s2)
    model.set_signal_processes('s*')
    if b > 0:
        hf_b = HistogramFunction()
        hf_b.set_nominal_histo(Histogram(0.0, 1.0, [float(b)]))
        model.set_histogram_function('obs', 'b', hf_b)
        if b_uncertainty > 0: model.add_lognormal_uncertainty('bunc', 1.0 * b_uncertainty / b, 'b')
    return model
    
# 2 signals in 2 bins (one in each), modeled as one observable
def counting_2signals(s1= 3.33, s2 = 10., nobs = [110.0, 120.0], b = [100.0, 100.0]):
    model = Model()
    model.set_data_histogram('obs', Histogram(0.0, 1.0, nobs))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, [float(s1), 0.0]))
    model.set_histogram_function('obs', 's1', hf_s)
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, [0.0, float(s2)]))
    model.set_histogram_function('obs', 's2', hf_s)
    model.set_signal_process_groups({'twos': {'beta1': ['s1'], 'beta2': ['s2']}})
    hf_b = HistogramFunction()
    hf_b.set_nominal_histo(Histogram(0.0, 1.0, b))
    model.set_histogram_function('obs', 'b', hf_b)
    return model
    
# signals are the signal yields, backgrounds are the background yields
# b_uncertainty1, b_uncertainty2 and b_uncertainty3 are either None or an array
# of *relative* background uncertainties in the channels.
def multichannel_counting(signals, n_obs = None, backgrounds = None, b_uncertainty1 = None, b_uncertainty2 = None, b_uncertainty3 = None, obsnames = None):
    n = len(signals)
    assert n_obs is None or n==len(n_obs)
    assert backgrounds is None or n==len(backgrounds)
    model = Model()
    for i in range(n):
        if obsnames is None: obsname = 'obs%d' % i
        else: obsname = obsnames[i]
        if n_obs is not None: model.set_data_histogram(obsname, Histogram(0.0, 1.0, [float(n_obs[i])]))
        hf_s = HistogramFunction()
        hf_s.set_nominal_histo(Histogram(0.0, 1.0, [float(signals[i])]))
        model.set_histogram_function(obsname, 's', hf_s)
        
        if backgrounds is None: continue
        hf_b = HistogramFunction()
        hf_b.set_nominal_histo(Histogram(0.0, 1.0, [float(backgrounds[i])]))
        model.set_histogram_function(obsname, 'b', hf_b)
        for pname, b_uncertainty in zip(('bunc1', 'bunc2', 'bunc3'), (b_uncertainty1, b_uncertainty2, b_uncertainty3)):
            if b_uncertainty is None: continue
            assert len(b_uncertainty) == n
            unc = float(b_uncertainty[i])
            if unc != 0.0:
                if unc == float("inf"):
                    model.add_lognormal_uncertainty(pname, 1.0, 'b', obsname)
                    model.distribution.set_distribution_parameters(pname, width = float("inf"))
                else:
                    model.add_lognormal_uncertainty(pname, unc, 'b', obsname)
    model.set_signal_processes('s*')
    return model


# returns a model in one bin with the given background. b_plus and b_minus are the background yields at +1sigma and -1sigma
def simple_counting_shape(s, n_obs = None, b=0.0, b_plus=0.0, b_minus=0.0):
    model = Model()
    if n_obs is not None:
        model.set_data_histogram('obs', Histogram(0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, [float(s)]))
    model.set_histogram_function('obs', 's', hf_s)
    model.set_signal_processes('s*')
    if b > 0:
        hf_b = HistogramFunction()
        hf_b.set_nominal_histo(Histogram(0.0, 1.0, [float(b)]))
        plus_histo = Histogram(0.0, 1.0, [float(b_plus)])
        minus_histo = Histogram(0.0, 1.0, [float(b_minus)])
        hf_b.set_syst_histos('bunc', plus_histo, minus_histo)
        model.distribution.set_distribution('bunc', 'gauss', mean = 0.0, width = 1.0, range = [-float("inf"), float("inf")])
        model.set_histogram_function('obs', 'b', hf_b)
    return model
    
# simple counting model with s signal events, no background, but an BB uncertainty on the signal of s_uncertainty (an absolute uncertainty)
def simple_counting_bb(s, s_uncertainty, n_obs):
    model = Model()
    model.set_data_histogram('obs', Histogram(0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, [float(s)], [float(s_uncertainty)]))
    model.set_histogram_function('obs', 's', hf_s)
    model.set_signal_processes('s*')
    model.bb_uncertainties = True
    return model

# counting model with s signal events, background, and b_uncertainty background uncertainty treated with BB light (an absolute uncertainty)
def counting_bb(s, b, b_uncertainty, n_obs):
    model = Model()
    model.set_data_histogram('obs', Histogram(0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, [float(s)]))
    hf_b = HistogramFunction()
    hf_b.set_nominal_histo(Histogram(0.0, 1.0, [float(b)], [float(b_uncertainty)]))
    model.set_histogram_function('obs', 's', hf_s)
    model.set_histogram_function('obs', 'b', hf_b)
    model.set_signal_processes('s*')
    model.bb_uncertainties = True
    return model
    
# in this case, s, s_uncertainty and n_obs should all be lists of the same length of floating point values.
def template_counting_bb(s, s_uncertainty, n_obs):
    assert len(s) == len(s_uncertainty) and len(s) == len(n_obs)
    model = Model()
    model.set_data_histogram('obs', Histogram(0.0, 1.0, n_obs))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(Histogram(0.0, 1.0, s, s_uncertainty))
    model.set_histogram_function('obs', 's', hf_s)
    model.set_signal_processes('s*')
    model.bb_uncertainties = True
    return model


# a gaussian "signal" distribution (mean s_mean, width s_sigma) over flat background on a range 0--100, with 100 bins; no observed data.
# b_uncertainty is the (absolute!) uncertainty on the background yield, which will be handeled with a log-normal.
def gaussoverflat(s, b, b_uncertainty = 0.0, s_mean = 50.0, s_sigma = 20., obs_suffix = ''):
    model = Model()    
    # something proportional to a Gauss:
    gauss = lambda x, mean, sigma: math.exp(-(x - mean)**2 / (2.0*sigma))
    hf_s = HistogramFunction()
    data = [gauss(x+0.5, s_mean, s_sigma) for x in range(100)]
    # normalize data to s:
    n_data = sum(data)
    data = [x*s/n_data for x in data]
    hf_s.set_nominal_histo(Histogram(0.0, 100.0, data))
    model.set_histogram_function('obs%s' % obs_suffix, 's', hf_s)
    
    hf_b = HistogramFunction()
    data = [b * 0.01] * 100
    hf_b.set_nominal_histo(Histogram(0.0, 100.0, data))
    model.set_histogram_function('obs%s' % obs_suffix, 'b', hf_b)
    if b_uncertainty > 0: model.add_lognormal_uncertainty('bunc', 1.0 * b_uncertainty / b, 'b')
    
    model.set_signal_processes('s*')
    return model
    

# create a model with one observable with n+1 "processes" (=n Berstein polynomial-shaped templates) with same yield. The yields
# are unknown, but sampled such that n_events_total is the total number of events.
def bernstein_model(n, n_events_total = 1000, nbins = 100):
    model = Model()
    bernstein = lambda x, nu, n: 1 if (nu,n) == (0,0) else (0 if nu > n or nu < 0 else (1 - x) * bernstein(x, nu, n-1) + x * bernstein(x, nu-1, n-1))
    for nu in range(n + 1):
        hf = HistogramFunction()
        data = [bernstein((i + 0.5) / nbins, nu, n) for i in range(nbins)]
        n_data = sum(data)
        data = [x * n_events_total / n / n_data for x in data]
        hf.set_nominal_histo(Histogram(0.0, 1.0, data))
        model.set_histogram_function('obs', 'proc%d' % nu, hf)
        model.add_lognormal_uncertainty('proc%d_unc' % nu, .1, 'proc%d' % nu)
        model.distribution.set_distribution_parameters('proc%d_unc' % nu, width = float("inf"))
    # we do not have any signal, so create one empty signal process group:
    model.signal_process_groups = {'': []}
    return model


# a 'large' model, i.e. many channels, bins and backgrounds. Used for testing the performance of
# e.g. the BB light method.
# creates a model with nchannels channels. Each channel is the same and has nbins bins from 0 to 100.
# nproc background processes are added, all have a normal shape with equidistant means at between 0 and 100 and standard deviation 20.
# The yield of each background process is given by bkgyield.
# Each process has a log-normal rate uncertainty in each channel of 10%, which amounts for nchannels * nproc nuisance parameters.
# If shpesyst is True, a total of 15 shape uncertainties are added which are treated via template morphing and slightly shift the peak
# position of all backgrounds up/down by 0.2, 0.4, ..., 3.0.
# The number of nuisance parameters excluding MC stat. treatment is 15 + nchannels * nproc.
# The relative MC stat. (bin-by-bin) uncertainty is mcstat. This uncertainty is either treated via template
# morphing (bblight=False) or via the Barlow-Beeston-light method (bblight = True).
#
# A signal process 's' is added in each channel at the bin x=40 with yield sigyield in each channel. It has no systematic uncertainties
# and is normalized to sigyield.
def large_model(nchannels = 10, nbins = 10, nproc = 5, bkgyield = 100.0, mcstat = 0.05, sigyield = 1.0, shapesyst = True, bblight = False):
    model = Model()
    inf = float("inf")
    binwidth = 100. / nbins
    # get a normal-shaped bin entries for  the range (0--100) with nbins bins:
    def get_gauss(mean, stddev, norm):
        gauss = lambda x, mean, sigma: math.exp(-(x - mean)**2 / (2.0*sigma))
        result = [gauss(binwidth * (i + 0.5), mean, stddev) for i in range(nbins)]
        n_result = sum(result)
        return [norm / n_result * r for r in result]
    # add backgrounds, same in each channel c:
    for c in range(nchannels):
        for p in range(nproc):
            hf_b = HistogramFunction()
            mean = (p + 0.5)* (100. / nproc)
            values = get_gauss(mean, 20., bkgyield)
            uncertainties = get_gauss(mean, 20., bkgyield * mcstat) if mcstat > 0.0 else None
            hf_b.set_nominal_histo(Histogram(0.0, 100.0, values, uncertainties))
            # add 15 shape uncertainties which shift the Gaussian mean up/down by up to 3.0:
            if shapesyst:
                for u in range(1, 16):
                    hf_b.set_syst_histos('shape%u' % u, Histogram(0.0, 100.0, get_gauss(mean + u*0.2, 20., bkgyield)), Histogram(0.0, 100.0, get_gauss(mean - u*0.2, 20., bkgyield)))
            model.set_histogram_function('channel%d' % c, 'bkg%d' % p, hf_b)
            model.add_lognormal_uncertainty('c%d_bkg%d_unc' % (c, p), 0.2, 'bkg%d' % p, 'channel%d' % c)
            if not bblight and mcstat > 0.0:
                 # add bin-by-bin uncertainties via template morphing. Note that this introduces a total of nbins * nproc * nchannels new nuisance parameters.
                 for i in range(nbins):
                     values_syst_plus = values[:]
                     values_syst_minus = values[:]
                     values_syst_plus[i] *= 1 + mcstat
                     values_syst_minus[i] *= 1 - mcstat
                     hf_b.set_syst_histos('mcstat_c%d_p%d_bin%d' % (c, p, i), Histogram(0.0, 100.0, values_syst_plus), Histogram(0.0, 100.0, values_syst_minus))
                     model.distribution.set_distribution('mcstat_c%d_p%d_bin%d' % (c, p, i), 'gauss', 0.0, 1.0, [-inf, inf])
        hf_s = HistogramFunction()
        hf_s.set_nominal_histo(Histogram(0.0, 100.0, get_gauss(40.0, 0.01*binwidth, sigyield)))
        model.set_histogram_function('channel%d' % c, 's', hf_s)
    # set the prior for the nuisance parameters used for the background shape systematics to a normal prior with mean 0, std.dev. 1:
    if shapesyst:
        for u in range(1, 16):
            model.distribution.set_distribution('shape%u' % u, 'gauss', 0.0, 1.0, [-inf, inf])
    model.bb_undertainties = bblight
    model.set_signal_processes('s*')
    return model
