# -*- coding: utf-8 -*-
import math, bisect
import scipy.special

from theta_interface import *
import frequentist


# container for toys made for the CLs construction
class truth_ts_values:
    def __init__(self):
        self.truth_to_ts_sb = {}
        self.truth_to_ts_b = {}

    def truth_values(self):
        result = set(self.truth_to_ts_sb.keys())
        if 0.0 in result: result.remove(0.0)
        return result

    def add_point_b(self, truth, ts):
        if truth not in self.truth_to_ts_b: self.truth_to_ts_b[truth] = []
        self.truth_to_ts_b[truth].append(ts)

    def add_point_sb(self, truth, ts):
        if truth not in self.truth_to_ts_sb: self.truth_to_ts_sb[truth] = []
        self.truth_to_ts_sb[truth].append(ts)


    def get_truth_to_ts(self, q):
        assert q > 0 and q < 1, str(q)
        result = {}
        for t in self.truth_values():
            b = sorted(self.truth_to_ts_b[t])
            ts = b[int(q*len(b))]
            result[t] = ts
        return result
        
    # return a tuple (cls+b, clb, cls) of plotutil instances which are
    # the CLs+b, clb and cls curves as a function of the "truth" (beta_signal)
    def get_cl_vs_truth(self, truth_to_ts):
        pd_clsb, pd_clb, pd_cls = plotutil.plotdata(as_function = True, legend = 'CLs+b', color = '#0000ff'), plotutil.plotdata(as_function = True, legend = 'CLb'), plotutil.plotdata(as_function = True, color = '#ff0000', legend = 'CLs')
        pd_clsb.x = sorted(list(self.truth_values()))
        pd_clb.x, pd_cls.x = pd_clsb.x, pd_clsb.x
        for t in pd_clb.x:
            ts = truth_to_ts[t]
            n_b = len(self.truth_to_ts_b[t])
            clb = len([x for x in self.truth_to_ts_b[t] if x < ts]) * 1.0 / n_b
            clb_error = max(math.sqrt(clb * (1 - clb) / n_b), 1.0 / n_b)
            pd_clb.y.append(clb)
            n_sb = len(self.truth_to_ts_sb[t])
            clsb = len([x for x in self.truth_to_ts_sb[t] if x < ts]) * 1.0 / n_sb
            clsb_error = max(math.sqrt(clsb * (1 - clsb) / n_sb), 1.0 / n_sb)
            pd_clsb.y.append(clsb)
            if clb==0: clb += clb_error
            cls = clsb / clb
            pd_cls.y.append(cls)
            cls_error = math.sqrt((clsb_error / clb)**2 + (clb_error * cls / clb)**2)
        return pd_clsb, pd_clb, pd_cls
            

# make some plots of the test statistic for the report, using the given db filename
#
# returns a list of four plotutil.plotdata instances which contain curves as a function of the truth value (beta_signal):
# * cls+b for data
# * clb for data
# * cls for data
# * expected CLs curve, including 1sigma and 2sigma bands
def debug_cls_plots(dbfile, ts_column = 'lr__nll_diff'):
    limits = sql(dbfile, 'select "index", "limit", limit_uncertainty from cls_limits')
    indices = [row[0] for row in limits]
    have_data = 0 in indices
    data = sql(dbfile, 'select runid, eventid, lr__poi, source__truth, "%s" from products order by runid, eventid' % ts_column)
    tts = truth_ts_values()
    truth_to_ts_data = {}
    for row in data:
        if row[0] == 0 and row[1] == 0:
            truth_to_ts_data[row[2]] = row[4]
            continue
        if row[3] > 0:  tts.add_point_sb(row[3], row[4])
        else: tts.add_point_b(row[2], row[4])
    plotsdir = os.path.join(config.workdir, 'plots')
    config.report.new_section('debug_cls for file %s' % dbfile)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts_data)
    
    expected_pd_cls = plotutil.plotdata(as_function = True)
    expected_pd_cls.x = pd_cls.x[:]
    # use median as "expected":
    truth_to_ts = tts.get_truth_to_ts(0.5)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    expected_pd_cls.y[:] = pd_cls.y
    #build the bands:
    band1s = [[], [], '#00ff00']
    band2s = [[], [], '#ffff00']
    truth_to_ts = tts.get_truth_to_ts(0.025)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band2s[0][:] = pd_cls.y[:]
    truth_to_ts = tts.get_truth_to_ts(0.975)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band2s[1][:] = pd_cls.y[:]
    truth_to_ts = tts.get_truth_to_ts(0.16)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band1s[0][:] = pd_cls.y[:]
    truth_to_ts = tts.get_truth_to_ts(0.84)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band1s[1][:] = pd_cls.y[:]        
    expected_pd_cls.bands = [band2s, band1s]    
    
    #draw all:
    pds = [pd_clsb, pd_clb, pd_cls, expected_pd_cls]
    plotname = 'debug_cl-vs-truth-data.png'
    plot(pds, 'beta_signal', 'p-value', os.path.join(plotsdir, plotname))
    config.report.add_html('<p>For data: <br/><img src="plots/%s" /></p>' % plotname)
    
    return pds
        

def asymptotic_cls_limits(model, use_data = True, signal_process_groups = None, beta_signal_expected = 0.0, bootstrap_model = True, input = None, n = 1, options = None, as_plotdata = True):
    """
    Calculate CLs limits using asymptotic formulae.
    
    Options:
    
    * ``use_data`` - if ``True``, also calculate observed limit.
    * ``beta_signal_expected`` - signal strength value to use to calculate the expected limit bands. The default of 0.0 corresponds to limits expected for background-only.
      If set to ``None``, no expected limit will be computed.
    * ``bootstrap_model`` - if this is set to ``True``  -- and ``use_data`` is ``True`` -- the parameter values are fitted to data first.
    * ``as_plotdata`` - if set to ``True``, returns the result as ``plotdata`` objects. Otherwise, returns the expected and observed limits dictionary ``signal -> list of values``.
    
    
    Also note that some parameters described in :ref:`common_parameters` have a special meaning here:
    
    * ``input`` - this is the data source to calculate the "observed" limit(s) for. The default of ``None`` is equivalent to "data" if ``use_data==True`` and to ``None`` (not computing any observed limit) if ``use_data=False``
    * ``n`` is the number of "observed" limits to calculate from the ``input`` data source. Only has effect if ``input`` is not the default ``None``.

    For ``signal_process_groups`` and ``options`` refer to :ref:`common_parameters`. Note that the common options ``signal_prior``, ``nuisance_prior``
    are missing on purpose as the asymptotic method implemented in theta only works for flat priors.
    
    Just like :meth:`cls_limits`, the return value is a two-tuple ``(pd_expected, pd_observed)`` of plotutil.plotdata instances that contain the
    expected and observed limits, including the 1sigma and 2sigma expected limit bands. If more than one "observed" limit is calculated,
    these limits are used to calculate 1sigma and 2sigma bands in ``pd_observed`` as well.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    use_data = use_data and model.has_data()
    if bootstrap_model and use_data:
        model = frequentist.get_bootstrapped_model(model, options = options)
    else:
        model = frequentist.frequentize_model(model)
    if input is None:
        if use_data:
            input = 'data'
            n = 1
        else:
            n = 0
    limits_expected = {}
    limits_observed = {}
    n_expected = 5
    if beta_signal_expected is None: n_expected = 0
    for spid, signal_processes in signal_process_groups.iteritems():
        r = AsymptoticClsMain(model, signal_processes, input = input, beta_signal_expected = beta_signal_expected, n = n)
        r.run_theta(options)
        limits = r.get_results('limits', ['limit'], 'index')['limit']
        limits_observed[spid] = limits[n_expected:] # can be empty in case of use_data = False
        limits_expected[spid] = limits[:n_expected]
    if not as_plotdata:
        return limits_expected, limits_observed
    # convert to plotdata:
    spids = signal_process_groups.keys()
    x_to_sp = get_x_to_sp(spids)
    pd_expected, pd_observed = None, None
    if beta_signal_expected is not None:
        pd_expected = plotdata(color = '#000000', as_function = True, legend = 'expected limit')
        pd_expected.x = sorted(list(x_to_sp.keys()))
        pd_expected.y  = []
        pd_expected.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
    if n > 0:
        pd_observed = plotdata(color = '#000000', as_function = True, legend = 'observed limit')
        pd_observed.x = sorted(list(x_to_sp.keys()))
        pd_observed.y = []
        if n > 1:
            pd_observed.bands = [([], [], '#0000ff'), ([], [], '#0000aa')]
    for x in sorted(x_to_sp.keys()):
        sp = x_to_sp[x]
        if beta_signal_expected is not None:
            pd_expected.bands[0][0].append(limits_expected[sp][0])
            pd_expected.bands[1][0].append(limits_expected[sp][1])
            pd_expected.y.append(limits_expected[sp][2])
            pd_expected.bands[1][1].append(limits_expected[sp][3])
            pd_expected.bands[0][1].append(limits_expected[sp][4])
        if n > 0:
            nobs = len(limits_observed[sp]) # note: can be < n, if some failed.
            lobs = sorted(limits_observed[sp])
            median = lobs[nobs/2]
            pd_observed.y.append(median)
            if nobs > 1:
                # first band is 2sigma band:
                pd_observed.bands[0][0].append(lobs[int(0.025 * nobs)])
                pd_observed.bands[0][1].append(lobs[int(0.975 * nobs)])
                # second is 1sigma:
                pd_observed.bands[1][0].append(lobs[int(0.16 * nobs)])
                pd_observed.bands[1][1].append(lobs[int(0.84 * nobs)])
    report_limit_band_plot(pd_expected, pd_observed, 'Asymptotic CLs', 'acls')
    return pd_expected, pd_observed

def cls_limits(model, use_data = True, signal_process_groups = None, nuisance_prior = None, frequentist_bootstrapping = False,
 cls_options = {}, seed = None, options = None):
    """
    Calculate CLs limits, based on toys.
    
    Options:
    
    * ``use_data`` - if ``True``, also calculate observed limits
    * ``frequentist_bootstrapping`` - if ``True``, do a fit to data first and use the parameter values at the best fit for the toys.
    * ``cls_options`` is a dictionary of CLs-specific options with the following keys:
       * "expected_bands" - number of toys to make for the expected limit bands (default: 2000)
       * "clb_cutoff" - the lowest allowed CLb value for the expected limit before giving up (default: 0.02)
       * "reltol_limit" - relative accuracy for the CLs limit: More toys will be done until this accracy is reached (default: 0.05)
       * "input_expected" - a input specification (see the discussion of the ``input`` parameter in :ref:`common_parameters`) to use for calculating the expected limit band. The default is "toys:0.0".
    * ``seed`` is a random seed. The default value ``None`` uses a different seed each time.
    
    For the options ``signal_process_groups``, ``nuisance_prior`` and ``options`` refer to :ref:`common_parameters`.
    
    Returns a tuple of two plotutil.plotdata instances. The first contains expected limit (including the bands) and the second the 'observed' limit.
    If ``use_data`` is ``False``, the second plotdata instance is ``None``.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    cls_options = dict(cls_options)
    cls_options['expected_bands'] = int(cls_options.get('expected_bands', 2000))
    cls_options['frequentist_bootstrapping'] = frequentist_bootstrapping
    cls_options['ts_column'] = 'dnll__nll_diff'    
    input = 'data' if use_data else None
    
    if frequentist_bootstrapping and use_data:
        model = frequentist.get_bootstrapped_model(model)
        
    # dictionaries from spid to
    # * tuple (limit, uncertainty) for observed
    # * list of limit for expected
    observed_limit, expected_limits = {}, {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = ClsMain(model, signal_processes, signal_prior = 'flat', input = input,
            producers = [DeltaNllHypotest(model, signal_processes, nuisance_prior, restrict_poi = 'beta_signal', signal_prior_sb = 'flat', signal_prior_b = 'flat')],
            seed = seed)
        r.set_cls_options(**cls_options)
        r.run_theta(options)
        data = r.get_results('cls_limits', ['index', 'limit', 'limit_uncertainty'])
        expected_limits[spid] = []
        for i in range(len(data['index'])):
            if data['index'][i] == 0: observed_limit[spid] = data['limit'][i], data['limit_uncertainty'][i]
            else: expected_limits[spid].append(data['limit'][i])
        expected_limits[spid].sort()
        if len(expected_limits[spid])==0: del expected_limits[spid]

    spids = signal_process_groups.keys()
    x_to_sp = get_x_to_sp(spids)
    pd_expected, pd_observed = None, None
    if len(expected_limits) > 0:
        pd_expected = plotdata(color = '#000000', as_function = True, legend = 'expected limit')
        pd_expected.x = sorted(list(x_to_sp.keys()))
        pd_expected.y  = []
        pd_expected.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
    if use_data:
        pd_observed = plotdata(color = '#000000', as_function = True, legend = 'observed limit')
        pd_observed.x = sorted(list(x_to_sp.keys()))
        pd_observed.y = []
        pd_observed.yerrors = []
    
    for x in sorted(x_to_sp.keys()):
        sp = x_to_sp[x]
        if pd_expected:
            limits = expected_limits[sp]
            n = len(limits)
            median, low_1s, high_1s, low_2s, high_2s = limits[int(0.5*n)], limits[int(0.16*n)], limits[int(0.84*n)], limits[int(0.05*n)], limits[int(0.95*n)]
            pd_expected.y.append(median)
            pd_expected.bands[1][0].append(low_1s)
            pd_expected.bands[1][1].append(high_1s)
            pd_expected.bands[0][0].append(low_2s)
            pd_expected.bands[0][1].append(high_2s)
        if use_data:
            observed, observed_unc = observed_limit[sp]
            pd_observed.y.append(observed)
            pd_observed.yerrors.append(observed_unc)
    report_limit_band_plot(pd_expected, pd_observed, 'CLs', 'cls')
    return pd_expected, pd_observed
