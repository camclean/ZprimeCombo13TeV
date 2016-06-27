# -*- coding: utf-8 -*-

from utils import *
from theta_interface import *
import plotutil
import itertools

def mle(model, input, n, with_error = True, with_covariance = False, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None,
 signal_prior = 'flat', ks = False, chi2 = False, pred = False, all_columns = False, options = None):
    """
    Find the maximum likelihood estimate for all model parameters.
    
    For the parameters ``model``, ``input``, ``n``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys``, ``signal_prior`` and ``options`` refer to :ref:`common_parameters`.
    
    Parameters:
    
    * ``with_error`` if ``True``, return the uncertainty for each parameter, as determined by migrad. If you do not need the uncertainty, you can set this to ``False`` which can
      increase speed and robustness of the minimizer.
    * ``with_covariance`` if ``True``, return the uncertainty for each parameter, as determined by migrad.
    * ``ks`` - calculate the Kolmogorov-Smirnov test statistic for each toy *after* the fit has run. Note that the ks test statistic will use the maximum value of all ks values
      over all channels which sometimes is not what you want. To access the ks test values, use ``key = "__ks"`` (see below).
    * ``chi2`` - calculate the chi-square test statistic for each toy *after* the fit has run; the chi-square value is generalized for Poisson statistics and is twice the logarithm
      of the Poisson likelihood ratio in each bin, using the model prediction and the observed data in the numerator and denominator for this ratio. This chi-square value is summed
      over all bin in all channels. To access the chi-square test values, use ``key = "__chi2"`` (see below).
    * ``pred`` - if ``True``, also return the model prediction as evaluated by theta at the parameter values maximizing the likelihood
    
    The return value respects the convention described in :ref:`return_values`. Specifically, this method returns a nested python dictionary with the following keys, in this order:
    
    * ``spid`` - the signal process id, i.e., the keys of ``model.signal_process_groups``, or (if given) the keys of ``signal_process_groups`` (see :ref:`what_is_signal` for details)
    * ``key`` - a parameter name or some special underscore names '__nll', '__ks' (only if ``ks = True``), '__chi2' (only if ``chi2 = True``), or '__cov' (only if ``with_covariance=True``)
    
    Using these two keys, the value is a list of length ``n`` of results. In the case ``key`` is a parameter name, each item of this list is
    a two-tuple ``(value, uncertainty)`` which are the parameter value at the minimum and the uncertainty (from migrad) for this parameter, resp. In case ``with_error``
    is ``False``, ``uncertainty`` is always ``None``. If ``key`` is one of the special underscore names, the value is a list where each entry in the list is one floating point value.
    
    The covariance matrix is a n*n matrix where n is the number of parameters the model depends on. To get the parameters for this matrix use::
    
     parameters = model.get_parameters(signal_process_group)
     
    where ``signal_process_group`` is a list of names of the signal processes (see :ref:`what_is_signal` for the definition of a signal process group).
    
    .. warning:: The covariance matrix is passed as evaluated by the underlying minimizer. This minimizer is ROOT's MINUIT2 by default which produces completely nonsensical results in some cases.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [MleProducer(model, signal_processes, nuisance_constraint, signal_prior, need_error = with_error, with_covariance = with_covariance, ks = ks, chi2 = chi2, write_prediction = pred)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products()
        parameters = [m.group(1) for m in map(lambda colname: re.match('mle__(.+)_error', colname), res.keys()) if m is not None]
        result[spid] = {}
        for p in parameters:
            if not with_error: result[spid][p] = list(itertools.izip_longest(res['mle__%s' % p], []))
            else: result[spid][p] = zip(res['mle__%s' % p], res['mle__%s_error' % p])
        result[spid]['__nll'] = res['mle__nll']
        if with_covariance: result[spid]['__cov'] = map(matrix_from_dbblob, res['mle__covariance'])
        if chi2: result[spid]['__chi2'] = res['mle__pchi2']
        if ks: result[spid]['__ks'] = res['mle__ks_ts']
        if pred:
            for obs in model.get_observables():
                result[spid]['__pred_%s' % obs] = map(histogram_from_dbblob, res['mle__pred_%s' % obs])
        """
        if all_columns:
            for key in res:
                if not key.startswith('mle__'):
                    result[spid]['__' + key] = res[key]
        """
    return result


def pl_interval(model, input, n, cls = [cl_1sigma, cl_2sigma], signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None,
parameter = 'beta_signal'):
    """
    Evaluate profile likelihood intervals for the given parameter. The profile liklihood function is defined by maxzimizing the "full" likelihood function w.r.t. all parameters
    but the parameter of interest (specified with ``parameter``). In the asymptotic case, Wilks' Theorem can be used to derive intervals for this parameters by finding the
    values of the parameter where the negative-log of the profile likelihood function takes a certain difference (which depends on the confidence level only)
    from the minimum value. This method is therefore also refered to as "Delta Log-Likelihood method" and is a generalization of the "Delta chi-square method".
    
    For the parameters ``model``, ``input``, ``n``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys``, ``signal_prior`` and ``options`` refer to :ref:`common_parameters`. Parameters specific to this method are:
    
    * ``cls`` - a list of confidence levels to calculate the interval for. The default is ``[cl_1sigma, cl_2sigma]`` which determines the 1sigma and 2sigma intervals.
    * ``parameter`` - a string which refers to the model parameter to evalaute the interval(s) for.

    The return value respects the convention described in :ref:`return_values`. Specifically, this method returns a nested python dictionary with the following keys, in this order:
    
    * ``spid`` - the signal process group id, i.e., the keys of ``model.signal_process_groups``, or (if given) the keys of ``signal_process_groups`` (see :ref:`what_is_signal` for details)
    * ``cl`` - the confidence level. Valid values are elements of ``cls`` and ``0.0`` which is always calculated; it is the maximum likelihood estimate for the parameter.
    
    Using these two keys, the value is a list of length ``n`` of results. Each item of this list is a two-tuple ``(lower, upper)`` which are the lower and upper end of the interval, respectively. The exception is ``cl=0.0`` where the list elements are floats, no two-tuples.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    colnames = ['pli__lower%05d' % int(cl*10000 + 0.5) for cl in cls] + ['pli__upper%05d' % int(cl*10000 + 0.5) for cl in cls] + ['pli__maxl']
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [PliProducer(model, signal_processes, nuisance_constraint, cls = cls, parameter = parameter, signal_prior = signal_prior)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(colnames)
        result[spid] = {0.0: res['pli__maxl']}
        for i, cl in enumerate(cls):
            result[spid][cl] = zip(res[colnames[i]], res[colnames[len(cls) + i]])
    return result



def pl_interval_coveragetest(model, spid, beta_signal_values = [0.1*i for i in range(21)], n = 1000, cl = cl_1sigma, nuisance_prior_toys = None, nuisance_constraint = None, options = None):
    """
    Performs toys with different values for ``beta_signal`` and determine the coverage probability of the profile liklihood interval.
    If the profile likelihood interval correctly covers, the coverage probability reported should be ``cl``. Note that if you are close to the
    physical boundary of ``beta_signal = 0.0``, it is expected that the method overcovers.
    
    For the parameters ``model``, ``nuisance_constraint``, ``nuisance_prior_toys``, and ``options`` refer to :ref:`common_parameters`.
    
    Parameters:
    * ``spid`` - the signal process group id, i.e., a key of ``model.signal_process_groups`` (see :ref:`what_is_signal` for details)
    * ``beta_signal_values`` - a list of floating point values to run the coverage test for
    * ``n`` - the number of toys to perform for each value of ``beta_signal``
    * ``cl`` - the confidence level
    
    The return value is a ``plotdata`` instance. The x-values are the beta_signal_values and the y-values are the
    estimated coverage probabilities for the confidence level, i.e., the fraction of toys in which the true value was contained in the interval.
    The yerrors are the uncertainty due to the limited number of toys performed at each point.
    """
    if options is None: options = Options()
    signal_process_groups = {spid: model.signal_process_groups[spid]}
    pd = plotdata(as_function = True)
    pd.yerrors = []
    for v in beta_signal_values:
        intervals = pl_interval(model, input = 'toys:%g' % v, n = n, cls = [cl], signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint, nuisance_prior_toys = nuisance_prior_toys)[spid][cl]
        pd.x.append(v)
        p_cov, p_cov_error = get_p(count(lambda x: x[0] <= v and x[1]>=v, intervals), len(intervals))
        pd.y.append(p_cov)
        pd.yerrors.append(p_cov_error)
    return pd


def nll_scan(model, input, n, npoints=100, range = [0.0, 3.0], adaptive_startvalues = True, parameter = 'beta_signal', signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None):
    """
    Evalaute the profile likelihood function in ``parameter`` at the values specified by ``npoints`` and ``range``.
    
    For the parameters ``model``, ``n``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys``, ``signal_prior``, and ``options`` refer to :ref:`common_parameters`.
    
    Parameters:
    
    * ``parameter`` - the model parameter the profile likelihood function is defined in; all other parameters will be "minimized out"
    * ``npoints``, ``range`` - the number of points and the range for ``parameter`` to evaluate the profile likelihood function at: ``parameter`` is fixed in turn to values
      ``range[0] + (range[1] - range[0]) / npoints * i`` -- where ``i`` runs from 0 to ``npoints - 1`` --, and all other parameters are "minimized out".
    * ``adaptive_startvalues`` - if ``True``, the minimizer will use the result of the previous fit as start value for the fit at the next scan point.
       This usually increases convergence speed, but might cause problems for some models in which the minimization stops too early. This can lead to an artificial structure
       of the profile likelihood function. If set to ``False``, the fit always starts at the same point (which corresponds to the most aprioi most probable parameter value).
    
    The return value is a nested dictionary: the first-level key is the signal process group id (see :ref:`what_is_signal` for a definition). The value
    is a list of length ``n`` of :class:`~theta_auto.plotutil.plotdata` instances, containing the negative profile log-likelihood values in the scan.
    In addition to those equidistant points of the scan, the point at the minimum of the negative log-likelihood is added, with a y-value of exactly 0.0.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [NllScanProducer(model, signal_processes, nuisance_constraint, npoints = npoints, range = range,
                 parameter = parameter, signal_prior = signal_prior, adaptive_startvalues = adaptive_startvalues)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(['nllscan__nll', 'nllscan__maxl'])
        histos = map(histogram_from_dbblob, res['nllscan__nll'])
        result[spid] = []
        for h, maxl in zip(histos, res['nllscan__maxl']):
            pd = plotutil.plotdata()
            pd.set_histogram(h)
            # insert the minimum:
            imin = 0
            while pd.x[imin] < maxl and imin < len(pd.x): imin += 1
            pd.x.insert(imin, maxl)
            pd.y.insert(imin, 0.0)
            result[spid].append(pd)
    return result
    
def zvalue_approx(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, options = None, eventid_info = False, signal_prior_sb = 'flat'):
    """
    Calculate Z values (significance in sigma), using Wilks' Theorem. The reported approximate Z-values are
    
    .. math::
    
      Z = \\sqrt{2 \\log \\frac{L_{s+b}}{L_b}}
      
    where both nominator and denominator are maximized: for L_s+b, all parameters are free, including the signal strength parameter ``beta_signal``. For L_b,
    ``beta_signal`` is set to 0.0. In both cases, the likelihood contains terms according to ``nuisance_contraint``.
    
    For the parameters ``model``, ``input``, ``n``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys``,
    and ``options`` refer to :ref:`common_parameters`.
    
    The return value respects the convention described in :ref:`return_values`. Specifically, this method returns a nested python dictionary with two levels;
    the keys are:
    
    * ``spid`` - the signal process group id, i.e., a key of ``model.signal_process_groups`` (see :ref:`what_is_signal` for details)
    * "Z", "__runid", or "__eventid". The latter two are only available in case ``eventid_info = True``
    
    The value is a list of the approximate Z values calculated as decribed above.
    
    To convert Z-values to p-values or vice versa, use :meth:`p_to_Z` and :meth:`Z_to_p`.
    
    .. note:: As a result of the definition as the above square root, the Z-value is never negative, and any p-value will be larger than 0.5.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        p = DeltaNllHypotest(model, signal_processes, nuisance_constraint, signal_prior_sb = signal_prior_sb)
        r = Run(model, signal_processes = signal_processes, signal_prior = 'flat', input = input, n=n, producers = [p], nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        result[spid] = {}
        data = r.get_products(['dnll__nll_diff', 'runid', 'eventid'])
        result[spid]['Z'] = [p_to_Z(scipy.stats.chi2.sf(2 * x, 1) / 2) for x in data['dnll__nll_diff']]
        if eventid_info:
            result[spid]['__runid'] = data['runid']
            result[spid]['__eventid'] = data['eventid']
    return result

