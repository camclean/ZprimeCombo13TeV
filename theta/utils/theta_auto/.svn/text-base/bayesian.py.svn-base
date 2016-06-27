# -*- coding: utf-8 -*-
import config, os.path, math
#from Distribution import Distribution
from theta_interface import *
import Report

from utils import *


## \brief Compute bayesian posterior quantiles for a given ensemble
#
# Uses the theta plugin mcmc_quantiles.
#
# options:
#   * 'toydata_seed' [default: 1] for another random seed used in toy data generation 
#   * 'run_theta' [default: True] whether to really run theta. Otherwise, only the cfg files are written. In case run_theta is False, None is returned
#   * 'mcmc_iterations' [default: 20000]
#   * 'mcmc_seed': the random sed to use for the mcmc_producer. The actual seed used for producer i is mcmc_seed + i
#      with i from 0 to mcmc_nproducers-1.
#
# Report: Writes a table containing the quantiles for each signal process to the report. For data, it will contain
# the best estimate and uncertainty of the determined limit (the uncertainty is only available for n > 1). For input!='data',
# the +-1sigma and +-2sigma (i.e., the centrgl 68% and central 84%) are given.
#
# \return A dictionary (signal process group id) --> 'quantiles' | 'accrates' --> (list of quantiles or acceptance rates).
def bayesian_quantiles(model, input = 'toys:0', n = 1000, signal_prior = 'flat', nuisance_prior = '', quantile = 0.95, write_report = True, signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    # allow n==0 to simplify some scripts ...
    if n==0: return {}
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ['@bayes_interval'], 'output_database': sqlite_database(), 'log-report': False}
    if int(options.get('n_threads', 0)) > 0:
        main['type'] = 'run_mt'
        main['n_threads'] = int(options['n_threads'])
    bayes_interval = {'type': 'mcmc_quantiles', 'name': 'bayes', 'parameter': 'beta_signal', 'diag': True, 're-init': 1,
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior"), 'quantiles': [quantile], 'iterations': 20000 }
    if 'mcmc_iterations' in options: bayes_interval['iterations'] = options['mcmc_iterations']
    toplevel_settings = {'signal_prior': signal_prior, 'main': main}
    options['load_root_plugins'] = False
    toplevel_settings.update(get_common_toplevel_settings(**options))
    if 'mcmc_seed' in options: bayes_interval['rnd_gen'] = {'seed': options['mcmc_seed'] }
    bayes_interval['name'] = 'bayes'
    toplevel_settings['bayes_interval'] = bayes_interval
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'bayesian_quantiles', input, additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')
    plotsdir = os.path.join(config.workdir, 'plots')
    if not os.path.exists(plotsdir): os.mkdir(plotsdir)
    
    result = {}
    for name in cfg_names_to_run:
        method, sp_id, dummy = name.split('-',2)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select bayes__quant%0.5d, bayes__accrate from products' % (quantile * 10000))
        result[sp_id] = {'quantiles': [row[0] for row in data], 'accrates': [row[1] for row in data]}
    return result


## \brief Populate an instance of plotutil.plotdata with results from bayesian_quantiles.
#
# 'quantiles' is the return value from \ref bayesian_quantiles.
#
# include_band is a boolean indicating whether to also include the +-1sigma and +-2sigma bands, or only the median line.
#
# options:
# - signalprocess_to_value: a dictionary mapping signal process name to values to be used as x axis for the band plot. As default, the first integer
#     in the signal process name is used for the x axis value.
#
# returns one plotutil.plotdata instance containing the 'observed' (or median expected) limit and 'expected' bands.
def limit_band_plot(quantiles, include_band, quantile = 0.95, **options):
    #expected results maps (process name) -> (median, band1, band2)
    # where band1 and band2 are tuples for the central 68 and 95%, resp.
    results = {}
    signal_processes = set()
    for sp in quantiles:
        # ignore uncertainties or other special entries:
        if '__' in sp: continue
        signal_processes.add(sp)
    # map process names and x values:
    x_to_sp = get_x_to_sp(signal_processes, **options)
    pd = plotdata()
    pd.color = '#000000'
    if include_band: pd.color = '#aaaaaa'
    pd.as_function = True
    pd.x = sorted(list(x_to_sp.keys()))
    pd.y = []
    if include_band:
        pd.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
    else: pd.yerrors = []
    for x in pd.x:
        sp = x_to_sp[x]
        data = sorted(quantiles[sp]['quantiles'])
        n = len(data)
        median, band1, band2 = (data[n / 2], (data[int(0.16 * n)], data[int(0.84 * n)]), (data[int(0.025 * n)], data[int(0.975 * n)]))
        if not include_band:
            mean, error = get_trunc_mean_width(data)
            pd.y.append(mean)
            pd.yerrors.append(error / math.sqrt(len(data)))
            continue
        pd.y.append(median)
        pd.bands[1][0].append(band1[0])
        pd.bands[1][1].append(band1[1])
        pd.bands[0][0].append(band2[0])
        pd.bands[0][1].append(band2[1])
    return pd

# TODO: move report_limit_band_plot out of bayesian: it is also used from cls_limits!

## \brief Make expected / observed limits plots and write them to the report
#
# name and shortname are used to distinguish different calls to this function. name is used ion the report
# and should be the 'human-firendly' version, while shortname is used for filenames and should be the 'computer-friendly'
# version.
def report_limit_band_plot(expected_limits, observed_limits, name, shortname, write_table = True):
    plotsdir = os.path.join(config.workdir, 'plots')
    plots = []
    extra_legend_items = []
    if expected_limits is not None:
        expected_limits.legend = 'median expected limit'
        extra_legend_items.append((expected_limits.bands[0][2], '$\\pm 2\\sigma$ expected limit'))
        extra_legend_items.append((expected_limits.bands[1][2], '$\\pm 1\\sigma$ expected limit'))
        plots.append(expected_limits)
    if observed_limits is not None:
        observed_limits.legend = 'observed limit'
        plots.append(observed_limits)
    if len(plots) == 0: return
    config.report.new_section('Limits %s' % name)
    if write_table:
        result_table = Report.table()
        result_table.add_column('process', 'signal process')
        if expected_limits is not None:
            result_table.add_column('exp', 'expected limit')
            result_table.add_column('exp1', 'expected limit (central 1sigma)')
            result_table.add_column('exp2', 'expected limit (central 2sigma)')
        if observed_limits is not None:
            result_table.add_column('obs', 'observed limit')
        x_values = []
        if expected_limits is not None: x_values = expected_limits.x
        else: x_values = observed_limits.x
        for i in range(len(x_values)):
            result_table.set_column('process', '%g' % x_values[i])
            if expected_limits is not None:
                result_table.set_column('exp', '%.3g' % expected_limits.y[i])
                result_table.set_column('exp1', '%.3g--%.3g' % (expected_limits.bands[1][0][i], expected_limits.bands[1][1][i]))
                result_table.set_column('exp2', '%.3g--%.3g' % (expected_limits.bands[0][0][i], expected_limits.bands[0][1][i]))
            if observed_limits is not None:
                result_table.set_column('obs', '%.3g +- %.3g' % (observed_limits.y[i], observed_limits.yerrors[i]))
            result_table.add_row()
        config.report.add_html(result_table.html())
    plot(plots, 'signal process', 'upper limit', os.path.join(plotsdir, 'limit_band_plot-%s.png' % shortname), extra_legend_items=extra_legend_items)
    plot(plots, 'signal process', 'upper limit', os.path.join(plotsdir, 'limit_band_plot-log-%s.png' % shortname), logy = True, extra_legend_items=extra_legend_items)
    config.report.add_html('<p><img src="plots/limit_band_plot-%s.png" /></p>' % shortname)
    config.report.add_html('<p><img src="plots/limit_band_plot-log-%s.png" /></p>' % shortname)
    return plots
    

## \brief Calculate Bayesian limits.
#
# This is a high-level interface to calculate expected and observed limits and make a limit band plot:
# it will calculate expected and observed limits (using \ref bayesian_quantiles), make a "limit vs. mass" band plot
# (using \ref limit_band_plot) and write the result to the report (using \ref report_limit_band_plot).
#
# The 'what' parameter controls which limits are computed. Valid vaues are:
# * 'observed': compute observed limit on data only
# * 'expected': compute +-1sigma, +-2sigma limit bands for background only
# * 'all': both 'data' and 'expected'
#
# options are as for bayesian_quantiles, with the modification that
# * there should be n_toy, n_obs  instead of n;  'n' will be ignored
# * 'input' will be ignored
#
# returns a tuple of two plotutil.plotdata instances. The first contains expected limit (including the band) and the second the 'observed' limit.
# If 'what' is not 'all', one of the tuple entries is None.
def bayesian_limits(model, what = 'all', **options):
    if 'n' in options: del options['n']
    if 'input' in options: del options['input']
    n = options.get('n_toy', 1000)
    plot_expected, plot_observed = None, None
    if what in ('expected', 'all'):
        expected_limits = bayesian_quantiles(model, input = 'toys:0', n = n, **options)
        if expected_limits is not None:
            plot_expected = limit_band_plot(expected_limits, True)
    if what in ('observed', 'all'):
        assert model.has_data()
        n = options.get('n_data', 10)
        observed_limits = bayesian_quantiles(model, input = 'data', n = n, **options)
        if observed_limits is not None:
            plot_observed = limit_band_plot(observed_limits, False)
    # catch the case where the routines return None (e.g., if run_theta = False)
    if plot_expected is not None and plot_observed is not None:
        report_limit_band_plot(plot_expected, plot_observed, 'Bayesian', 'bayesian')
    return (plot_expected, plot_observed)


#
# sp_group is a list of strings: the signal process names to scale with beta_signal
# options: same as bayesian_quantiles (in particular signal_processes!):
# * 'input' is ignored / replaced by settings according to beta_signal_range, beta_signal_n
#def bayesian_limits_coveragetest(model, beta_signal_range = [0.0, 10.0], beta_signal_n = 10, **options): pass
    

## \brief Calculate the marginal posterior of the given parameters
#
#
# histogram_specs is a dictionary of (parameter name) -> tuple(int nbins, float xmin, float max) and determines for which
# parameters the poterior is computed on which range and binning. Note that the computational complexity is roughly
# proportional to nbins * mcmc_iterations. Therefore, use large values only if you have to / for the final result. It is suggested
# to use 30 bins as a start and mcmc_iterations = 10000.
#
# returns a dictionary (signal process group id) --> (parameter name) -> (list of plotutil.plotdata)
#
# Writes the posterior plots to the report (TODO: more info ...)
def posteriors(model, histogram_specs, input = 'data', n = 3, signal_prior = 'flat', nuisance_prior = '', signal_processes = None, mcmc_iterations = 10000, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ('@posteriors',), 'output_database': sqlite_database(), 'log-report': False}
    posteriors = {'type': 'mcmc_posterior_histo', 'name': 'posteriors', 'parameters': [],
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
       'smooth': options.get('smooth', True), 'iterations': mcmc_iterations }
    for par in histogram_specs:
        posteriors['parameters'].append(par)
        nbins, xmin, xmax = histogram_specs[par]
        posteriors['histo_%s' % par] = {'range': [float(xmin), float(xmax)], 'nbins': nbins}
    toplevel_settings = {'signal_prior': signal_prior, 'posteriors': posteriors, 'main': main}
    options['load_root_plugins'] = False
    toplevel_settings.update(get_common_toplevel_settings(**options))
    cfg_names_to_run = []
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'posteriors', input, additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    
    cachedir = os.path.join(config.workdir, 'cache')
    plotsdir = os.path.join(config.workdir, 'plots')
    
    result = {}
    config.report.new_section('Posteriors %s' % input)
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    parameters = sorted([par for par in histogram_specs])
    for par in parameters:
        result_table.add_column('maximum posterior %s' % par)
    for name in cfg_names_to_run:
        method, sp, dummy = name.split('-',2)
        config.report.add_html('<h2>For signal "%s"</h2>' % sp)
        result_table.set_column('process', sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['posteriors__posterior_%s' % par for par in parameters]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        data = [map(plotdata_from_histoColumn, row) for row in data]
        result[sp] = {}
        i = 0
        for par in parameters:
            result[sp][par] = [row[i] for row in data]
            for pd in result[sp][par]: pd.as_function = True
            plot(result[sp][par], par, 'posterior density', os.path.join(plotsdir, '%s-%s.png' % (name, par)))
            config.report.add_html('<p>%s:<br/><img src="plots/%s-%s.png" /></p>' % (par, name, par))
            i += 1
            maxima = sorted(map(argmax, result[sp][par]))
            result_table.set_column('maximum posterior %s' % par, '%.3g' % maxima[int(0.5 * len(maxima))])
        result_table.add_row()
    config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_prior)))
    config.report.add_html(result_table.html())
    return result

