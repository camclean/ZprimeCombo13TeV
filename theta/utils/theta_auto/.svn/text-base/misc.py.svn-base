# -*- coding: utf-8 -*-

from utils import *
from theta_interface import *
import Report
from Model import *

# base class for producers, providing consistent approach to override_parameter_distribution
# and additional_nll_term
class ProducerBase:
    def gef_cfg(self, parameters): pass
    
    # model_parameters is the list of parameters (not necesarily those of additional_nll_term)
    def get_cfg_base(self, model, signal_processes):
        result = {'name': self.name}
        if self.override_parameter_distribution is not None:
            parameters = set(model.get_parameters(signal_processes))
            if model.additional_nll_term is not None:
                parameters.update(model.additional_nll_term.get_parameters())
            result['override-parameter-distribution'] = self.override_parameter_distribution.get_cfg(parameters)
        return result
    
    def __init__(self, parameter_dist, name):
        self.override_parameter_distribution = parameter_dist
        self.name = name
    

class MleProducer(ProducerBase):
    def __init__(self, parameter_dist, name = 'mle'):
        ProducerBase.__init__(self, parameter_dist, name)
        
    # parameters_write is the list of parameters to write the mle for in the db. The default (None) means
    # to use all parameters.
    def get_cfg(self, model, signal_processes, parameters_write = None, **options):
        model_parameters = model.get_parameters(signal_processes)
        if parameters_write is None: parameters_write = model_parameters
        result = {'type': 'mle', 'minimizer': minimizer(**options), 'parameters': list(parameters_write)}
        result.update(self.get_cfg_base(model, signal_processes))
        return result


class PliProducer(ProducerBase):
    def __init__(self, parameter_dist, cls, name = 'pli', **options):
        ProducerBase.__init__(self, parameter_dist, name)
        self.cls = cls[:]
        self.parameter = 'beta_signal'
        self.options = options

    def get_cfg(self, model, signal_processes, parameters_write = None):
        result = {'type': 'deltanll_intervals', 'minimizer': minimizer(need_error = False, **self.options), 'parameter': self.parameter, 'clevels': self.cls}
        result.update(self.get_cfg_base(model, signal_processes))
        return result
        
class NllScanProducer(ProducerBase):
    def __init__(self, parameter_dist, parameter = 'beta_signal', name = 'nll_scan', range = [0.0, 3.0], npoints = 101):
        ProducerBase.__init__(self, parameter_dist, name)
        self.parameter = parameter
        self.range = range
        self.npoints = npoints

    def get_cfg(self, model, signal_processes, parameters_write = None):
        result = {'type': 'nll_scan', 'minimizer': minimizer(need_error = False), 'parameter': self.parameter,
               'parameter-values': {'start': self.range[0], 'stop': self.range[1], 'n-steps': self.npoints}}
        result.update(self.get_cfg_base(model, signal_processes))
        return result


#the class for the configuration of the "main" path and common options.
class MainBase:
    def __init__(self, root_plugins = True):
        self.cfg_options = {'plugin_files': ['$THETA_DIR/lib/core-plugins.so']}
        if root_plugins: self.cfg_options['plugin_files'].append('$THETA_DIR/lib/root.so')
        
class Run(MainBase):
    def __init__(self, n_events, data_source_dict, model_dist, root_plugins = True):
        MainBase.__init__(self, root_plugins)
        self.n_events = n_events
        self.producers = []
        self.data_source_dict = data_source_dict
        self.model_dist = model_dist
        
    def add_producer(self, producer):
        self.producers.append(producer)
        
    def get_cfg(self, model, signal_processes, **options):
        main = {'n-events': self.n_events, 'producers': [], 'log-report': False, 'output_database': sqlite_database(), 'model': "@model", 'data_source': self.data_source_dict}
        toplevel_settings = {'main': main, 'options': self.cfg_options, 'model': model.get_cfg(signal_processes)}
        model_parameters = model.get_parameters(signal_processes)
        toplevel_settings['model']['parameter-distribution'] = Distribution.merge(model.distribution, self.model_dist).get_cfg(model_parameters)
        for p in self.producers:
            toplevel_settings['p%s' % p.name] = p.get_cfg(model, signal_processes, **options)
            main['producers'].append('@p%s' % p.name)
        return toplevel_settings


def model_signal_prior_dist(input_spec):
    if input_spec.startswith('toys:') or input_spec.startswith('toys-asimov:'):
        beta_signal_value = float(input_spec[input_spec.find(':') + 1:])
    else: beta_signal_value = 0.0
    result = Distribution()
    result.set_distribution('beta_signal', 'gauss', beta_signal_value, 0.0, [beta_signal_value, beta_signal_value])
    return result


def default_signal_processes(model, signal_processes):
    if signal_processes is not None: return signal_processes
    signal_processes = [[sp] for sp in model.signal_processes]
    if len(signal_processes)==0: signal_processes.append('')
    return signal_processes

def signal_prior_dist(spec):
    result = Distribution()
    if type(spec) == str:
        if spec.startswith('flat'):
            if spec.startswith('flat:'):
                res = re.match('flat:\[([^,]+),(.*)\]', spec)
                if res is None: raise RuntimeError, "signal_prior specification '%s' invalid (does not match range specification syntax)" % spec
                xmin, xmax = float(res.group(1)), float(res.group(2))
            else:
                if spec!='flat': raise RuntimeError, "signal_prior specification '%s' invalid" % spec
                xmin, xmax = 0.0, float("inf")
            value = 0.5 * (xmax - xmin)
            if value==float("inf"): value = 1.0
            result.set_distribution('beta_signal', 'gauss', value, float("inf"), [xmin, xmax])
            return result
        elif spec.startswith('fix:'):
            v = float(spec[4:])
            result.set_distribution('beta_signal', 'gauss', v, 0.0, [v,v])
            return result
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if not isinstance(spec,Distribution): raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a Distribution instance!"
        return spec

# TODO: does not support rvobservables fully
def write_cfg2(main, model, signal_processes, method, input, id = None, **options):
    all_parameters = model.get_parameters(signal_processes)
    all_parameters = sorted(list(all_parameters))
    rvobservables = model.rvobs_distribution.get_parameters()
    theta_cfg = "parameters = " + settingvalue_to_cfg(all_parameters, 0, ['parameters']) + ";\n"
    if len(rvobservables) > 0:
        config = ""
        config += "rvobservables = " + settingvalue_to_cfg(rvobservables, 0, ['rvobservables']) + ";\n"
    obs = {}
    for o in model.observables:
        xmin, xmax, nbins = model.observables[o]
        obs[o] = {'range': [xmin, xmax], 'nbins': nbins}
    theta_cfg += "observables = " + settingvalue_to_cfg(obs, 0, ['observables']) + ";\n"
    cfg = main.get_cfg(model, signal_processes, **options)
    cfg['main']['output_database']['filename'] = '@output_name';
    for s in cfg:
        theta_cfg += s + " = " + settingvalue_to_cfg(cfg[s]) + ";\n"
    m = hashlib.md5()
    m.update(theta_cfg)
    hash = m.hexdigest()[:10]
    if id is None:
        name = '%s-%s-%s-%s' % (method, ''.join(signal_processes), input, hash)
    else:
        name = '%s-%s-%s-%s-%s' % (method, ''.join(signal_processes), input, id, hash)
    f = open(os.path.join(global_config.workdir, name + '.cfg'), 'w')
    print >>f, theta_cfg
    print >>f, 'output_name = "%s.db";\n' % name
    f.close()
    return name

def ml_fit2(model, input = 'data', signal_prior = 'flat', nuisance_constraint = 'shape:fix', signal_processes = None, n = 1, **options):
    signal_processes = default_signal_processes(model, signal_processes)
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dist(signal_prior)
    model_signal_prior = model_signal_prior_dist(input)
    data_source_dict, model_dist_signal_dict = utils.data_source_dict(model, input, **options)
    cfg_names_to_run = []
    for sp in signal_processes:
        main = Run(n, data_source_dict, model_signal_prior)
        mle = MleProducer(Distribution.merge(signal_prior, nuisance_constraint))
        main.add_producer(mle)
        name = write_cfg2(main, model, sp, 'ml_fit', input, **options)
        cfg_names_to_run.append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_: run_theta(cfg_names_to_run)
    else: return None

    cachedir = os.path.join(config.workdir, 'cache')    
    result = {}
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    
    nuisance_parameters = sorted(list(model.get_parameters('')))
    for p in nuisance_parameters:
        suffix = ''
        if nuisance_constraint.get_distribution(p)['width'] == 0.0: suffix = ' (fixed)'
        result_table.add_column(p, '%s%s' % (p, suffix))
    suffix = ''
    if signal_prior_spec.startswith('fix:'): suffix = ' (fixed)'
    result_table.add_column('beta_signal', 'beta_signal%s' % suffix)
    for icfg in range(len(cfg_names_to_run)):
        sp = signal_processes[icfg]
        name = cfg_names_to_run[icfg]
        method, sp_id, dummy = name.split('-',2)
        result[sp_id] = {}
        result_table.set_column('process', sp_id)
        parameters = set(model.get_parameters(sp))
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['mle__%s, mle__%s_error' % (p, p) for p in parameters]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
        if input != 'data':
            cols = ['source__%s' % p for p in parameters]
            source_par_values = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        i = 0
        for p in parameters:
            result[sp_id][p] = [(row[2*i], row[2*i+1]) for row in data]
            if input != 'data': result[sp_id]['source__%s' % p] = [row[i] for row in source_par_values]
            i += 1
            sorted_res = sorted([res[0] for res in result[sp_id][p]])
            n = len(sorted_res)
            if n >= 10:
                result_table.set_column(p, '%.3g (%.3g, %.3g)' % (sorted_res[int(0.5*n)], sorted_res[int(0.16*n)], sorted_res[int(0.84*n)]))
            else: result_table.set_column(p, '%.3g' % sorted_res[int(0.5*n)])
        nll_values = sql(sqlfile, 'select eventid, mle__nll from products')
        result[sp_id]['nll'] = [row[1] for row in nll_values]
        result[sp_id]['eventid'] = [row[0] for row in nll_values]
        for p in nuisance_parameters + ['beta_signal']:
            if p in parameters: continue
            result_table.set_column(p, 'n/a')
        result_table.add_row()
    config.report.new_section("Maximum Likelihood fit on ensemble '%s'" % input)
    config.report.add_p('The table entries give the median (and, if n>=10 the central 68%) of the parameter values at the found maximum of the likelihood function.')
    config.report.add_html(result_table.html())
    return result


## \brief Quantify the approximate impact of an unertainty on the result
#
# This will run
#    method(model, input = 'toys-asimov:1.0', n = 1, signal_processes = ..., nuisance_constraint = <see below>, **method_options)
# for each systematic uncertainty parameter twice, setting this parameter to +1sigma and -1sigma resp, 
# fixing all other parameters to their nominal value.
#
# It is possible to change this default of scanning -1sigma and +1sigma by setting sigma_factors.
#
# nuisance_constraint is set to a copy of model.distribution; if method_options contain
# 'nuisance_constraint', this will be considered as usual (i.e., you can set it to 'shape:fix', ...).
#
#
# returns a dictionary
# (spid) --> (parameter name) --> (sigma factor) --> (method result)
#
#
# Example: to see how the result changes if "switching on" the uncertainty in the toys and fixing the same parameter
# in the evaluation, you can use something along these line:
# for p in parameters:
#    fixp = Distribution()
#    fixp.set_distribution(p, 'gauss', 0.0, 0.0, [0., 0.])
#    individual_uncertainties(model, ..., nuisance_constraint = fixp, parameters = [p])
def individual_uncertainties(model, method, signal_processes = None, sigma_factors = [-1.0, 1.0], parameters = None,  **method_options):
    assert 'signal_processes' not in method_options
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    dist_for_method = copy.deepcopy(model.distribution)
    if 'nuisance_constraint' in method_options:
        constr = nuisance_prior_distribution(model, method_options['nuisance_constraint'])
        dist_for_method = Distribution.merge(dist_for_method, constr)
        del method_options['nuisance_constraint']
    model_dist_orig = copy.deepcopy(model.distribution)
    result = {}
    if 'n' not in method_options: method_options['n'] = 1
    for sp in signal_processes:
        spid = ''.join(sp)
        result[spid] = {}
        if parameters is None: pars = model.get_parameters(sp)
        else: pars = parameters
        for par in pars:
            if par == 'beta_signal': continue
            result[spid][par] = {}
            # fix all parameters but par:
            model.distribution = copy.deepcopy(model_dist_orig)
            for par2 in model.get_parameters(sp):
                if par == par2 or par2 == 'beta_signal': continue
                model.distribution.set_distribution_parameters(par2, width=0.0)
            # find out mean and width of par:
            dist_pars = model.distribution.get_distribution(par)
            mean, width = dist_pars['mean'], dist_pars['width']
            for sf in sigma_factors:
                model.distribution.set_distribution_parameters(par, mean = mean + sf * width, width=0.0, range = (mean + sf * width, mean + sf * width))
                res = method(model, input = 'toys-asimov:1.0', nuisance_constraint = dist_for_method, signal_processes = [sp], **method_options)
                result[spid][par][sf] = res
    model.distribution = model_dist_orig
    return result


## \brief Perform a maximum likelihood fit
#
# Finds the parameter values for all model parameters at the maximum of the likelihood, using
# the theta plugin opf type 'mle'
#
# The parameters \c input, \c signal_prior, \c nuisance_constraint, and \c signal_processes are documented in more detail
# in \ref theta_auto_params.
#
# In the report, a table with the parameter values at the maximum is created, with one entry per signal process group.
#
# \param n The number of fits to run. Values &gtr; 1 only make sense for input="toys"
# \return A dictionary with the signal process group id as key (see \ref theta_auto_signal on the definition of signal process group id).
# The value is a dictionary with the parameter name as key and a list of length n with two-tuple of floats, (parameter value, uncertainty) as value.
# Note that in case the maximization does not converge, these toys will be omitted from the result and the lists of two-tuples have length
# smaller than \c n.
def ml_fit(model, input = 'data', signal_prior = 'flat', nuisance_constraint = 'shape:fix', signal_processes = None, n = 1, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    beta_signal_value = 0.0
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dict(signal_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ('@mle',), 'output_database': sqlite_database(), 'log-report': False}
    mle = {'type': 'mle', 'name': 'mle', 'parameters': None,
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_constraint"),
       'minimizer': minimizer()}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'model-distribution-signal': delta_distribution(beta_signal = beta_signal_value), 'mle': mle, 'main': main, 'signal_prior':
        signal_prior_dict(signal_prior),  'options': cfg_options}
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input, **options)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = sorted(list(model.get_parameters(sp)))
        mle['parameters'] = model_parameters
        if 'beta_signal' in model_parameters: mle['override-parameter-distribution'] = product_distribution("@signal_prior", "@nuisance_constraint")
        else: mle['override-parameter-distribution'] = "@nuisance_constraint"
        toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'ml_fit', input, additional_settings = toplevel_settings)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')
    
    result = {}
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    
    nuisance_parameters = sorted(list(model.get_parameters('')))
    for p in nuisance_parameters:
        suffix = ''
        if nuisance_constraint.get_distribution(p)['width'] == 0.0: suffix = ' (fixed)'
        result_table.add_column(p, '%s%s' % (p, suffix))
    suffix = ''
    if signal_prior_spec.startswith('fix:'): suffix = ' (fixed)'
    result_table.add_column('beta_signal', 'beta_signal%s' % suffix)
    for i in range(len(cfg_names_to_run)):
        sp = signal_processes[i]
        name = cfg_names_to_run[i]
        method, sp_id, dummy = name.split('-',2)
        result[sp_id] = {}
        result_table.set_column('process', sp_id)
        model_parameters = model.get_parameters(sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['mle__%s, mle__%s_error' % (p, p) for p in model_parameters]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
        i = 0
        for p in model_parameters:
            result[sp_id][p] = [(row[2*i], row[2*i+1]) for row in data]
            i += 1
            sorted_res = sorted([res[0] for res in result[sp_id][p]])
            n = len(sorted_res)
            if n >= 10:
                result_table.set_column(p, '%.3g (%.3g, %.3g)' % (sorted_res[int(0.5*n)], sorted_res[int(0.16*n)], sorted_res[int(0.84*n)]))
            else: result_table.set_column(p, '%.3g' % sorted_res[int(0.5*n)])
        for p in nuisance_parameters + ['beta_signal']:
            if p in model_parameters: continue
            result_table.set_column(p, 'n/a')
        result_table.add_row()
    config.report.new_section("Maximum Likelihood fit on ensemble '%s'" % input)
    config.report.add_p('The table entries give the median (and, if n>=10 the central 68%) of the parameter values at the found maximum of the likelihood function.')
    #config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_constraint)))
    config.report.add_html(result_table.html())
    return result
    
## \brief Perform a maximum likelihood fit and get the coefficient function values for all processes / channels
#
# This is useful to scale templates to make plots. The result should not be used
# to extract cross sections or other results.
#
# options: see ml_fit
#
# returns a dictionary
# (signal process id) --> (observable name) --> (process name) --> (factor)
#
# Does not write anything to the report.
def ml_fit_coefficients(model, signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    spids = [''.join(sps) for sps in signal_processes]
    result = {}
    res = ml_fit(model, signal_processes = signal_processes, **options)
    for i in range(len(spids)):
        sp = spids[i]
        values = {}
        for param in res[sp]:
            values[param] = res[sp][param][0][0]
        result[sp] = {}
        for obs in model.observables:
            result[sp][obs] = {}
            for proc in model.get_processes(obs):
                # skip signal processes we are not interested in:
                if proc in model.signal_processes and proc not in signal_processes[i]: continue
                result[sp][obs][proc] = model.get_coeff(obs, proc).get_value(values)
                if proc in model.signal_processes:
                    result[sp][obs][proc] *= values['beta_signal']
    return result
            

## \brief Perform a KS-test on th whole model
#
# Perform a KS-test by (i) dicing toy data from the model (including nuisance parameters according to model.distribution),
# (ii) performing a maximum likelihood fit using nuisance_constraint and (iii) calculating the Kolmogorov-Smirnov test-statistic value comparing the
# prediction and data after the fit.
#
# This method always uses the background-only model. Therefore, it has no signal_priori or signal_processes parameters.
# If you want to include signal in the KS-test, you have to define it as background (via Model.set_signal_processes) before calling\
# this function.
#
# Returns the p-value for data. Does not write anything to the report.
def ks_test(model, n = 5000, nuisance_constraint = 'shape:fix', **options):
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    main = {'n-events': n, 'model': '@model', 'producers': ('@mle',), 'output_database': sqlite_database(), 'log-report': False}
    mle = {'type': 'mle', 'name': 'mle', 'parameters': [], 'write_ks_ts': True,
       'override-parameter-distribution': product_distribution("@nuisance_constraint"),
       'minimizer': minimizer(need_error = False)}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'mle': mle, 'main': main, 'options': cfg_options}
    cfg_names_to_run = []
    model_parameters = sorted(list(model.get_parameters('')))
    toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
    id = None
    if 'id' in options: id = options['id']
    # for toys:
    input = 'toys:0'
    main['data_source'], dummy = data_source_dict(model, input)
    name = write_cfg(model, '', 'ks_test', input, id=id, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    # for data:
    input = 'data'
    main['data_source'], dummy = data_source_dict(model, input)
    main['n-events'] = 1
    name = write_cfg(model, '', 'ks_test', input, id=id, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')

    name_toys = cfg_names_to_run[0]
    sqlfile = os.path.join(cachedir, '%s.db' % name_toys)
    data = sql(sqlfile, 'select mle__ks_ts from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    ks_values_bkg = [row[0] for row in data]
    
    name_data = cfg_names_to_run[1]
    sqlfile = os.path.join(cachedir, '%s.db' % name_data)
    data = sql(sqlfile, 'select mle__ks_ts from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    ks_value_data = data[0][0]
    p = len([x for x in ks_values_bkg if x >= ks_value_data]) * 1.0 / len(ks_values_bkg)
    return p
    
## \brief perform a chi2-based goodness-of-fit test
#
# This will toss n toys according to "input" (you should set this to 'toys:1.0' to include signal or "toys:0.0" for background-only).
# For each toy, a maximum likelihood fit is performed, according to the model and the given nuisance_constraint. After the fit,
# the chi2 value is computed. Comparing the chi2 distribution for these toys with the one computed after the fit on data yields a p-value
# which is returned by this function.
def chi2_test(model, signal_process_group, input = 'toys:1.0', n = 5000, signal_prior = 'flat', nuisance_constraint = '', **options):
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    main = {'n-events': n, 'model': '@model', 'producers': ('@mle',), 'output_database': sqlite_database(), 'log-report': False}
    mle = {'type': 'mle', 'name': 'mle', 'parameters': [], 'write_pchi2': True,
       'override-parameter-distribution': product_distribution("@nuisance_constraint", "@signal_prior"),
       'minimizer': minimizer(need_error = False)}
    signal_prior = signal_prior_dict(signal_prior)
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'mle': mle, 'main': main, 'options': cfg_options, 'signal_prior': signal_prior}
    cfg_names_to_run = []
    model_parameters = sorted(list(model.get_parameters(signal_process_group)))
    toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
    # for toys:
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    name = write_cfg(model, signal_process_group, 'chi2_test', input, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    # for data:
    input = 'data'
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    main['n-events'] = 1
    name = write_cfg(model, signal_process_group, 'chi2_test', input, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run, **options)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')

    name_toys = cfg_names_to_run[0]
    sqlfile = os.path.join(cachedir, '%s.db' % name_toys)
    data = sql(sqlfile, 'select mle__pchi2 from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    chi2_values_bkg = [row[0] for row in data]
    #print chi2_values_bkg
    
    name_data = cfg_names_to_run[1]
    sqlfile = os.path.join(cachedir, '%s.db' % name_data)
    data = sql(sqlfile, 'select mle__pchi2 from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    chi2_value_data = data[0][0]
    #print chi2_value_data
    p = len([x for x in chi2_values_bkg if x >= chi2_value_data]) * 1.0 / len(chi2_values_bkg)
    return p


# range and npoints define the scan points for beta_signal
#
# returns a map
# (signal process id) --> (list of pd instances)
def nll_scan(model, input='data', n = 1, nuisance_constraint = '', signal_prior = 'flat', signal_processes = None, par_range = [0.0, 3.0], npoints = 101, **options):
    signal_processes = default_signal_processes(model, signal_processes)
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dist(signal_prior)
    model_signal_prior = model_signal_prior_dist(input)
    data_source_dict, model_dist_signal_dict = utils.data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        main = Run(n, data_source_dict, model_signal_prior)
        nll = NllScanProducer(Distribution.merge(signal_prior, nuisance_constraint), range = par_range, npoints = npoints)
        main.add_producer(nll)
        name = write_cfg2(main, model, sp, 'nll_scan', input)
        cfg_names_to_run.append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_: run_theta(cfg_names_to_run)
    else: return None
    
    cachedir = os.path.join(config.workdir, 'cache')
    result = {}
    for i in range(len(cfg_names_to_run)):
        name = cfg_names_to_run[i]
        method, sp_id, dummy = name.split('-',2)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select nll_scan__nll from products')
        pds = [plotdata_from_histoColumn(row[0]) for row in data]
        for pd in pds: pd.as_function = True
        result[sp_id] = pds
    return result


## \brief Perform a KS compatibility test for all individual channels of the given model
#
# Each channel is treated independently. For each channel in the model, a "restricted" model is built
# which contains each only that channel. This restricted model and all options are passed to to ks_test, see there
# for valid options.
#
# If fit_yield is True, one parameter is added to the model which is used to scale all processes
# at the same time (=the fraction between processes is not modified, at least not by this parameter).
# In this case, after the fit, the normalisation of simulation and data coincides by construction
# and the KS-test effectively compares the shapes only, not the rate.
# Note that even if fit_yield is False, there is still a maximum likelihood fit done which finds the
# parameter values at the maximum of the likelihood function, using nuisance_constraints (details see ks_test).
# Depending on the model, this can already mean that the yield is fitted.
#
# Returns a dictionary (channel name) --> (p-value).
#
# Does not write anything to the report.
def ks_test_individual_channels(model, fit_yield = False, **options):
    observables = model.observables.keys()
    result = {}
    for obs in observables:
        model_obs = copy.deepcopy(model)
        model_obs.restrict_to_observables([obs])
        if fit_yield:
            model_obs.distribution.set_distribution('scale_', 'gauss', 1.0, inf, (0, inf))
            procs = model_obs.get_processes(obs)
            for p in procs:
                model_obs.get_coeff(obs, p).add_factor('id', parameter = 'scale_')
        options['id'] = obs
        result[obs] = ks_test(model_obs, **options)
    return result
    

## \brief Calculate profile likelihood intervals
#
# Calculate profile likelihood intervals using the deltanll_intervals plugin in %theta.
#
# For the \c input, \c nuisance_constraint, \c signal_processes parameters, see \ref theta_auto_params.
#
# \param n The number of toys. A value n &gt; 1 only makes sense for input='toys'
# \param cls A list of confidence levels
# \param write_report If True, a summary will be written to the report. If None, the report will be written if and only if input is 'data'
#
# \return A dictionary with the signal process group id as key. The value is a dictionary with the confidence levels as keys; the values are lists of
#  two-tuples (lower interval border, upper interval border). In addition to the confidence levels, there is a special key 'mle' which contains a list
#  of beta_signal values at the maximum.
#
# For example, if the only signal process is called 'signal' and cl=[0.68, 0.95], the 1sigma intervals are
# \code
#  result['signal'][0.68][i] # i=0..n-1
# \endcode
# and the 2signma intervals are
# \code
#  result['signal'][0.95][i] # i=0..n-1
# \endcode
# and the beta_signal values at the maximum are 
# \code
#  result['signal']['mle'][i] # i=0..n-1
# \endcode
#
# In case the minimization fails, the lists have less than n entries.
#
# If write_report is true, it will only write out the result of the first fit. This usually makes only sense for data,
# not for toy MC.
def pl_intervals(model, input = 'toys:0', n = 100, signal_prior = 'flat', nuisance_constraint = '', cls = [0.90], signal_processes = None, write_report = None, **options):
    signal_processes = default_signal_processes(model, signal_processes)
    if write_report is None: write_report = input == 'data'
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dist(signal_prior)
    model_signal_prior = model_signal_prior_dist(input)
    data_source_dict, model_dist_signal_dict = utils.data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        main = Run(n, data_source_dict, model_signal_prior)
        pl = PliProducer(Distribution.merge(signal_prior, nuisance_constraint), cls, **options)
        main.add_producer(pl)
        name = write_cfg2(main, model, sp, 'pl_intervals', input)
        cfg_names_to_run.append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_: run_theta(cfg_names_to_run)
    else: return None

    result_table = Report.table()
    result_table.add_column('signal process')
    result_table.add_column('mle', 'mle')
    for cl in cls: result_table.add_column('cl%g' % cl, 'confidence level %g' % cl)
    cachedir = os.path.join(config.workdir, 'cache')
    col_suffixes= ['%05d' % (cl*10000) for cl in cls]
    result = {}
    for i in range(len(cfg_names_to_run)):
        name = cfg_names_to_run[i]
        method, sp_id, dummy = name.split('-',2)
        result_table.set_column('signal process', sp_id)
        result[sp_id] = {'mle': []}
        for cl in cls: result[sp_id][cl] = []
        model_parameters = model.get_parameters(sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        colnames = []
        for cs in col_suffixes:
            colnames.append('pli__lower%s' % cs)
            colnames.append('pli__upper%s' % cs)
        data = sql(sqlfile, 'select pli__maxl, %s from products' % (', '.join(colnames)))
        if len(data)==0: raise RuntimeError, "no result (fit not coverged?)"
        first_row = True
        for row in data:
            result[sp_id]['mle'].append(row[0])
            if first_row: result_table.set_column('mle', '%g' % row[0])
            for icol in range(len(cls)):
                interval = (row[2*icol+1], row[2*icol+2])
                result[sp_id][cls[icol]].append(interval)
                if first_row:
                    result_table.set_column('cl%g' % cls[icol], '(%g, %g)' % interval)
            first_row = False
        result_table.add_row()
    if write_report:
        config.report.new_section("deltanll intervals")
        config.report.add_html(result_table.html())
    return result


# calculate the expected limits using asimov data (i.e., without dicing nuisance parameters and without dicing Poisson)
# returns a dictionary (signal process id) --> (expected limit on beta_signal).
#
# Note that the PL method is an asymptotic method, so this method is mainly useful to compare different scenarios, but not for the final result.
def get_expected_pl_limits(model, input = 'toys-asimov:0.0', signal_processes = None, **options):
    dist_orig = model.distribution
    model.distribution = get_fixed_dist(dist_orig)
    res = pl_intervals(model, input = input, n=1, signal_processes = signal_processes, nuisance_constraint = dist_orig, **options)
    result = {}
    for spid in res:
        result[spid] = res[spid][0.9][0][1]
    model.distribution = dist_orig
    return result

# runs deltanll_intervals and measures coverage as a function of the true beta_signal
# Returns: dictionary: (spid) --> (true beta signal) --> |  (cl) --> coverage
#                                                        |  'successrate' --> fit success rate
def pl_coveragetest(model, beta_signal_values = [0.2*i for i in range(10)], n = 1000, write_report = True, **deltanll_options):
    result = {}
    result_tables = {}
    for beta_signal in beta_signal_values:
        res = pl_intervals(model, input='toys:%g' % beta_signal, n = n, **deltanll_options)
        for spid in res:
            cls = [k for k in res[spid].keys() if type(k)==float]
            if spid not in result: result[spid] = {}
            if spid not in result_tables:
                 result_tables[spid] = Report.table()
                 result_tables[spid].add_column('beta_signal', 'true beta_signal')
                 for cl in cls: result_tables[spid].add_column('coverage %g' % cl, 'Coverage for cl=%g' % cl)
                 result_tables[spid].add_column('fit success fraction')
            result[spid][beta_signal] = {}
            result_tables[spid].set_column('beta_signal', '%g' % beta_signal)
            icl = 0
            for cl in cls:
                n_covered = len([1 for i in range(len(res[spid][cl])) if res[spid][cl][i][0] <= beta_signal and res[spid][cl][i][1] >= beta_signal])
                n_total = len(res[spid][cl])
                coverage = n_covered*1.0 / n_total
                result[spid][beta_signal][cl] = coverage
                successrate = (n_total*1.0 / n)
                result[spid][beta_signal]['successrate'] = successrate
                result_tables[spid].set_column('coverage %g' % cl, '%g' % coverage)
                result_tables[spid].set_column('fit success fraction', '%g' % successrate)
            result_tables[spid].add_row()
    if write_report:
        for spid in result_tables:
            config.report.new_section("deltanll interval coverage test for signal process '%s'" % spid)
            config.report.add_html(result_tables[spid].html())
    return result
