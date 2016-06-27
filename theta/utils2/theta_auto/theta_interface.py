# -*- coding: utf-8 -*-

import hashlib, io, os, time, re, numpy, shutil, termios, threading, StringIO
from Model import *
from utils import *
import config

# In general, each (plugin) class in theta (which in turn corresponds to a setting group with a "type" setting) corresponds to a class here.
#
# Conventions:
# * each class representing a concrete theta plugin should implement a constructor which takes all required arguments to build the class completely.
# * All classes are immutable in the sense that no change should be done after they are constructed. In case manipulation is needed,
#   use methods which return modified copies.
# * Note that any model instances, etc. passed to the modules here is assumed to be unmodified until theta is run.
#   This allows the modules to just save references instead of performing copies.
# * the theta configuration should be returned as dictionary by a method get_cfg(options) where options is an instance of Options class



# Conventions for theta config file generation:
# * there is only one model relevant to a config file, which is written to the path 'model', s.t. other modules
#   can refer to it via "@model"
# * the default nusiance parameter distribution (for toy generation and likelihood construction) is the one
#   from the model (model.distribution), but it can be overridden (even for a subset of parameters) for both toy generation and for
#   the likelihood construction in the producers.


from ConfigParser import SafeConfigParser


# Options to control certain details of the theta configuration generation which do not change the actual modules run etc. but generate
# somewhat different theta configurations. Examples are
# * number of threads to use
# * default minimizer options, such as whether or not to use mcmc_minimizer and how
# * whether or not to use llvm
# * some fine-tuning options for cls_limits, such as the number of toys for background-only and clb_cutoff, etc.
#
# The reason not to implement these as options for the methods are that the Options instance can be passed around more easily,
# even if new options are added, or old ones removed.
#
# Options affecting the actual statistical / physical _meaning_ of what should be done (and not
# just some details of how the theta config is generated)
# should NOT go here. Rather, they should be standard arguments for the respective methods.
#
class Options(SafeConfigParser):
    def __init__(self):
        SafeConfigParser.__init__(self)
        self.default_config = """
# global options concerning all / multiple modules or theta itself
[global]
debug = False
check_cache = True
        
# minimizer options
[minimizer]
always_mcmc = False
mcmc_iterations = 1000
bootstrap_mcmcpars = 0
strategy = fast
minuit_tolerance_factor = 1

[newton]
debug = 0
#-1 means default from theta plugin
step_cov = -1
maxit = -1
par_eps = -1
force_cov_positive = False
second_pass = False
use_nll_der = False
       
[cls_limits]
write_debuglog = True
        
[model]
use_llvm = False
use_tbb = False
tbb_nthreads = 0
robust_nll = False
        
[main]
n_threads = 1

[mcmc]
strategy = asimov_widths
stepsize_factor = None
"""
        self.readfp(io.BytesIO(self.default_config))
        
    def get_workdir(self):
        return os.path.realpath(config.workdir)
    
    def copy(self):
        s = StringIO.StringIO()
        self.write(s)
        result = Options()
        result.read(s)
        return result

# each class representing a theta module should inherit from ModuleBase
class ModuleBase:
    def __init__(self):
        self.submodules = []
        self.required_plugins = frozenset(['core-plugins.so'])
    
    # derived classes should usually set self.required_plugins instead of re-implementing this function.
    # One exception is the case in which the plugins required depend on options.
    def get_required_plugins(self, options):
        plugins = set(self.required_plugins)
        for m in self.submodules:
            plugins.update(m.get_required_plugins(options))
        return plugins
        
    # declare module to be a submodule of the current module. This is used to track plugin
    # dependencies.
    def add_submodule(self, module):
        self.submodules.append(module)


class DataSource(ModuleBase):
    # input_string is either 'toys-asimov:XXX', 'toys:XXX', (where 'XXX' is a floating point value, the value of beta_signal), 'data',
    # or 'replay-toys:XXX' where XXX is the filename. It is also valid to directly specify the filename of the .db file.
    #
    # override_distribution is a Distribution which can override the default choice in the 'toys' case. If None, the default is to
    # * use model.distribution in case of input_string = 'toys:XXX'
    # * use the model.distrribution fixed to the most probable parameter values in case of input_string = 'toys-asimov:XXX'
    def __init__(self, input_string, model, signal_processes, override_distribution = None, name = 'source', seed = None):
        ModuleBase.__init__(self)
        flt_regex = r'([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
        if input_string.endswith('.db') and not input_string.startswith('replay-toys'):
            input_string = 'replay-toys:' + input_string
        self.name = name
        self.seed = seed
        self._id = input_string
        self.signal_processes = signal_processes
        self.model = model
        if input_string.startswith('toys'):
            self.mode = 'toys'
            self.beta_signal = None
            self.asimov = False
            match_toys = re.match('toys:%s' % flt_regex, input_string)
            if match_toys is not None:
                self.beta_signal = float(match_toys.group(1))
                if override_distribution is None:
                    self.distribution = model.distribution
                else:
                    self.distribution = Distribution.merge(model.distribution, override_distribution)
            else:
                match_toys_asimov = re.match('toys-asimov:%s' % flt_regex, input_string)
                if match_toys_asimov is None: raise RuntimeError, "invalid input specification '%s'" % input_string
                self.asimov = True
                self.beta_signal = float(match_toys_asimov.group(1))
                if override_distribution is None:
                    self.distribution = get_fixed_dist(model.distribution)
                else:
                    self.distribution = Distribution.merge(get_fixed_dist(model.distribution), override_distribution)
        elif input_string.startswith('replay-toys:'):
            self.mode = 'replay-toys'
            self.fname = input_string[len('replay-toys:'):]
            if not os.path.exists(self.fname): raise RuntimeError, "specified '%s' as data source, but file '%s' does not exist!" % (input_string, self.fname)
            self._id = 'replay-toys:' + self.fname[self.fname.rfind('-')+1:self.fname.rfind('.db')]
        else:
            if input_string != 'data': raise RuntimeError, "invalid input specification '%s'" % input_string
            self.mode = 'data'
            self.data_histos = {}
            for o in model.get_observables():
                self.data_histos[o] = model.get_data_histogram(o)
            self.data_rvobsvalues = model.get_data_rvobsvalues()
            
            
    def get_id(self): return self._id
    
    def get_cfg(self, options):
        result = {'name': self.name}
        if self.mode == 'toys':
            seed = self.seed
            if seed is None: seed = -1
            parameters = self.model.get_parameters(self.signal_processes)
            bkg_parameters = set(parameters)
            bkg_parameters.discard('beta_signal')
            dist = {'type': 'product_distribution', 'distributions': [self.distribution.get_cfg(bkg_parameters)]}
            if 'beta_signal' in parameters: dist['distributions'].append({'type': 'delta_distribution', 'beta_signal' : self.beta_signal})
            result.update({'type': 'model_source', 'model': '@model', 'override-parameter-distribution': dist, 'rnd_gen': {'seed': seed}})
            if self.asimov:
                result['dice_poisson'] = False
                result['dice_template_uncertainties'] = False
                result['dice_rvobs'] = False
            return result
        elif self.mode == 'replay-toys':
            result.update({'type': 'replay_toys', 'input_database': {'type': 'sqlite_database_in', 'filename': self.fname}})
            return result
        else:
            result['type'] = 'histo_source'
            for o in self.data_histos: result[o] = self.data_histos[o].get_cfg(False)
            if len(self.data_rvobsvalues) > 0:
                result['rvobs-values'] = dict(self.data_rvobsvalues)
            return result


# a data source which adds many root_ntuple_sources. This is suitable for drawing toy datasets which are random subsets of
# events of a TTree which in turn is mainly useful for correlation studies ...
#
# Usage::
#
#  source = RootNtupleSource(model)
#  # Usually, the branch name is assumed to be the same as the observable name. If this is not the case, you have to explicitly state that:
#  source.define_observable_branchname(obsname = 'nnout', branchname = 'discriminator')
#  # add some files. Usually, there is one file per 'process'. You can specify optional total_nevents and weight_branchname to use
#  source.add_file('ttbar.root', mean_nevents = 107.2, relweight_branchname = 'weight2')
#
# The default for mean_nevents is the sum of 
#
class RootNtupleSource(ModuleBase):
    def __init__(self, model, name = 'source'):
        ModuleBase.__init__(self)
        self.required_plugins = frozenset(['root.so', 'core-plugins.so'])
        self.name = name
        # each entry in files is a dictionary with the keys 'filename', 'mean_nevents', etc.
        self.files = []
        # the 'observables' configuration for the root_ntuple_source:
        self.observables_cfg = {}
        for obs in model.get_observables():
            xmin, xmax, nbins = model.get_range_nbins(obs)
            self.observables_cfg[obs] = {'branchname': obs, 'nbins': nbins, 'range': (xmin, xmax)}
        
        
    def define_observable_branchname(self, obsname, branchname):
        if obsname not in self.observables_cfg: raise RuntimeError, "unknown observable '%s'" % obsname
        self.observables_cfg[obsname]['branchname'] = branchname
    
    # using None will use the defaults from the root_ntuple_source plugin.
    def add_file(self, filename, mean_nevents = None, relweight_branchname = None, treename = None, seed = None):
        d = {'filename': filename}
        if mean_nevents is not None: d['mean_nevents'] = float(mean_nevents)
        if relweight_branchname is not None: d['relweight_branchname'] = str(relweight_branchname)
        if treename is not None: d['treename'] = str(treename)
        if seed is not None: d['seed'] = int(seed)
        
        
    def get_cfg(self, options):
        result = {'name': self.name, 'type': 'add_sources'}
        result['sources'] = []
        for i, f in enumerate(self.files):
            s_cfg = {'name': '%s%d' % (self.name, i), 'type': 'root_ntuple_source', 'observables': self.observables_cfg}
            s_cfg['filename'] = f['filename']
            for key in ('treename', 'mean_nevents', 'relweight_branchname'):
                if key in f: s_cfg[key] = f[key]
            if 'seed' in f: s_cfg['rnd_gen'] = {'seed': f['seed']}
            result.sources.append(s_cfg)

# base class for producers, providing consistent approach to override_parameter_distribution and name
class ProducerBase(ModuleBase):
    def get_cfg_base(self, options):
        result = {'name': self.name}
        if self.override_distribution_cfg is not None:
            result['override-parameter-distribution'] = self.override_distribution_cfg
        return result
    
    
    # note that signal_prior is only respected if override_distribution is not None.
    # if model is None, override_distribution and signal_prior are ignored.
    def __init__(self, model, signal_processes, override_distribution, name, signal_prior):
        ModuleBase.__init__(self)
        assert type(name) == str
        self.override_distribution_cfg = None
        self.name = name
        if model is not None:
            parameters = set(model.get_parameters(signal_processes))
            if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
            else: dist = model.distribution
            if 'beta_signal' in parameters:
                signal_prior_dist = _signal_prior_dist(signal_prior)
                dist = dist.copy()
                dist = Distribution.merge(dist, signal_prior_dist)
                self.override_distribution_cfg = dist.get_cfg(parameters)
            else:
                self.override_distribution_cfg = dist.get_cfg(parameters)
    

class MleProducer(ProducerBase):
    # parameters_write is the list of parameters to write the mle for in the db. The default (None) means
    # to use all parameters.
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'mle', need_error = True, with_covariance = False, parameters_write = None, ks = False, chi2 = False, write_prediction = False):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.minimizer = Minimizer(need_error, with_covariance)
        self.add_submodule(self.minimizer)
        self.ks = ks
        self.chi2 = chi2
        self.with_covariance = with_covariance
        self.write_prediction = write_prediction
        if parameters_write is None: self.parameters_write = sorted(list(model.get_parameters(signal_processes)))
        else: self.parameters_write = sorted(list(parameters_write))
        
    def get_cfg(self, options):
        result = {'type': 'mle', 'minimizer': self.minimizer.get_cfg(options), 'parameters': list(self.parameters_write)}
        if self.with_covariance: result['write_covariance'] = True
        if self.ks: result['write_ks_ts'] = True
        if self.chi2: result['write_pchi2'] = True
        if self.write_prediction: result['write_prediction'] = True
        result.update(self.get_cfg_base(options))
        return result
        
        
class QuantilesProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'quant', parameter = 'beta_signal', quantiles = [0.16, 0.5, 0.84], iterations = 10000, seed = 0):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.parameter = parameter
        self.quantiles = quantiles
        self.iterations = iterations
        self.seed = seed
        
    def get_cfg(self, options):
        strategy = options.get('mcmc', 'strategy')
        stepsize_factor = options.get('mcmc', 'stepsize_factor')
        result = {'type': 'mcmc_quantiles', 'parameter': self.parameter, 'quantiles': self.quantiles}
        result.update(self.get_cfg_base(options))
        result['mcmc_strategy'] = {'type': strategy, 'name': self.name + "_mcs", 'iterations': self.iterations}
        if stepsize_factor != 'None':  result['mcmc_strategy']['factor'] = float(stepsize_factor)
        if self.seed != 0: result['mcmc_strategy']['rnd_gen'] = {'seed': self.seed}
        return result
        
        
        
class PosteriorProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution, histogram_specs, signal_prior = 'flat', name = 'post', iterations = 10000, smooth = True):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.histogram_specs = histogram_specs
        self.iterations = iterations
        self.smooth = smooth
        
    def get_cfg(self, options):
        result = {'type': 'mcmc_posterior_histo', 'iterations': self.iterations, 'smooth': self.smooth, 'parameters': self.histogram_specs.keys()}
        for (p, (nbins, xmin, xmax)) in self.histogram_specs.iteritems():
            result['histo_%s' % p] = {'range': [float(xmin), float(xmax)], 'nbins': nbins}
        result.update(self.get_cfg_base(options))
        return result
    
class MCMCMeanPredictionProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'mp', iterations = 10000):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.iterations = iterations
        self.observables = sorted(list(model.get_observables()))
        
    def get_cfg(self, options):
        result = {'type': 'mcmc_mean_prediction', 'iterations': self.iterations, 'observables': self.observables}
        result.update(self.get_cfg_base(options))
        return result
    


class PDWriter(ProducerBase):
    def __init__(self, name = 'pdw'):
        ProducerBase.__init__(self, model = None, signal_processes = None, override_distribution = None, name = name, signal_prior = None)
        
    def get_cfg(self, options):
        result = self.get_cfg_base(options)
        result['type'] = 'pseudodata_writer'
        result['write-data'] = True
        return result


class PliProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, cls = [0.6827, 0.95], name = 'pli', parameter = 'beta_signal', signal_prior = None):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.cls = [float(s) for s in cls]
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)
        self.parameter = parameter

    def get_cfg(self, options):
        result = {'type': 'deltanll_intervals', 'minimizer': self.minimizer.get_cfg(options), 'parameter': self.parameter, 'clevels': self.cls}
        result.update(self.get_cfg_base(options))
        return result
        
        
class DeltaNllHypotest(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, name = 'dnll', restrict_poi = None, restrict_poi_value = None, signal_prior_sb = 'flat', signal_prior_b = 'fix:0.0', write_pchi2 = False):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = None, name = name, signal_prior = None)
        self.minimizer = Minimizer(need_error = False)
        self.restrict_poi = restrict_poi
        self.write_pchi2 = write_pchi2
        self.restrict_poi_value = restrict_poi_value
        self.add_submodule(self.minimizer)
        if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
        else: dist = model.distribution
        sb_parameters = set(model.get_parameters(signal_processes))
        b_parameters = set(model.get_parameters(''))
        means = dist.get_means()
        # the "signal parameters": those which the model depends only for the signal ...
        means_spar = {}
        for p in means:
            if p in sb_parameters and p not in b_parameters: means_spar[p] = means[p]
        dist_bkg = Distribution.merge(dist, get_fixed_dist_at_values(means_spar))
        self.sb_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_sb)]}
        self.b_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist_bkg.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_b)]}
        
        
    def get_cfg(self, options):
        result = {'type': 'deltanll_hypotest', 'minimizer': self.minimizer.get_cfg(options),
           'background-only-distribution': self.b_distribution_cfg, 'signal-plus-background-distribution': self.sb_distribution_cfg}
        if self.write_pchi2: result['write_pchi2'] = True
        result.update(self.get_cfg_base(options))
        if 'override-parameter-distribution' in result: del result['override-parameter-distribution']
        if self.restrict_poi is not None: result['restrict_poi'] = self.restrict_poi
        if self.restrict_poi_value is not None: result['default_poi_value'] = self.restrict_poi_value
        return result
    
class NllDerProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, name = 'nll_der', signal_prior_sb = 'flat', signal_prior_b = 'fix:0.0', parameter = 'beta_signal'):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = None, name = name, signal_prior = None)
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)
        self.parameter = parameter
        if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
        else: dist = model.distribution
        sb_parameters = set(model.get_parameters(signal_processes))
        assert self.parameter in sb_parameters, "Parameter %s invalid: model does not depend on it" % parameter
        b_parameters = set(model.get_parameters(''))
        means = dist.get_means()
        # the "signal parameters": those which the model depends only for the signal ...
        means_spar = {}
        for p in means:
            if p in sb_parameters and p not in b_parameters: means_spar[p] = means[p]
        dist_bkg = Distribution.merge(dist, get_fixed_dist_at_values(means_spar))
        self.sb_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_sb)]}
        self.b_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist_bkg.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_b)]}

    def get_cfg(self, options):
        result = {'type': 'nll_der', 'minimizer': self.minimizer.get_cfg(options),
           'background-only-distribution': self.b_distribution_cfg, 'signal-plus-background-distribution': self.sb_distribution_cfg, 
           'parameter': self.parameter}
        result.update(self.get_cfg_base(options))
        if 'override-parameter-distribution' in result: del result['override-parameter-distribution']
        return result

class MCMCRatioProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, name = 'mcmcratio', signal_prior_sb = 'fix:1.0', signal_prior_b = 'fix:0.0', iterations = 10000):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = None, name = name, signal_prior = None)
        if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
        else: dist = model.distribution
        sb_parameters = set(model.get_parameters(signal_processes))
        b_parameters = set(model.get_parameters(''))
        means = dist.get_means()
        # the "signal parameters": those which the model depends only for the signal ...
        means_spar = {}
        for p in means:
            if p in sb_parameters and p not in b_parameters: means_spar[p] = means[p]
        dist_bkg = Distribution.merge(dist, get_fixed_dist_at_values(means_spar))
        self.sb_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_sb)]}
        self.b_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist_bkg.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_b)]}
        self.iterations = iterations        
        
    def get_cfg(self, options):
        result = {'type': 'mcmc_posterior_ratio', 'background-only-distribution': self.b_distribution_cfg, 'signal-plus-background-distribution': self.sb_distribution_cfg, 'iterations': self.iterations}
        result.update(self.get_cfg_base(options))
        if 'override-parameter-distribution' in result: del result['override-parameter-distribution']
        return result
       
class NllScanProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, signal_prior = 'flat', name = 'nllscan', parameter = 'beta_signal', range = [0.0, 3.0], npoints = 101, adaptive_startvalues = False):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = override_distribution, name = name, signal_prior = signal_prior)
        self.parameter = parameter
        self.range = range
        self.npoints = npoints
        self.adaptive_startvalues = adaptive_startvalues
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)

    def get_cfg(self, options):
        result = {'type': 'nll_scan', 'minimizer': self.minimizer.get_cfg(options), 'parameter': self.parameter,
               'parameter-values': {'start': self.range[0], 'stop': self.range[1], 'n-steps': self.npoints}, 'adaptive_startvalues': self.adaptive_startvalues}
        result.update(self.get_cfg_base(options))
        return result


class Minimizer(ModuleBase):
    def __init__(self, need_error, need_covariance = None):
        ModuleBase.__init__(self)
        self.need_error = need_error
        if need_covariance is None:
            self.need_covariance = self.need_error
        else:
            self.need_covariance = need_covariance
        
    def get_required_plugins(self, options):
        plugins = set(['core-plugins.so'])
        strategy = options.get('minimizer', 'strategy')
        if strategy == 'newton_vanilla':  pass
        else: plugins.add('root.so')
        return plugins
        
    
    def get_cfg(self, options):
        always_mcmc = options.getboolean('minimizer', 'always_mcmc')
        mcmc_iterations = options.getint('minimizer', 'mcmc_iterations')
        strategy = options.get('minimizer', 'strategy')
        tolerance_factor = options.getfloat('minimizer', 'minuit_tolerance_factor')
        assert strategy in ('fast', 'robust', 'minuit_vanilla', 'newton_vanilla', 'lbfgs_vanilla', 'tminuit')
        if strategy == 'tminuit':
            result = {'type': 'root_minuit1'}
            if self.need_covariance: result['hesse'] = True
            return result
        elif strategy == 'minuit_vanilla':
            result = {'type': 'root_minuit'}
        elif strategy == 'newton_vanilla':
            result = {'type': 'newton_minimizer'}
            if self.need_covariance: result['improve_cov'] = True
            debug = options.getint('newton', 'debug')
            if debug > 0: result['debug'] = debug
            step_cov = options.getfloat('newton', 'step_cov')
            if step_cov > 0: result['step_cov'] = step_cov
            par_eps = options.getfloat('newton', 'par_eps')    
            if par_eps > 0: result['par_eps'] = par_eps
            maxit = options.getint('newton', 'maxit')
            if maxit > 0: result['maxit'] = maxit
            force_cov_positive = options.getboolean('newton', 'force_cov_positive')
            if force_cov_positive: result['force_cov_positive'] = True
            use_nll_der = options.getboolean('newton', 'use_nll_der')
            if use_nll_der: result['use_nll_der'] = True
            second_pass = options.getboolean('newton', 'second_pass')
            if second_pass: result['second_pass'] = True
            return result
        elif strategy == 'lbfgs_vanilla':
            result = {'type': 'lbfgs_minimizer'}
            return result
        else:
            minimizers = []
            #try, in this order: migrad, mcmc+migrad, simplex, mcmc+simplex, more mcmc+simplex
            if not always_mcmc: minimizers.append({'type': 'root_minuit'})
            minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min0', 'rnd_gen': {'seed' : 1}, 'iterations': mcmc_iterations, 'after_minimizer': {'type': 'root_minuit'}})
            if not always_mcmc: minimizers.append({'type': 'root_minuit', 'method': 'simplex'})
            minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min1', 'rnd_gen': {'seed' : 1}, 'iterations': mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
            if strategy == 'robust':
                minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min2', 'rnd_gen': {'seed' : 1}, 'iterations': 50 * mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
                minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min3', 'rnd_gen': {'seed' : 1}, 'iterations': 500 * mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
            bootstrap_mcmcpars = options.getint('minimizer', 'bootstrap_mcmcpars')
            if bootstrap_mcmcpars > 0:
                for m in minimizers:
                    if m['type'] == 'mcmc_minimizer': m['bootstrap_mcmcpars'] = bootstrap_mcmcpars
            result = {'type': 'minimizer_chain', 'minimizers': minimizers}
            if self.need_error: result['last_minimizer'] = {'type': 'root_minuit'}
            
        def apply_tolfactor(d):
            if type(d) in (list, tuple):
                for v in d: apply_tolfactor(v)
                return
            if type(d)==dict:
                if 'type' in d and d['type'] == 'root_minuit': d['tolerance_factor'] = tolerance_factor
                for key, value in d.iteritems():
                    apply_tolfactor(value)
        if tolerance_factor != 1: apply_tolfactor(result)
        return result


def _settingvalue_to_cfg(value, indent=0, current_path = [], value_is_toplevel_dict = False):
    if type(value) == numpy.float64: value = float(value)
    tval = type(value)
    if tval == array.array:
        return "(" + ",".join(map(lambda f: "%.5e" % f, value)) + ")"
    if tval == str: return '"%s"' % value
    if tval == bool: return 'true' if value else 'false'
    if tval == int: return '%d' % value
    if tval == float:
        if value == inf: return '"inf"'
        elif value == -inf: return '"-inf"'
        return '%.5e' % value
    if tval == list or tval == tuple:
        return "(" + ",".join([_settingvalue_to_cfg(value[i], indent + 4, current_path + [str(i)]) for i in range(len(value))]) + ')'
    if tval == dict:
        result = ''
        if not value_is_toplevel_dict: result += "{\n"
        new_indent = (indent + 4) if not value_is_toplevel_dict else 0
        # sort keys to make config files reproducible:
        for s in sorted(value.keys()):
            result += ' ' * new_indent + s + " = " + _settingvalue_to_cfg(value[s], new_indent, current_path + [s]) + ";\n"
        if not value_is_toplevel_dict: result += ' ' * indent + "}"
        return result
    raise RuntimeError, "Cannot convert type %s to theta cfg in path '%s'" % (type(value), '.'.join(current_path))




# return a config dictionary for the signal prior. s can be a string
# such as "flat", "fix:XXX", etc., or a dictionary in which case
# it is returned unmodified
#
# None is equivalent to 'flat'
def _signal_prior_dict(spec):
    if spec is None: spec = 'flat'
    if type(spec) == str:
        if spec.startswith('flat'):
	    value = None
            if spec.startswith('flat:'):
                res = re.match('flat:\[([^,]+),(.*)\](:.*)?', spec)
                if res is None: raise RuntimeError, "signal_prior specification '%s' invalid (does not match range specification syntax)" % spec
                xmin, xmax = float(res.group(1)), float(res.group(2))
                svalue = res.group(3)
		if svalue is not None: value = float(svalue[1:])
            else:
                if spec!='flat': raise RuntimeError, "signal_prior specification '%s' invalid" % spec
                xmin, xmax = 0.0, float("inf")
            if value is None: value = 0.5 * (xmax - xmin)
            if value==float("inf"): value = 1.0
            signal_prior_dict = {'type': 'flat_distribution', 'beta_signal': {'range': [xmin, xmax], 'fix-sample-value': value}}
        elif spec.startswith('fix:'):
            v = float(spec[4:])
            signal_prior_dict = {'type': 'delta_distribution', 'beta_signal': v}
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if type(spec) != dict: raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a dictionary!"
        signal_prior_dict = spec
    return signal_prior_dict

# return Distribution object instead of dictionary
def _signal_prior_dist(spec):
    if spec is None: spec = 'flat'
    dist = Distribution()
    if spec.startswith('flat'):
	value = None
        if spec.startswith('flat:'):
            res = re.match('flat:\[([^,]+),(.*)\](:.*)?', spec)
            if res is None: raise RuntimeError, "signal_prior specification '%s' invalid (does not match range specification syntax)" % spec
            xmin, xmax = float(res.group(1)), float(res.group(2))
            svalue = res.group(3)
            if svalue is not None: value = float(svalue[1:])
        else:
            if spec!='flat': raise RuntimeError, "signal_prior specification '%s' invalid" % spec
            xmin, xmax = 0.0, float("inf")
        if value is None: value = 0.5 * (xmax - xmin)
        if value==float("inf"): value = 1.0
        dist.set_distribution('beta_signal', 'gauss', value, width = inf, range = [xmin, xmax])
    elif spec.startswith('fix:'):
        v = float(spec[4:])
        dist.set_distribution('beta_signal', 'gauss', mean = v, width = 0.0, range = [v,v])
    return dist


class DbResult(object):
    def __init__(self, theta_db_fname = None):
        self._fname = theta_db_fname
        self._init_conn()
        
    def set_db_fname(self, fname):
        self._fname = fname
        self._init_conn()

    def _init_conn(self):
        if self._fname is None: return
        if not os.path.exists(self._fname): raise RuntimeError, "The file '%s' does not exist!" % self._fname
        self.conn = sqlite3.connect(self._fname)
    
    # returns a cursor
    def _query(self, sql):
        c = self.conn.cursor()
        try: c.execute(sql)
        except Exception, ex:
            print "exception executing %s on file %s: %s" % (sql, self._fname, str(ex))
            raise ex
        return c
        
    def _get_tables(self):
        c = self._query("select name from sqlite_master where type='table'")
        rows = c.fetchall()
        c.close()
        table_names = [r[0] for r in rows]
        return table_names
    
    
    def _get_columns(self, table):
        c = self._query("PRAGMA table_info(\"%s\")" % table)
        columns = c.fetchall()
        colnames = [c[1] for c in columns]
        return colnames

    def get_results(self, table_name, columns, order_by = None):
        """
        Get the content of the given table
        
        Parameters:
        
        * ``table_name`` is the name of the table
        * ``columns`` should either be a list of column names or the string "*" to select all available columns.
        * ``order_by`` is an optional column name which will be used to sort the result.
        
        The return value is a dictionary. The keys are the column names; each dictionary value is a list of values
        for this column.
        """
        if self._fname is None: return None
        all_tables = self._get_tables()
        tables = [table_name] + [t for t in all_tables if t.startswith(table_name + '__')]
        order_by_sql = ''
        if len(tables) == 1:
            if type(columns)!=str:
                # make sure all columns asked for exist (if we don't do this, then sqlite returns the column name in case a
                # column does not exist, which is annoying and a hard to detect bug).
                cols = self._get_columns(tables[0])
                for c in columns:
                    if c not in cols: raise RuntimeError, "Asked for column '%s' but it is not present in table '%s' (columns: %s)" % (c, tables[0], cols)
                columns_sql = '"' + '", "'.join(columns) + '"'
            else:
                assert columns == '*'
                columns_sql = '*'
            if order_by is not None:
                order_by_sql = ' order by "%s"' % order_by
        else:
            # always select all columns in this case, as otherwise, we would need to make a mapping column names <-> tables,
            # which is probably not worth the effort. The columns will be filtered below.
            columns_sql = '*'
            # order_by is not supported in this case:
            if order_by: raise RuntimeError, "order_by not supported for large tables (specified order_by='%s' for table='%s' matching actual sqlite tables %s)" %(order_by, table_name, str(tables))
        result = {}
        for table_name in tables:
            c = self._query('select %s from "%s"%s' % (columns_sql, table_name, order_by_sql))
            data = c.fetchall()
            res_columns = c.description
            c.close()
            for i in range(len(res_columns)):
                colname = res_columns[i][0].strip('"')
                if columns != '*' and colname not in columns:
                    continue
                result[colname] = [row[i] for row in data]
        return result
        
    def get_products(self, columns = '*', fail_if_empty = True):
        """
        Get the result (in the products table) of this invocation of theta as dictionary::
        
           column_name -> list of values
           
        Which columns are returned can be controlled by the ``columns`` argument which is either a list of column names
        or as special case the string '*'. The default is '*' which will return all columns.
        
        Which columns are available depends on the producers :program:`theta` has run. Refer to the doxygen documentation
        of the producers for the list of columns.
    
        If ``fail_if_empty`` is ``True`` and no rows are found, some information from the log table is printed and a `RuntimeError` is raised.
        This can happen e.g. in case there was one attempt to fit on data which failed due to a minimizer which failed.
        """
        res = self.get_results('products', columns)
        if fail_if_empty and len(res[res.keys()[0]])==0:
            self.print_logtable_errors()
            raise RuntimeError, "No result is available, see errors above"
        return res
        
    def print_logtable_errors(self, limit = 10):
        if not self._fname: return
        data = sql_singlefile(self._fname, 'select "eventid", "message" from "log" where "severity"=0 limit %d' % int(limit))
        if len(data)==0: return
        print "There have been errors for some toys: eventid, message"
        for i in range(len(data)):
            print data[i][0], data[i][1]
        
    def get_products_nevents(self):
        data = sql_singlefile(self._fname, 'select count(*) from products')
        return data[0][0]
    
    def get_db_fname(self):
        """
        Return the path to the db file created by theta. Is ``None`` if no db file is available (yet).
        """
        return self._fname


class MainBase(ModuleBase, DbResult):
    """
    Base class corresponding to a theta `Main` plugin. One instance of this class corresponds to
    one theta .cfg file and one execution of theta, and one db result file. This class contains management
    code for config file creation, theta execution, and (via `DbResult`) reading the resulting db file. It
    also takes care of the theta-auto cache management.
    
    Derived classes must implement a method `get_cfg(self, options)`.
    """
    
    def __init__(self, fname_fragment):
        """
        `name` is used in the created cnofiguration file name
        """
        ModuleBase.__init__(self)
        DbResult.__init__(self)
        self.result_available = False
        # full path + filename of theta configuration file, if created:
        self.theta_cfg_fname = None
        self.thread = None
        self.cfg_fname_fragment = fname_fragment
        
    
    # return a dictionary with settings for the theta configuration for the plugins, parameters, rvobservables, and observables.
    def common_main_cfg(self, model, signal_processes, options):
        plugins = self.get_required_plugins(options)
        result = {'options': {'plugin_files': ['$THETA_DIR/lib/'+s for s in sorted(list(plugins))]}}
        result['parameters'] = sorted(list(model.get_parameters(signal_processes)))
        rvobservables = sorted(list(model.rvobs_distribution.get_parameters()))
        if len(rvobservables) > 0:
            result['rvobservables'] = rvobservables
        result['observables'] = {}
        for obs in sorted(list(model.get_observables())):
            xmin, xmax, nbins = model.observables[obs]
            result['observables'][obs] = {'range': [xmin, xmax], 'nbins': nbins}
        return result
    
    
    def _check_cache(self, options):
        """
        Check in the cache directory whether result is already available there
        returns True if it is, false otherwise.
        Note that get_result_available will return `True` if the result is in the cache
        only after calling check_cache.
        """
        workdir = options.get_workdir()
        cache_dir = os.path.join(workdir, 'cache')
        cfgfile_full = self.get_configfile(options)
        cfgfile_cache = os.path.join(cache_dir, os.path.basename(cfgfile_full))
        dbfile_cache = cfgfile_cache.replace('.cfg', '.db')
        if os.path.exists(cfgfile_cache) and os.path.exists(dbfile_cache) and  open(cfgfile_cache, 'r').read() == open(cfgfile_full, 'r').read():
            self.set_db_fname(dbfile_cache)
            self.result_available = True
            return True
        return False
    
    
    def get_configfile(self, options):
        """
        Returns the path to the :program:`theta` configuration file, creating it if necessary.
        """
        if self.theta_cfg_fname is not None: return self.theta_cfg_fname
        cfg_dict = self.get_cfg(options)
        cfg_string = _settingvalue_to_cfg(cfg_dict, value_is_toplevel_dict = True)            
        output_name =  self.cfg_fname_fragment + '-' + hashlib.md5(cfg_string).hexdigest()[:10]
        if '@output_name' in cfg_string: cfg_string += "output_name = \"%s.db\";\n" % output_name
        if '@debuglog_name' in cfg_string: cfg_string += 'debuglog_name = "%s-debuglog.txt";\n' % output_name
        workdir = options.get_workdir()
        assert os.path.exists(workdir)
        self.theta_cfg_fname = os.path.join(workdir, output_name + '.cfg')
        f = open(self.theta_cfg_fname, 'w')
        f.write(cfg_string)
        f.close()
        return self.theta_cfg_fname
    

    # can be run in a different thread
    def _exec(self, cmd):
        self.error = None
        ret = os.system(cmd)
        if ret != 0:
            self.error = "executing theta ('%s') failed with exit code %d" % (cmd, ret)
            if os.isatty(1):
                attr = termios.tcgetattr(1)
                attr[3] |= termios.ECHO
                termios.tcsetattr(1, termios.TCSANOW, attr)
                        
    # should always be run in the main thread, after _exec() returns
    def _cleanup_exec(self):
        cfgfile = self.theta_cfg_fname
        assert cfgfile is not None
        dbfile = os.path.basename(cfgfile).replace('.cfg', '.db')
        if self.error is None:
            cache_dir = os.path.join(os.path.dirname(cfgfile), 'cache')
            if not os.path.exists(cache_dir): os.mkdir(cache_dir)
            cfgfile_cache = os.path.join(cache_dir, os.path.basename(cfgfile))
            dbfile_cache = cfgfile_cache.replace('.cfg', '.db')
            shutil.move(dbfile, dbfile_cache)
            shutil.copy(cfgfile, cfgfile_cache)
            self.result_available = True
            self.set_db_fname(dbfile_cache)
        else:
            if os.path.exists(dbfile) and not self.debug: os.unlink(dbfile)
            error = self.error
            self.error = None
            raise RuntimeError, error
            
    def run_theta(self, options, in_background_thread = False):
        """
        Execute the :program:`theta` program on the current .cfg file locally, creating the config file if necessary.
        If the database file is already in the cache -- and cache checking is not diabled via ``options`` -- :program:`theta`
        is not actually executed.
        
        Parameters:
        
        * ``options`` - an instance of :class:`Options` which is passed to the submodules to control some details of the configuration
        * ``in_background_thread`` - if ``True``, this methods returns immediately and :program:`theta` is executed in a background thread. Otherwise, the method blocks until
          :program:`theta` execution is done.
        
        Communicating errors in :program:`theta` execution are communicated via a RuntimeException.
        The point at which it is raised depends on the value of ``in_background_thread``:
        
         * if ``in_background_thread = False``, this method will raise the exception
         * if ``in_background_thread = True``, this method returns immediately and cannot know about any error in the :program:`theta` execution. Instead, the eception will be raised in
           :meth:`wait_for_result_available`.
                 
        Executing multiple instances of the theta program in parallel is usually done like this::
           
           runs = [] # the list of Run instances
           ...
           for r in runs: r.run_theta(options, in_background_thread = True)
           for r inruns: r.wait_for_result_available()
        
        The first for-loop spawns the threads which execute theta, so theta is executed in parallel. The second loop waits for all theta executions
        to terminate (either normally or abnormally).
        """
        check_cache = options.getboolean('global', 'check_cache')
        if check_cache: self._check_cache(options)
        cfgfile = self.get_configfile(options)
        self.debug = options.getboolean('global', 'debug')
        if self.result_available:
            if self.debug: info("Found config file '%s' in cache, not running again." % cfgfile)
            return
        workdir = options.get_workdir()
        cache_dir = os.path.join(workdir, 'cache')
        theta = os.path.realpath(os.path.join(config.theta_dir, 'bin', 'theta'))
        cmd = theta + " " + cfgfile
        if self.debug: cmd += " --print-time"
        else: cmd += " --nowarn"
        info("Running '%s'" % cmd)
        to_execute = lambda : self._exec(cmd)
        if in_background_thread:
            assert self.thread is None
            self.thread = threading.Thread(target = to_execute)
            self.thread.start()
        else:
            to_execute()
            self._cleanup_exec()    
    
    def wait_for_result_available(self):
        """
        Waits untils :program:`theta` execution terminates.
        
        This method is only useful after calling :meth:`run_theta` with ``in_background_thread = True``.
        
        Throws and exception in case :program:`theta`  exited with an error.
        """
        if self.result_available: return
        if self.thread is None: return
        self.thread.join()
        self.thread = None
        self._cleanup_exec()
    
class Run(MainBase):
    def __init__(self, model, signal_processes, signal_prior, input, n, producers, nuisance_prior_toys = None, seed = None):
        """
        Note that you usually do not need to construct instances of this class directly. Rather, construct them via
        methods such as :meth:`deltanll` using the parameter ``run_theta = False`` which will give you a list of properly
        configured :class:`Run` instances.
        
        Parameters:
        
         * ``model`` - the statistical model, a :class:`Model` instance
         * ``signal_processes`` - a list of strings, the processes to be considered as signal: they will be scaled together with the parameter ``beta_signal``
         * ``producers`` - a list of producers to run, of class :class:`ProducerBase`. This corresponds to C++ plugins of type ``Producer`` in theta.
         * ``nuisance_prior_toy`` - the nuisance parameter prior to use for toy data generation. Either a :class:`Distribution` instance or ``None`` is use ``model.distribution``
         * ``seed`` - the random seed for toy data generation. Set to ``None`` to use a different seed for each :program:`theta` execution
         
        For ``signal_prior``, ``input``, ``n``, see :ref:`common_parameters`.
        """
        input_fragment = str(input)
        if '/' in input_fragment: input_fragment = input_fragment[input_fragment.rfind('/')+1:]
        fragment = ''.join([p.name for p in producers]) + '-' + input_fragment + '-' + ''.join(signal_processes)
        MainBase.__init__(self, fragment)
        self.model = model
        self.signal_processes = signal_processes
        self.n_events = n
        self.producers = producers
        self.nuisance_prior_toys = nuisance_prior_toys
        for p in producers: self.add_submodule(p)
        self.seed = seed
        if type(input) == str:
            self.data_source = DataSource(input, model, signal_processes, override_distribution = nuisance_prior_toys, seed = seed)
        else:
            assert isinstance(input, DataSource), "unexpected type for 'input' argument: %s" % type(input)
            self.data_source = input
        self.signal_prior_cfg = _signal_prior_dict(signal_prior)
        
    def get_cfg(self, options):
        result = self.common_main_cfg(self.model, self.signal_processes, options)
        main = {'n-events': self.n_events, 'producers': [], 'log-report': False, 'output_database': {'type': 'sqlite_database', 'filename': '@output_name'},
            'model': "@model", 'data_source': self.data_source.get_cfg(options)}
        n_threads = options.getint('main', 'n_threads')
        if n_threads > 1:
            main['type'] = 'run_mt'
            main['n_threads'] = n_threads
        for p in self.producers:
            result['p%s' % p.name] = p.get_cfg(options)
            main['producers'].append('@p%s' % p.name)
        result.update({'main': main, 'model': self.model.get_cfg(self.signal_processes, self.signal_prior_cfg, options)})
        robust_nll = options.getboolean('model', 'robust_nll')
        if robust_nll: result['model']['robust_nll'] = True
        use_llvm = options.getboolean('model', 'use_llvm')
        use_tbb = options.getboolean('model', 'use_tbb')
        if use_tbb and use_llvm:
            raise RuntimeError, "specified both options: model use_llvm and model use_tbb. This is not allowed."
        if use_llvm:
            print "Using llvm. This is EXPERIMENTAL. Use at your own risk"
            result['model']['type'] = 'llvm_model'
            result['options']['plugin_files'].append('$THETA_DIR/lib/llvm-plugins.so')
        elif use_tbb:
            print "Using tbb. This is EXPERIMENTAL. Use at your own risk"
            result['model']['type'] = 'tbb_model'
            result['options']['plugin_files'].append('$THETA_DIR/lib/tbb-plugins.so')
            n_threads = options.getint('model', 'tbb_nthreads')
            if n_threads > 0:
                result['model']['n_threads'] = n_threads
        return result

class ClsMain(MainBase):
    def __init__(self, model, signal_processes, signal_prior, input, producers, nuisance_prior_toys = None, seed = None):
        """
        Note that `input` is used for the "observed" cls limit; only the first (pseudo-)dataset of this DataSource will be used.
        """
        MainBase.__init__(self, 'cls')
        assert len(producers)==1
        self.model = model
        self.signal_processes = signal_processes
        self.producers = producers
        self.nuisance_prior_toys = nuisance_prior_toys
        for p in producers: self.add_submodule(p)
        self.seed = seed
        if input is None: self.data_source = None
        else:
            if type(input) == str:
                self.data_source = DataSource(input, model, signal_processes, override_distribution = nuisance_prior_toys, seed = seed)
            else:
                assert type(input) in (RootNtupleSource, DataSource), "unexpected type for 'input' argument: %s" % type(input)
                self.data_source = input
        self.signal_prior_cfg = _signal_prior_dict(signal_prior)
        self.cls_options = {'ts_column': self.producers[0].name + '__nll_diff'}
        self.minimizer = Minimizer(need_error = True, need_covariance = False)
        self.add_submodule(self.minimizer)
        
    def set_cls_options(self, **args):
        allowed_keys = ['ts_column', 'expected_bands', 'clb_cutoff', 'tol_cls', 'truth_max', 'reltol_limit', 'frequentist_bootstrapping', 'input_expected']
        for k in allowed_keys:
            if k in args:
                self.cls_options[k] = args[k]
                del args[k]
        if len(args) > 0: raise RuntimeError, "unrecognized cls options: %s" % str(args.keys())
    
    def get_cfg(self, options):
        assert 'ts_column' in self.cls_options
        result = self.common_main_cfg(self.model, self.signal_processes, options)
        seed = self.seed
        if seed is None: seed = -1
        main = {'type': 'cls_limits', 'producer': self.producers[0].get_cfg(options), 'output_database': {'type': 'sqlite_database', 'filename': '@output_name'},
            'truth_parameter': 'beta_signal', 'tol_cls': self.cls_options.get('tol_cls', 0.025),
            'clb_cutoff': self.cls_options.get('clb_cutoff', 0.02), 'model': '@model',
            'debuglog': '@debuglog_name', 'rnd_gen': {'seed': seed }, 'ts_column': self.cls_options['ts_column'],
            'minimizer': self.minimizer.get_cfg(options), 'expected_bands': self.cls_options.get('expected_bands', 2000)}
        if not options.getboolean('cls_limits', 'write_debuglog'): del main['debuglog']
        if self.cls_options.get('frequentist_bootstrapping', False):  main['nuisancevalues-for-toys'] = 'datafit'
        if 'truth_max' in self.cls_options: main['truth_max'] = float(self.cls_options['truth_max'])
        if 'reltol_limit' in self.cls_options: main['reltol_limit'] = float(self.cls_options['reltol_limit'])
        if self.data_source is not None: main['data_source'] = self.data_source.get_cfg(options)
        if 'input_expected' in self.cls_options:
            print "using input = %s for expected limits" % self.cls_options['input_expected']
            ds = DataSource(self.cls_options['input_expected'], self.model, self.signal_processes,
                override_distribution = self.nuisance_prior_toys)
            main['data_source_expected'] = ds.get_cfg(options)
        result.update({'main': main, 'model': self.model.get_cfg(self.signal_processes, self.signal_prior_cfg, options)})
        return result
    
    
class AsymptoticClsMain(MainBase):
    def __init__(self, model, signal_processes, input, n = 1, beta_signal_expected = 0.0):
        """
        `input` should be either 'data' or None.
        """
        MainBase.__init__(self, 'acls')
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)
        self.model = model
        self.n = n
        self.signal_processes = signal_processes
        self.beta_signal_expected = beta_signal_expected
        if input is None: self.input = None
        else:
            self.input = DataSource(input, self.model, self.signal_processes)
        self.signal_prior_cfg = _signal_prior_dict('flat')
        
    def get_cfg(self, options):
        result = self.common_main_cfg(self.model, self.signal_processes, options)
        main = {'type': 'asymptotic_cls', 'model': '@model', 'parameter': 'beta_signal', 'minimizer': self.minimizer.get_cfg(options),
            'output_database': {'type': 'sqlite_database', 'filename': '@output_name'}}
        if self.input is not None:
            main['data'] = self.input.get_cfg(options)
            if self.n != 1: main['n'] = self.n
        if self.beta_signal_expected is None:
            main['quantiles_expected'] = []
        elif self.beta_signal_expected != 0.0:
            main['parameter_value_expected'] = self.beta_signal_expected
        result.update({'main': main, 'model': self.model.get_cfg(self.signal_processes, self.signal_prior_cfg, options)})
        return result


# create a Histogram instance from the blob data in a sqlite db
def histogram_from_dbblob(blob_data, blob_uncertainties = None):
    a = array.array('d')
    a.fromstring(blob_data)
    uncs = None
    if blob_uncertainties is not None:
        a_uncs = array.array('d')
        a_uncs.fromstring(blob_uncertainties)
        uncs = a_uncs[3:-1]
    return Histogram(xmin = a[0], xmax = a[1], values = a[3:-1], uncertainties = uncs)

def matrix_from_dbblob(blob_data):
    n = int(math.sqrt(len(blob_data) / 8 - 4) + 0.4)
    assert (n*n + 4) * 8 == len(blob_data)
    return numpy.ndarray(shape = (n, n), buffer = blob_data[3*8:-8])
