# -*- coding: utf-8 -*-
import ROOT, re, fnmatch, math
import os, os.path
from theta_interface import *
from utils import *


class rootfile:
    
    # these are caching dictionaries indexing by filename
    tfiles = {}
    all_templates = {}
    
    def __init__(self, filename):
        assert os.path.isfile(filename), "File %s not found (cwd: %s)" % (filename, os.getcwd())
        self.filename = filename
        if self.filename not in rootfile.tfiles: rootfile.tfiles[self.filename] = ROOT.TFile(filename, "read")
        self.tfile = rootfile.tfiles[self.filename]
        
    @staticmethod
    def th1_to_histo(th1, include_uncertainties):
        xmin, xmax, nbins = th1.GetXaxis().GetXmin(), th1.GetXaxis().GetXmax(), th1.GetNbinsX()
        #data = [th1.GetBinContent(i) for i in range(1, nbins+1)]
        values = array.array('d', [th1.GetBinContent(i) for i in range(1, nbins+1)])
        uncertainties = array.array('d', [th1.GetBinError(i) for i in range(1, nbins+1)]) if include_uncertainties else None
        return Histogram(xmin, xmax, values, uncertainties, th1.GetName())
  
    # get all templates as dictionary (histogram name) --> (float xmin, float xmax, list float data)
    # only checks type of histogram, not naming convention
    def get_all_templates(self, include_uncertainties, warn = True):
        if (self.filename, include_uncertainties) in rootfile.all_templates: return rootfile.all_templates[self.filename, include_uncertainties]
        result = {}
        l = self.tfile.GetListOfKeys()
        for key in l:
            clas = key.GetClassName()
            if clas == 'TDirectoryFile': continue
            if clas not in ('TH1F', 'TH1D') and warn:
                print "WARNING: ignoring key %s in input file %s because it is of ROOT class %s, not TH1F / TH1D" % (key.GetName(), self.filename, clas)
                continue
            th1 = key.ReadObj()
            result[str(key.GetName())] = rootfile.th1_to_histo(th1, include_uncertainties)
        rootfile.all_templates[self.filename, include_uncertainties] = result
        return result
        
    def get_histogram(self, hname, include_uncertainties, fail_with_exception = False):
        if (self.filename, include_uncertainties) not in rootfile.all_templates: self.get_all_templates(include_uncertainties, False)
        if hname not in rootfile.all_templates[(self.filename, include_uncertainties)]:
            if fail_with_exception: raise RuntimeError, "histogram '%s' in root file '%s' not found!"
            else: return None
        return rootfile.all_templates[(self.filename, include_uncertainties)][hname].copy()


# par_values is a dictionary (parameter name) --> (floating point value)
#
# If "beta_signal" is in par_values, it will be included as factor for all signal processes. Otherwise,
# it is not mutliplied, effectively assuming beta_signal = 1.0.
#
# returns a dictionary
# (observable name) --> (process name) --> (xmin, xmax, data)
def get_shifted_templates(model, par_values, include_signal = True, observables = None, eval2 = False):
    result = {}
    if observables is None: observables = model.get_observables()
    for obs in observables:
        result[obs] = {}
        for p in model.get_processes(obs):
            if p in model.signal_processes and not include_signal: continue
            coeff = model.get_coeff(obs, p).get_value(par_values)
            if p in model.signal_processes: coeff *= par_values['beta_signal']
            if eval2:
                histo = model.get_histogram_function(obs, p).evaluate2(par_values).scale(coeff)
            else:
                histo = model.get_histogram_function(obs, p).evaluate(par_values)
                for i in range(len(histo[2])):
                    histo[2][i] *= coeff
            result[obs][p] = histo
    return result


## \brief A statistical Model
#
# This class contains all information about a statistical model, including:
# * data histograms (if any)
# * additional likelihood terms
# * observables, with xmin, xmax, ranges
# * which processes are potentially signal
#
# In general, more than one process is marked as 'signal' (such as different masses), but only
# a subset of these, the signal process group, is actually used to build a concrete, usable model. In this
# sense, the Model class represents a collection of <em>potential</em> models. In order to write a concrete
# model configuration, some additional informations are required, namely:
# * which processes are to be considered as signal; these will all be scaled by the 'beta_signal' parameter
# * the prior distribution for beta_signal for model.parameter_distribution
#
#
# notes:
# * the coefficients should be chosen such that their mean / most probable value is 1.0 as this is assumed e.g. in model_summary
class Model:
    def __init__(self):
        # observables is a dictionary str name ->  (float xmin, float xmax, int nbins)
        self.observables = {}
        self.processes = set() # all processe
        self.signal_processes = set() # the subset of processes considered signal
        # in addition to the signal_processes, we need to define which processes to consider together
        # as signal. This is defined in signal_process_groups. It is a dictionary from the signal process group is to
        # the list of signal processes to consider.
        self.signal_process_groups = {}
        # observable_to_pred is a dictionary from (str observable) --> dictionary (str process)
        #  --> | 'histogram'            --> HistogramFunction instance
        #      | 'coefficient-function' --> Function instance
        self.observable_to_pred = {}
        # like model.parameter-distribution (excluding beta_signal):
        self.distribution = Distribution()
        # data histograms: dictionary str obs -> histo
        self.data_histos = {}
        # real-valued observable values: dictionary str -> float
        self.data_rvobsvalues = {}
        # a FunctionBase instance or None:
        self.additional_nll_term = None
        self.bb_uncertainties = False
        self.rvobs_distribution = Distribution()
    
    def reset_binning(self, obs, xmin, xmax, nbins):
        assert obs in self.observables
        self.observables[obs] = (xmin, xmax, nbins)
    
    # combines the current with the other_model by including the observables and processes from the other_model.
    # Note that:
    # * other_model must not include an observable of the current model
    # * other_model should define the same signal_processes (this is enforced if strict=true)
    # * other_model.distribution can include the same nuisance parameters (=same name).
    #   For shared nuisance parameters, the prior for self.distribution and other.distribution must be identical.
    def combine(self, other_model, strict=True):
        assert len(set(self.observables.keys()).intersection(set(other_model.observables.keys())))==0, "models to be combined share observables, but they must not!"
        if strict:
            assert self.signal_processes == other_model.signal_processes, "signal processes not equal: left-right=%s; right-left=%s;" \
                       % (str(self.signal_processes.difference(other_model.signal_processes)), str(other_model.signal_processes.difference(self.signal_processes)))
        self.distribution = Distribution.merge(self.distribution, other_model.distribution, False)
        self.rvobs_distribution = Distribution.merge(self.rvobs_distribution, other_model.rvobs_distribution, False)
        self.observables.update(other_model.observables)
        self.processes.update(other_model.processes)
        self.data_histos.update(other_model.data_histos)
        self.observable_to_pred.update(other_model.observable_to_pred)
        self.bb_uncertainties = self.bb_uncertainties or other_model.bb_uncertainties
        if other_model.additional_nll_term is not None:
            if self.additional_nll_term is None: self.additional_nll_term = other_model.additional_nll_term
            else: self.additional_nll_term = self.additional_nll_term + other_model.additional_nll_term

    def rename_parameter(self, current_name, new_name):
        self.distribution.rename_parameter(current_name, new_name)
        for o in self.observable_to_pred:
            for p in self.observable_to_pred[o]:
                self.observable_to_pred[o][p]['histogram'].rename_parameter(current_name, new_name)
                self.observable_to_pred[o][p]['coefficient-function'].rename_parameter(current_name, new_name)
    
    # modify the model to only contain a subset of the current observables.
    # The parameter observables must be convertible to a set of strings
    def restrict_to_observables(self, observables):
        observables = set(observables)
        model_observables = set(self.observables.keys())
        assert observables.issubset(model_observables), "observables must be a subset of model.observables!"
        for obs in model_observables:
            if obs in observables: continue
            del self.observables[obs]
            del self.observable_to_pred[obs]
            if obs in self.data_histos: del self.data_histos[obs]
        # in theory, we could also update self.processes / self.signal_processes and
        # self.distribution (if the model now has fewer dependencies), but this is not necessary.
       
    def set_data_histogram(self, obsname, histo, reset_binning = False):
        xmin, xmax, nbins = histo[0], histo[1], len(histo[2])
        if obsname not in self.observables:
            self.observables[obsname] = xmin, xmax, nbins
            self.observable_to_pred[obsname] = {}
        if reset_binning:
            self.observables[obsname] = xmin, xmax, nbins
        xmin2, xmax2, nbins2 = self.observables[obsname]
        assert (xmin, xmax, nbins) == (xmin2, xmax2, nbins2)
        self.data_histos[obsname] = histo

    def get_data_histogram(self, obsname):
        if obsname in self.data_histos: return self.data_histos[obsname]
        else: return None
    
    def get_data_rvobsvalues(self):
        return self.data_rvobsvalues

    def get_observables(self):
        return set(self.observables.keys())
    
    def has_data(self):
        for o in self.observables:
            if o not in self.data_histos: return False
        return True
    
    #important: always call set_histogram_function first, then get_coeff!
    def set_histogram_function(self, obsname, procname, histo):
        xmin, xmax, nbins = histo.nominal_histo[0], histo.nominal_histo[1], len(histo.nominal_histo[2])
        if obsname not in self.observables:
            self.observables[obsname] = xmin, xmax, nbins
            self.observable_to_pred[obsname] = {}
        xmin2, xmax2, nbins2 = self.observables[obsname]
        assert (xmin, xmax, nbins) == (xmin2, xmax2, nbins2), "detected inconsistent binning setting histogram for (obs, proc) = (%s, %s)" % (obsname, procname)
        self.processes.add(procname)
        if procname not in self.observable_to_pred[obsname]: self.observable_to_pred[obsname][procname] = {}
        self.observable_to_pred[obsname][procname]['histogram'] = histo
        self.observable_to_pred[obsname][procname]['coefficient-function'] = Function()
        
    def get_coeff(self, obsname, procname):
        return self.observable_to_pred[obsname][procname]['coefficient-function']
        
    def get_processes(self, obsname):
        return self.observable_to_pred[obsname].keys()
        
    def get_histogram_function(self, obsname, procname):
        if procname not in self.observable_to_pred[obsname]: return None
        return self.observable_to_pred[obsname][procname]['histogram']
    
    # make sure all histogram bins of histogram h contain at least
    # epsilon * (average bin content of histogram h)
    #
    # The special case epsion = None will leave the prediction at the value of 0.0
    # but will make sure that the uncertainty in each bin is at least the one corresponding 1 MC event.
    # The average weight for the MC sample is estimated from the uncertainties given in the histogram.
    def fill_histogram_zerobins(self, epsilon = 0.001):
        for o in self.observable_to_pred:
            for p in self.observable_to_pred[o]:
                hf = self.observable_to_pred[o][p]['histogram']
                if epsilon is None:
                    h = hf.get_nominal_histo()
                    unc = h.get_value_sum_uncertainty()
                    if unc is None: continue
                    s = h.get_value_sum()
                    # the number of effective MC events is estimated via the MC stat. uncertainty:
                    n_eff = s**2 / unc**2
                    # the average MC event weight; this is used as minimum uncertainty ...
                    w = s / n_eff
                    #print o, p, n_eff, w
                    uncs = h.get_uncertainties()
                    new_uncs = [max(u, w) for u in uncs]
                    hf.set_nominal_histo(Histogram(h.get_xmin(), h.get_xmax(), h.get_values(), new_uncs, h.get_name()))
                else:
                    histos = [hf.get_nominal_histo()]
                    for par in hf.syst_histos: histos.extend([hf.syst_histos[par][0], hf.syst_histos[par][1]])
                    for h in histos:
                        s = sum(h[2])
                        nbins = len(h[2])
                        h[2][:] = array.array('d', [max(epsilon * s / nbins, y) for y in h[2]])

    # obsname and procname can be '*' to match all observables / processes.
    def scale_predictions(self, factor, procname = '*', obsname = '*'):
        found_match = False
        for o in self.observable_to_pred:
            if obsname != '*' and o!=obsname: continue
            for p in self.observable_to_pred[o]:
                if procname != '*' and procname != p: continue
                found_match = True
                self.observable_to_pred[o][p]['histogram'].scale(factor)
        if not found_match: raise RuntimeError, 'did not find obname, procname = %s, %s' % (obsname, procname)

    def rebin(self, obsname, rebin_factor):
        xmin, xmax, nbins_old = self.observables[obsname]
        assert nbins_old % rebin_factor == 0
        self.observables[obsname] = (xmin, xmax, nbins_old / rebin_factor)
        pred = self.observable_to_pred[obsname]
        for p in pred:
            pred[p]['histogram'].rebin(rebin_factor)
        if obsname in self.data_histos:
            rebin_hlist(self.data_histos[obsname][2], rebin_factor)
    
    ## \brief define which processes should be considered as signal processes
    #
    # Any process not defined explicitly as signal is considered to be background.
    # Assumes that each signal process is independent, and sets the signal_process_groups accordingly.
    # For more control in situations where the signal has multiple histograms to be used simultaneously,
    # use the set_signal_process_groups.
    #
    # procs is a list / set of glob patterns (or a single pattern).
    def set_signal_processes(self, procs):
        if type(procs)==str: procs = [procs]
        self.signal_processes = set()
        for pattern in procs:
            found_match = False
            for p in self.processes:
                if fnmatch.fnmatch(p, pattern):
                    found_match = True
                    self.signal_processes.add(p)
            if not found_match: raise RuntimeError, "no match found for pattern '%s'" % pattern
        self.signal_process_groups = {}
        for p in self.signal_processes: self.signal_process_groups[p] = [p]
        
    ## \brief Define the signal process groups
    #
    # groups is a dictionary (id) --> (list of processes)
    def set_signal_process_groups(self, groups):
        for spid in groups:
            for p in groups[p]: assert type(p)==str and p in self.processes, "unknown process '%s'" % p
        self.signal_processes = set()
        for spid in groups:
            self.signal_processes.update(groups[spid])
        self.signal_process_groups = copy.deepcopy(groups)

    ## \brief Add a new, rate-only uncertainty to the model
    #
    # Adds a new parameter with name \c u_name tio the distribution (unless it already exists) with
    # a Gauss prior around 0.0 with width 1.0 and adds a factor exp(u_name * rel_uncertainty) for the 
    # processes and observables specified by procname and obsname. If u_name has a Gaussian prior with
    # width 1.0 and mean 0.0 this is effectively a log-normal uncertainty.
    #
    # note that rel_uncertainty_plus and rel_uncertainty_minus usually have the same sign(!), unless the change in acceptance
    # goes in the same direction for both directions of the underlying uncertainty parameters.
    def add_asymmetric_lognormal_uncertainty(self, u_name, rel_uncertainty_minus, rel_uncertainty_plus, procname, obsname='*'):
        found_match = False
        par_name = u_name
        if par_name not in self.distribution.get_parameters():
            self.distribution.set_distribution(par_name, 'gauss', mean = 0.0, width = 1.0, range = [-float("inf"), float("inf")])
        for o in self.observable_to_pred:
            if obsname != '*' and o!=obsname: continue
            for p in self.observable_to_pred[o]:
                if procname != '*' and procname != p: continue
                self.observable_to_pred[o][p]['coefficient-function'].add_factor('exp', parameter = par_name, lambda_plus = rel_uncertainty_plus, lambda_minus = rel_uncertainty_minus)
                found_match = True
        if not found_match: raise RuntimeError, 'did not find obname, procname = %s, %s' % (obsname, procname)
    

    # shortcut for add_asymmetric_lognormal_uncertainty for the symmetric case.
    def add_lognormal_uncertainty(self, u_name, rel_uncertainty, procname, obsname='*'):
        self.add_asymmetric_lognormal_uncertainty(u_name, rel_uncertainty, rel_uncertainty, procname, obsname)

        
    # get the set of parameters the model predictions depends on, given the signal_processes.
    # This does not cover additional_nll_terms; it is useful to pass the result to
    # Distribution.get_cfg.
    def get_parameters(self, signal_processes):
        result = set()
        for sp in signal_processes:
            assert sp in self.signal_processes
        for o in self.observable_to_pred:
            pred = self.observable_to_pred[o]
            for proc in pred:
                if proc in self.signal_processes and proc not in signal_processes: continue
                histo_pars = pred[proc]['histogram'].get_parameters()
                coeff_pars = pred[proc]['coefficient-function'].get_parameters()
                for par in histo_pars: result.add(par)
                for par in coeff_pars: result.add(par)
        if len(signal_processes) > 0: result.add('beta_signal')
        if self.additional_nll_term is not None:
            result.update(self.additional_nll_term.get_parameters())
        return result
    
    # returns two sets: the rate and shape changing parameters (can overlap!)
    # does not include beta_signal
    def get_rate_shape_parameters(self):
        rc, sc = set(), set()
        for o in self.observable_to_pred:
            pred = self.observable_to_pred[o]
            for proc in pred:
                sc.update(pred[proc]['histogram'].get_parameters())
                rc.update(pred[proc]['coefficient-function'].get_parameters())
        return rc, sc
    
    # options supported: use_llvm (default: False)
    #
    # signal_prior_cfg is the theta config dictionary
    def get_cfg2(self, signal_processes, signal_prior_cfg, options):
        result = {}
        if options.getboolean('model', 'use_llvm'): result['type'] = 'llvm_model'
        for sp in signal_processes:
            assert sp in self.signal_processes
        for o in self.observable_to_pred:
            result[o] = {}
            for proc in self.observable_to_pred[o]:
                if proc in self.signal_processes and proc not in signal_processes: continue
                result[o][proc] = {'histogram': self.observable_to_pred[o][proc]['histogram'].get_cfg(),
                     'coefficient-function': self.observable_to_pred[o][proc]['coefficient-function'].get_cfg()}
                if proc in signal_processes:
                    result[o][proc]['coefficient-function']['factors'].append('beta_signal')
        parameters = self.get_parameters(signal_processes)
        result['parameter-distribution'] = product_distribution(self.distribution.get_cfg(parameters), signal_prior_cfg)
        #rv observables:
        rvobservables = self.rvobs_distribution.get_parameters()
        if len(rvobservables) > 0:
            result['rvobs-distribution'] = self.rvobs_distribution.get_cfg(rvobservables)
        if self.bb_uncertainties: result['bb_uncertainties'] = True
        if self.additional_nll_term: result['additional_nll_term'] = self.additional_nll_term.get_cfg()
        return result
        
        
    def get_cfg(self, signal_processes = [], **options):
        result = {}
        if options.get('use_llvm', False): result['type'] = 'llvm_model'
        for sp in signal_processes:
            assert sp in self.signal_processes
        for o in self.observable_to_pred:
            result[o] = {}
            for proc in self.observable_to_pred[o]:
                if proc in self.signal_processes and proc not in signal_processes: continue
                result[o][proc] = {'histogram': self.observable_to_pred[o][proc]['histogram'].get_cfg(),
                     'coefficient-function': self.observable_to_pred[o][proc]['coefficient-function'].get_cfg()}
                if proc in signal_processes:
                    result[o][proc]['coefficient-function']['factors'].append('beta_signal')
        if self.bb_uncertainties: result['bb_uncertainties'] = True
        if self.additional_nll_term: result['additional_nll_term'] = self.additional_nll_term.get_cfg()
        return result


# a function multiplying 1d functions.
# TODO: Replace by FunctionBase ...
class Function:
    def __init__(self):
        self.value = 1.0
        self.factors = {} # map par_name -> theta cfg dictionary
    
    # supported types are: 'exp', 'id', 'constant'
    # pars:
    # * for typ=='exp': 'parameter', either 'lmbda' or 'lambda_plus' and 'lambda_minus'
    # * for typ=='id': 'parameter'
    # * for typ=='constant': value
    def add_factor(self, typ, **pars):
        assert typ in ('exp', 'id', 'constant')
        if typ=='constant':
            self.value *= pars['value']
        elif typ=='exp':
            if 'lmdba' in pars:
                self.factors[pars['parameter']] = {'type': 'exp_function', 'parameter': pars['parameter'], 'lambda_plus': pars['lmbda'], 'lambda_minus': float(pars['lmbda'])}
            else:
                self.factors[pars['parameter']] = {'type': 'exp_function', 'parameter': pars['parameter'], 'lambda_minus': float(pars['lambda_minus']), 'lambda_plus':float(pars['lambda_plus'])}
        elif typ=='id':
            p = pars['parameter']
            self.factors[p] = p
            
    def remove_parameter(self, par_name):
        if par_name in self.factors: del self.factors[par_name]

    # get a tuple (lambda_minus, lambda_plus) for the parameter par_name. Will throw an exception if par_name does not correspond to a 
    # exp function factor
    def get_exp_coeffs(self, par_name):
        assert self.factors[par_name]['type'] == 'exp_function'
        return self.factors[par_name]['lambda_minus'], self.factors[par_name]['lambda_plus']

    def set_exp_coeffs(self, par_name, lambda_minus, lambda_plus):
        assert self.factors[par_name]['type'] == 'exp_function'
        self.factors[par_name]['lambda_minus'], self.factors[par_name]['lambda_plus'] = lambda_minus, lambda_plus

    
    def rename_parameter(self, current_name, new_name):
        if current_name in self.factors:
            self.factors[new_name] = self.factors[current_name]
            del self.factors[current_name]
            if type(self.factors[new_name])==str:
                assert(self.factors[new_name] == current_name)
                self.factors[new_name] = new_name
            elif type(self.factors[new_name]) == dict:
                self.factors[new_name]['parameter'] = new_name
            else: assert(type(self.factors[new_name])==float)

    def get_value(self, par_values):
        result = self.value
        for p in self.factors:
            if self.factors[p] is p: result *= par_values[p]
            elif self.factors[p]['type'] == 'exp_function':
                if par_values[p] > 0: result *= math.exp(self.factors[p]['lambda_plus'] * par_values[p])
                else: result *= math.exp(self.factors[p]['lambda_minus'] * par_values[p])
            else: raise RuntimeError, 'unknown factor!'
        return result
    
    def get_cfg(self, optimize = True):
        result = {'type': 'multiply', 'factors': self.factors.values()}
        #print result
        # optimize by using several exponentials together:
        if optimize:
            result['factors'] = []
            parameters = [];
            lambdas_plus = [];
            lambdas_minus = [];
            for p in self.factors:
                if type(self.factors[p])!=dict or self.factors[p]['type'] != 'exp_function':
                    result['factors'].append(self.factors[p])
                    continue
                parameters.append(p)
                lambdas_plus.append(self.factors[p]['lambda_plus'])
                lambdas_minus.append(self.factors[p]['lambda_minus'])
            if len(parameters) > 0:
                 result['factors'].append({'type': 'exp_function', 'parameters': parameters, 'lambdas_plus': lambdas_plus, 'lambdas_minus': lambdas_minus})
        #print result
        if self.value != 1.0: result['factors'].append(self.value)
        return result
    
    def get_parameters(self):
        return self.factors.keys()

# modifies l(!)
def rebin_hlist(l, factor):
    assert len(l) % factor == 0
    new_len = len(l) / factor
    new_l = [sum(l[factor*i:factor*(i+1)]) for i in range(new_len)]
    while len(l) > new_len: del l[new_len]
    l[0:new_len] = new_l

def close_to(x1, x2):
    # if one is zero, both should be:
    if x1*x2 == 0: return x1==0.0 and x2==0.0
    # otherwise: relative difference small:
    return abs(x1-x2) / max(abs(x1), abs(x2)) < 1e-10

# returns a new h_triple. xmin, xmax might be slightly different from the requested values;
# the returned values are guaranteed to coincide with the bin borders
def set_hrange(h_triple, new_xmin, new_xmax):
    xmin, xmax, data = h_triple
    binsize = (xmax - xmin) / len(data)
    imin = int(new_xmin - xmin / binsize)
    imax = int(new_xmax - xmin / binsize)
    actual_new_xmin = xmin + imin * binsize
    actual_new_xmax = xmin + imax * binsize
    assert close_to(actual_new_xmin, new_xmin)
    assert close_to(actual_new_xmax, new_xmax)
    return (actual_new_xmin, actual_new_xmax, data[imin:imax])


# l+= coeff * l2 where l, l2 are lists of floats
def add_inplace(l, l2, coeff = 1.0):
    assert len(l) == len(l2)
    for i in range(len(l)): l[i] += coeff * l2[i]

# returns l + coeff * l2, without modifying l, l2
def add(l, l2, coeff = 1.0):
    assert len(l) == len(l2)
    return [l[i] + coeff * l2[i] for i in range(len(l))]



# a histogram class which stores the x range and 1D data.
# note that for backward-compatibility, this emulates a three-tuple
# (xmin, xmax, values) where xmin, xmax are floats and values is
# list-like (=int-indexed) container of floats
#
# It extends this by
# * adding uncertainties
# * providing some convenience functions like rebinning, scaling
#
# Note that those functions return a new Histogram instance.
class Histogram:
    def __init__(self, xmin, xmax, values, uncertainties = None, name = None):
        self.xmin = xmin
        self.xmax = xmax
        self.values = values
        self.name = name
        if uncertainties is not None: assert len(values) == len(uncertainties)
        self.uncertainties = uncertainties
        
    def get_uncertainties(self): return self.uncertainties
    
    def get_values(self): return self.values
    
    def get_value_sum(self): return sum(self.values)
    
    def get_value_sum_uncertainty(self):
        if self.uncertainties is None: return None
        return math.sqrt(sum([x**2 for x in self.uncertainties]))
        
    def get_xmin(self): return self.xmin
    def get_xmax(self): return self.xmax
    
    def get_cfg(self, include_error = True):
        result = {'type': 'direct_data_histo', 'range': [self.xmin, self.xmax], 'nbins': len(self.values), 'data': self.values}
        if self.uncertainties is not None and include_error: result['uncertainties'] = self.uncertainties
        return result
    
    def get_name(self): return self.name
    
    def copy(self):
        uncs = None if self.uncertainties is None else self.uncertainties[:]
        return Histogram(self.xmin, self.xmax, self.values[:], uncs, self.name)
    
    def strip_uncertainties(self): return Histogram(self.xmin, self.xmax, self.values, None, self.name)
    
    def scale(self, factor, new_name = None):
        uncs = None if self.uncertainties is None else array.array('d', [v * abs(factor) for v in self.uncertainties])
        return Histogram(self.xmin, self.xmax, array.array('d', [v * factor for v in self.values]), uncs, new_name)
        
    # calculate self + coeff * other_h and return the result as new Histogram (does not modify self).
    def add(self, coeff, other_h):
        assert (self.xmin, self.xmax, len(self.values)) == (other_h.xmin, other_h.xmax, len(other_h.values))
        result = Histogram(self.xmin, self.xmax, self.values[:])
        for i in range(len(self.values)):
            result.values[i] += coeff * other_h.values[i]
        if self.uncertainties is None:
            if other_h.uncertainties is None: pass
            else: result.uncertainties = array.array('d', [v * abs(coeff) for v in other_h.uncertainties])
        else:
            if other_h.uncertainties is None: result.uncertainties = self.uncertainties[:]
            else: result.uncertainties = array.array('d', [math.sqrt(u**2 + v**2 * coeff**2) for u,v in zip(self.uncertainties, other_h.uncertainties)])
        return result
        

    def rebin(self, rebin_factor, new_name = None):
        if len(self.values) % rebin_factor != 0: raise RuntimeError, "tried to rebin %d bins with factor %d" % (len(self.values), rebin_factor)
        newlen = len(self.values) / rebin_factor
        new_vals = array.array('d', [0.0] * newlen)
        uncs2 = [0.0] * newlen
        for i in range(len(self.values)):
            new_vals[i / rebin_factor] += self.values[i]
            if self.uncertainties is not None: uncs2[i / rebin_factor] += self.uncertainties[i]**2
        new_uncs = None
        if self.uncertainties is not None: new_uncs = array.array('d', map(math.sqrt, uncs2))
        return Histogram(self.xmin, self.xmax,  new_vals, new_uncs, new_name)
        
        
    def __len__(self): return 3
    
    def __getitem__(self, index):
        if index==0: return self.xmin
        if index==1: return self.xmax
        if index==2: return self.values
        raise KeyError
        
    def __iter__(self):
        class it:
            def __init__(self, h):
                self.index = 0
                self.h = h
            def __iter__(self): return self
            def next(self):
                if self.index==3: raise StopIteration
                res = self.h[self.index]
                self.index += 1
                return res
        return it(self)



# for morphing histos:
class HistogramFunction:
    def __init__(self, typ = 'cubiclinear'):
        assert typ in ('cubiclinear',)
        self.typ = typ
        self.parameters = set()
        self.nominal_histo = None
        self.factors = {} # parameter -> factor
        self.normalize_to_nominal = False
        self.syst_histos = {} # map par_name -> (plus histo, minus_histo)
        self.histrb = None # xmin, xmax nbins
        
    def get_nominal_histo(self): return self.nominal_histo
    def get_plus_histo(self, par): return self.syst_histos[par][0]
    def get_minus_histo(self, par): return self.syst_histos[par][1]
        

    def rename_parameter(self, current_name, new_name):
        if current_name not in self.parameters: return
        assert current_name != new_name
        self.parameters.remove(current_name)
        self.parameters.add(new_name)
        self.factors[new_name] = self.factors[current_name]
        del self.factors[current_name]
        self.syst_histos[new_name] = self.syst_histos[current_name]
        del self.syst_histos[current_name]
    
    # return a histogram triplet (xmin, xmax, data)
    def evaluate(self, par_values):
        if self.typ == 'cubiclinear':
            result = copy.deepcopy(self.nominal_histo)
            for p in self.parameters:
                delta = par_values[p] * self.factors[p]
                if abs(delta) > 1:
                    if delta > 0: h = add(self.syst_histos[p][0][2],self.nominal_histo[2], -1.0)
                    else: h = add(self.syst_histos[p][1][2],self.nominal_histo[2], -1.0)
                    add_inplace(result[2], h, abs(delta))
                else:
                    hsum = add(self.syst_histos[p][0][2], self.syst_histos[p][1][2])
                    add_inplace(hsum, self.nominal_histo[2], -2.0)
                    hdiff = add(self.syst_histos[p][0][2], self.syst_histos[p][1][2], -1.0)
                    diff = [0.5 * delta * x for x in hdiff]
                    add_inplace(diff, hsum, delta**2 - 0.5 * abs(delta)**3)
                    add_inplace(result[2], diff)
                 
            for i in range(len(result[2])):
                result[2][i] = max(0.0, result[2][i])
            if self.normalize_to_nominal:
                factor = sum(self.nominal_histo[2]) / sum(result[2])
                for i in range(len(result[2])): result[2][i] *= factor
            return result
        raise RuntimeError, "unknown typ '%s'" % self.typ
    
    # return a Histogram instance
    def evaluate2(self, par_values):
        if self.typ == 'cubiclinear':
            result = self.nominal_histo.copy()
            for p in self.parameters:
                delta = par_values[p] * self.factors[p]
                if abs(delta) > 1:
                    if delta > 0: h = self.nominal_histo.add(-1.0, self.syst_histos[p][0].strip_uncertainties())
                    else: h = self.nominal_histo.add(-1.0, self.syst_histos[p][1].strip_uncertainties())
                    result = result.add(abs(delta), h)
                else:
                    hsum = self.syst_histos[p][0].strip_uncertainties().add(1.0, self.syst_histos[p][1].strip_uncertainties())
                    hsum = hsum.add(-2.0, self.nominal_histo)
                    hdiff = self.syst_histos[p][0].strip_uncertainties().add(-1.0, self.syst_histos[p][1])
                    diff = hdiff.scale(0.5 * delta)
                    diff = diff.add(delta**2 - 0.5 * abs(delta)**3, hsum)
                    result = result.add(diff)
            for i in range(len(result[2])):
                result[2][i] = max(0.0, result[2][i])
            if self.normalize_to_nominal: result = result.scale(self.nominal_histo.get_value_sum() / result.get_value_sum())
            return result
        raise RuntimeError, "unknown typ '%s'" % self.typ
    
    
    # nominal_histo is a triple (xmin, xmax, data)
    def set_nominal_histo(self, nominal_histo, reset_binning = False):
        if reset_binning:
            self.histrb = None
            assert len(self.syst_histos) == 0
        histrb = nominal_histo[0], nominal_histo[1], len(nominal_histo[2])
        if self.histrb is None: self.histrb = histrb
        assert histrb == self.histrb, "histogram range / binning inconsistent!"
        self.nominal_histo = nominal_histo
    
    def set_syst_histos(self, par_name, plus_histo, minus_histo, factor = 1.0):
        self.parameters.add(par_name)
        self.factors[par_name] = factor
        histrb_plus = plus_histo[0], plus_histo[1], len(plus_histo[2])
        histrb_minus = minus_histo[0], minus_histo[1], len(minus_histo[2])
        assert histrb_plus == histrb_minus, "histogram range / binning inconsistent between plus / minus histo"
        if self.histrb is None: self.histrb = histrb_plus
        assert histrb_plus == self.histrb, "histogram range / binning inconsistent!"
        self.syst_histos[par_name] = (plus_histo, minus_histo)
    
    def scale(self, factor):
        self.nominal_histo = self.nominal_histo.scale(factor)
        for s in self.syst_histos:
            self.syst_histos[s] = self.syst_histos[s][0].scale(factor), self.syst_histos[s][1].scale(factor)
    """    
    def rebin(self, rebin_factor):
        rebin_hlist(self.nominal_histo[2], rebin_factor)
        for hp, hm in self.syst_histos.itervalues():
            rebin_hlist(hp, rebin_factor)
            rebin_hlist(hm, rebin_factor)
    """

    def get_cfg(self):
        if len(self.syst_histos) == 0:
            return self.nominal_histo.get_cfg()
        result = {'type': 'cubiclinear_histomorph', 'parameters': sorted(list(self.parameters)),
           'nominal-histogram': self.nominal_histo.get_cfg(), 'normalize_to_nominal': self.normalize_to_nominal}
        if set(self.factors.values()) != set([1.0]):
            result['parameter_factors'] = []
            for p in result['parameters']:
                result['parameter_factors'].append(self.factors[p])
        for p in self.parameters:
            result['%s-plus-histogram' % p] = self.syst_histos[p][0].get_cfg(False)
            result['%s-minus-histogram' % p] = self.syst_histos[p][1].get_cfg(False)
        return result

    def get_parameters(self): return self.parameters
    
    
# product of 1d distributions
class Distribution:
    def __init__(self):
        self.distributions = {} # map par_name -> dictionary with parameters 'mean', 'width', 'range', 'typ'
    
    # supported types: 'gauss', 'gamma'
    # note that width can be infinity or 0.0 to get flat and delta distributions, resp. In this case, the
    # type does not matter
    def set_distribution(self, par_name, typ, mean, width, range):
        assert typ in ('gauss', 'gamma')
        assert range[0] <= range[1]
        if type(mean) != str: # otherwise, the mean is a parameter ...
            mean = float(mean)
            assert range[0] <= mean and mean <= range[1] and width >= 0.0
        else: assert typ == 'gauss', "using a parameter as mean is only supported for gauss"
        self.distributions[par_name] = {'typ': typ, 'mean': mean, 'width': float(width), 'range': [float(range[0]), float(range[1])]}
        
    # Changes parameters of an existing distribution. pars can contain 'typ', 'mean', 'width', 'range'. Anything
    # not specified will be unchanged
    def set_distribution_parameters(self, par_name, **pars):
        assert par_name in self.distributions
        assert set(pars.keys()).issubset(set(['typ', 'mean', 'width', 'range']))
        self.distributions[par_name].update(pars)
        
    def get_distribution(self, par_name): return copy.deepcopy(self.distributions[par_name])
    
    def get_parameters(self): return self.distributions.keys()

    def rename_parameter(self, current_name, new_name):
        if current_name not in self.distributions: return
        self.distributions[new_name] = self.distributions[current_name]
        del self.distributions[current_name]

    def remove_parameter(self, par_name):
        del self.distributions[par_name]
    
    
    def get_means(self):
        result = {}
        for p in self.distributions:
            result[p] = self.distributions[p]['mean']
        return result
            
    
    # merged two distribution by preferring the distribution from dist1 over those from dist0.
    # override controls how merging of the distribution for a parameter in both dist1 and dist2 is done:
    #  * if override is True, dist1 takes preference over dist0
    #  * if override is False, dist0 and dist1 must define the same distribution (otherwise, an exception is thrown)
    @staticmethod
    def merge(dist0, dist1, override=True):
        result = Distribution()
        all_pars = set(dist0.get_parameters() + dist1.get_parameters())
        for p in all_pars:
            if p in dist1.distributions and p not in dist0.distributions: result.distributions[p] = dist1.get_distribution(p)
            elif p not in dist1.distributions and p in dist0.distributions: result.distributions[p] = dist0.get_distribution(p)
            else:
                if override: result.distributions[p] = dist1.get_distribution(p)
                else:
                    d0 = dist0.get_distribution(p)
                    assert d0 == dist1.get_distribution(p), "distributions for parameter '%s' not the same" % p
                    result.distributions[p] = d0
        return result
    
    def get_cfg(self, parameters):
        result = {'type': 'product_distribution', 'distributions': []}
        flat_dist = {'type': 'flat_distribution'}
        delta_dist = {'type': 'delta_distribution'}
        set_parameters = set(parameters)
        set_parameters.discard('beta_signal')
        set_self_parameters = set(self.distributions.keys())
        assert set_parameters.issubset(set_self_parameters), "Requested more parameters than distribution" + \
             " defined for: requested %s, got %s (too much: %s)" % (set_parameters, set_self_parameters, set_parameters.difference(set_self_parameters))
        for p in self.distributions:
            if p not in parameters: continue
            d = self.distributions[p]
            if d['width'] == 0.0:
                delta_dist[p] = d['mean']
            elif d['width'] == float("inf"):
                flat_dist[p] = {'range': d['range'], 'fix-sample-value': d['mean']}
            else:
                theta_type = {'gauss': 'gauss1d', 'gamma': 'gamma_distribution'}[d['typ']]
                result['distributions'].append({'type': theta_type, 'parameter': p, 'mean': d['mean'], 'width': d['width'], 'range': d['range']})
        if len(flat_dist) > 1:
            result['distributions'].append(flat_dist)
        if len(delta_dist) > 1:
            result['distributions'].append(delta_dist)
        return result


# return a Distribution object in which all parameters are fixed to their default values, using the Distribution
# template_dist.
def get_fixed_dist(template_dist):
    result = Distribution()
    for p in template_dist.get_parameters():
        val = template_dist.get_distribution(p)['mean']
        result.set_distribution(p, 'gauss', val, 0.0, [val, val])
    return result

# return a Distribution in which all parameters are fixed to the value given in par_values; par_values
# is a dictionary (parameter name) -> (parameter value)
def get_fixed_dist_at_values(par_values):
    result = Distribution()
    for p in par_values:
        val = par_values[p]
        result.set_distribution(p, 'gauss', val, 0.0, [val, val])
    return result


## \brief Build a multi-channel model based on template morphing from histograms in a root file
#
# This root file is expected to contain all the templates of the model adhering to a certain naming scheme:
# <observable>__<process>     for the "nominal" templates (=not affect by any uncertainty) and
# <observable>__<process>__<uncertainty>__(plus,minus)  for the "shifted" templates to be used for template morphing.
#
# ("observable" in theta is a synonym for what other call "channel": a statistically independent selection).
#
# <observable>, <process> and <uncertainty> are names you can choose at will as long as it does not contain '__'. You are encouraged
# to choose sensible names as these names are used in the output a lot.
#
# For example, if you want to make a combined statistical evaluation of a muon+jets and an electron+jets ttbar cross section measurement,
# you can name the observables "mu" and "ele"; the processes might be "ttbar", "w", "nonw", the uncertainties "jes", "q2". Provided
# all uncertainties affect all template shapes, you would supply 6 nominal and 24 "uncertainty" templates:
# The 6 nominal would be: mu__ttbar, mu__w, mu__nonw, ele__ttbar, ele__w, ele__nonw
# Some of the 24 "uncertainty" histograms would be: mu__ttbar__jes__plus, mu__ttbar__jes__minus, ..., ele__nonw__q2__minus
#
# theta-auto.py expects that all templates of one observable have the same range and binning. All templates should be normalized
# to the same luminosity (although normalization can be changed from the analysis python script later, this is generally not recommended, unless
# scaling everything to a different lumi).
#
# It is possible to omit some of the systematic templates completely. In this case, it is assumed
# that the presence of that uncertainty has no influence on this process in this observable.
#
# Data has the special process name "DATA" (all capitals!), so for each observable, there should be exactly one "<observable>_DATA"
# histogram, if you have data at all. If you do not have data, just omit this; the methods will be limited to calculating the expected
# result.
#
# To identify which process should be considered as signal, use as process name
# "signal". If you have multiple, uncorrelated signal scenarios (such as a new type of particle
# with different masses / widths), call it "signal<something>" where <something> is some name you can choose.
# Note that for multiple signals, use a number if possible. This number will be used in summary plots if investigating multiple signals.
# 
# For example, if you want to search for Higgs with masses 100, 110, 120, ..., you can call the processes "signal100",
# "signal110".
# Note that this naming convention for signal can be overriden by analysis.py via Model.set_signal_processes.
#
# The model built is based on the given templates where the systematic uncertainties are fully correlated across different
# observables and processes, i.e., the same parameter is used to interpolate between the nominal and shifted templates
# if the name of the uncertainty is the same. Two different systematic uncertainties (=with different names) are assumed to be uncorrelated.
# Each parameter has a Gaussian prior with width 1.0 and mean 0.0 and has the same name as the uncertainty. You can use
# the functions in Distribution (e.g., via model.distribution) to override this prior. This is useful if the "plus" and "minus" templates
# are not the +-1sigma deviations, but, say, the +-2sigma in which case you can use a prior with width 0.5.
#
# * histogram_filter is a function which -- given a histogram name as in the root file --
#   returns either True (=keep histogram) or False (do not keep histogram). the default is to keep all histograms.
# * root_hname_to_convention maps histogram names (as in the root file) to histogram names as expected by the
#    naming convention as described above. The default is to not modify the names.
# * transform_histo should transform the tuple (name, xmin, xmax, data) to a new tuple (name, xmin, xmax, data). name is not allowed to change.
#   name is the "internal" name  (=according to convention, after root_hname_to_convention).
#   This can be used e.g. to rebin histograms, set the range differently, etc.
# * include_mc_uncertainties is a boolean. If set to true, the bin-by-bin uncertainties will be included (and treated with the "Barlow-Beeston light" method).
#
# filenames can be either a string or a list of strings
#
# note: if using include_mc_uncertainties = True, transform_histo will only be called for the actual histogram values, not for the "uncertainty" histogram.
# Therefore, rebinning via transform_histo will not work if using include_mc_uncertainties.
def build_model_from_rootfile(filenames, histogram_filter = lambda s: True, root_hname_to_convention = lambda s: s, transform_histo = lambda h: h, include_mc_uncertainties = False):
    if type(filenames)==str: filenames = [filenames]
    result = Model()
    histos = {}
    observables, processes, uncertainties = set(), set(), set()

    for fname in filenames:
        rf = rootfile(fname)
        templates = rf.get_all_templates(include_mc_uncertainties)
        for hexternal in templates:
            if not histogram_filter(hexternal): continue
            hname_theta = root_hname_to_convention(hexternal)
            h_rootfile = templates[hexternal]
            h_mine = h_rootfile.copy()
            h_mine.name = hname_theta
            h_mine = transform_histo(h_mine)
            assert h_mine.get_name() == hname_theta, "transform_histo changed the name. This is not allowed; use root_hname_to_convention!"
            l = hname_theta.split('__')
            observable, process, uncertainty, direction = [None]*4
            if len(l)==2:
                observable, process = map(utils.transform_name_to_theta, l)
                observables.add(observable)
                processes.add(process)
            elif len(l)==4:
                observable, process, uncertainty, direction = l
                observable, process = map(utils.transform_name_to_theta, [observable, process])
            else:
                print "Warning: ignoring template %s (was: %s) which does not obey naming convention!" % (hname_theta, hexternal)
                continue
            if direction not in (None, 'plus', 'minus', 'up', 'down'):
                print "Warning: ignoring template %s (was: %s) which does not obey naming convention!" % (hname_theta, hexternal)
                continue
            if process == 'DATA':
                assert len(l)==2
                result.set_data_histogram(observable, h_mine.strip_uncertainties())
                continue
            if uncertainty is not None: uncertainties.add(uncertainty)
            if direction=='up': direction='plus'
            if direction=='down': direction='minus'
            if uncertainty is not None: h_new = '%s__%s__%s__%s' % (observable, process, uncertainty, direction)
            else: h_new = '%s__%s' % (observable, process)
            histos[h_new] = h_mine

    # build histogram functions from templates, and make some sanity checks:
    for o in observables:
        for p in processes:
            hname_nominal = '%s__%s' % (o, p)
            if hname_nominal not in histos: continue
            hf = HistogramFunction()
            h = histos[hname_nominal]
            if not include_mc_uncertainties: h =  h.strip_uncertainties()
            hf.set_nominal_histo(h)
            for u in uncertainties:
                n_syst = 0
                if ('%s__%s__%s__plus' % (o, p, u)) in histos: n_syst += 1
                if ('%s__%s__%s__minus' % (o, p, u)) in histos: n_syst += 1
                if n_syst == 0: continue
                if n_syst != 2: raise RuntimeError, "only one direction given for (observable, process, uncertainty) = (%s, %s, %s)" % (o, p, u)
                hf.set_syst_histos('%s' % u, histos['%s__%s__%s__plus' % (o, p, u)], histos['%s__%s__%s__minus' % (o, p, u)])
            result.set_histogram_function(o, p, hf)
    for u in uncertainties:
        result.distribution.set_distribution('%s' % u, 'gauss', mean = 0.0, width = 1.0, range = (-float("inf"), float("inf")))
    result.bb_uncertainties = include_mc_uncertainties
    return result

