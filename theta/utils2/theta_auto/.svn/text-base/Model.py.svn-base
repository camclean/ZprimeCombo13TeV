# -*- coding: utf-8 -*-
import re, fnmatch, math, copy
import os, os.path
import array
import utils



def evaluate_prediction(model, par_values, include_signal = True, observables = None):
    """
    Get the model histograms, applying template morphing and scaling according to the given parameter values.
    
    * ``model`` is the :class:`Model` you want to evlaute
    * ``par_values`` is a dictionary with the parameter names as key and the parameter values as dictionary values
    * ``include_signal`` - if ``True``, all signal processes are included as well. Note that this requires a parameter value "beta_signal" in ``par_values``
    * ``observables`` - is a list of observables (=channels) to include. The default value ``None`` includes all observables in the model.
    
    The return value is a nested dictionary with the observable name as first-level key, the process name as the second level key.
    You can pass the result to :meth:`write_histograms_to_rootfile` to write these to a root file.
    """
    result = {}
    if observables is None: observables = model.get_observables()
    for obs in observables:
        result[obs] = {}
        for p in model.get_processes(obs):
            if p in model.signal_processes and not include_signal: continue
            coeff = model.get_coeff(obs, p).get_value(par_values)
            if p in model.signal_processes: coeff *= par_values['beta_signal']
            result[obs][p] = model.get_histogram_function(obs, p).evaluate(par_values).scale(coeff)
    return result


class Model(utils.Copyable):
    """
    The statistical model specifies the predicted event yields *λ_ci* for all channels *c* and bins *i*, as a function of the model
    parameters. This class also contains information about the observed data. See :ref:`model_intro` for details. 
    The Model can be seen as multiple mappings: each pair (channel, process) is assigned one coefficient, the :ref:`Function`
    and one :class:`HistogramFunction` object. In addition, `model.distribution` is a :class:`Distribution` instance contains information
    about the prior of the nuisance parameters.
    
    Usually, a Model instance is not constructed directly but rather via the methods :func:`build_model_from_rootfile`
    or :func:`higgs_datacard.build_model`.
    """
    
    def __init__(self):
        # observables is a dictionary str name ->  (float xmin, float xmax, int nbins)
        self.observables = {}
        self.processes = set() # all processe
        self.signal_processes = set() # the subset of processes considered signal
        # in addition to the signal_processes, we need to define which processes to consider together
        # as signal. This is defined in signal_process_groups. It is a dictionary from the signal process group is to
        # the list of signal processes to consider.
        self.signal_process_groups = {}
        # observable_to_pred is a nested dictionary (str observable) --> (str process)
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


    def filter_processes(self, filter_function):
        """
        Filter the process list: keep exactly those processes for which filter_function(p) returns True. The method
        will also clean up the observable list by removing observables without any process.
        """
        all_to_delete = set() # process list
        observables = list(self.observable_to_pred.keys())
        for o in observables:
            to_delete = set()
            for p in self.observable_to_pred[o]:
                keep = filter_function(p)
                if keep: continue
                to_delete.add(p)
                all_to_delete.add(p)
            for p in to_delete: del self.observable_to_pred[o][p]
            # delete observable o, if empty:
            if len(self.observable_to_pred[o]) == 0:
                del self.observable_to_pred[o]
                del self.observables[o]
                if o in self.data_histos: del self.data_histos[o]
        self.processes = self.processes.difference(all_to_delete)
        self.signal_processes = self.signal_processes.difference(all_to_delete)
    

    def reset_binning(self, obs, xmin, xmax, nbins):
        assert obs in self.observables
        self.observables[obs] = (xmin, xmax, nbins)
    
    def rename_observable(self, current_name, new_name):
        """
        Rename the observables (channel) to the new name.
        """
        assert current_name in self.observables
        self.observables[new_name] = self.observables[current_name]
        del self.observables[current_name]
        self.observable_to_pred[new_name] = self.observable_to_pred[current_name]
        del self.observable_to_pred[current_name]
        if current_name in self.data_histos:
            self.data_histos[new_name] = self.data_histos[current_name]
            del self.data_histos[current_name]
    
    def combine(self, other_model, strict=True):
        """
        Combines the current with the Model ``other_model`` by adding all channels (and processes) from ``other_model``.
        
        Note that:
        
        * ``other_model`` must not include an observable of the current model
        * ``other_model`` should define the same signal_processes (this is enforced if ``strict=True``)
        * ``other_model.distribution`` can include the same nuisance parameters (=same name).
    
        For shared nuisance parameters, the prior for self.distribution and other.distribution must be identical.
        """
        my_obs = set(self.observables.keys())
        other_obs = set(other_model.observables.keys())
        if len(my_obs.intersection(other_obs))!=0:
            print "models to be combined share observables, but they must not!"
            print "shared observables: ", my_obs.intersection(other_obs)
            raise RuntimeError, "models share observables"
        if strict:
            assert self.signal_processes == other_model.signal_processes, "signal processes not equal: left-right=%s; right-left=%s;" \
                       % (str(self.signal_processes.difference(other_model.signal_processes)), str(other_model.signal_processes.difference(self.signal_processes)))
        else:
            self.signal_processes.update(other_model.signal_processes)
            self.signal_process_groups.update(other_model.signal_process_groups)
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
        """
        Rename a nuisance parameter in the model. On general, this affects coefficiencts, histogram functions, and the prior parameter distribution
        in ``self.distrbution``.
        """

        self.distribution.rename_parameter(current_name, new_name)
        for o in self.observable_to_pred:
            for p in self.observable_to_pred[o]:
                self.observable_to_pred[o][p]['histogram'].rename_parameter(current_name, new_name)
                self.observable_to_pred[o][p]['coefficient-function'].rename_parameter(current_name, new_name)
    
    def restrict_to_observables(self, observables):
        """
        Modify ``self`` to be a model of a subset of the current channels (=obervables).
        
        The parameter ``observables`` should be a list or set of the observable names you want to keep. It must be a subset
        of the current set of observables.
        """
        observables = set(observables)
        model_observables = set(self.observables.keys())
        assert observables.issubset(model_observables), "observables must be a subset of model.observables!"
        for obs in model_observables:
            if obs in observables: continue
            del self.observables[obs]
            del self.observable_to_pred[obs]
            if obs in self.data_histos: del self.data_histos[obs]
        # in theory, we could also update self.processes / self.signal_processes and
        # self.distribution (if the model now has fewer dependencies), but this is not necessary,
        # as the creation of the .cfg file will only ask for the parameters the current model depends on anyway.
       
    def set_data_histogram(self, obsname, histo, reset_binning = False):
        """
        Set the data histogram to the :class:`Histogram` histo for the observable given in ``obsname``. If
        ``reset_binning`` is ``True``, the binning for the observable ``obsname`` will be set to the binning of the
        given ``histo``.
        """
        #TODO: the whole reset_binning should be replaced by a more transparent method of handling
        # the binning (?)
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
        """
        Get the data histogram for the observable ``obsname`` previously set by :meth:Model.set_data_histogram.
        Returns ``None`` if no data histogram was set for this observable (=channel).
        """
        if obsname in self.data_histos: return self.data_histos[obsname]
        else: return None
        
        
    def has_data(self):
        """
        Retuns ``True`` if and only if a data histogram has been speficied for all observables (=channels).
        """
        for o in self.observables:
            if o not in self.data_histos: return False
        return True
    
    def get_data_rvobsvalues(self):
        return self.data_rvobsvalues

    def get_observables(self):
        """
        Get the set of observables (=channels).
        """
        return set(self.observables.keys())
        
    def get_range_nbins(self, obsname):
        """
        return triple ``(xmin, xmax, nbins)`` for the given observable
        """
        return self.observables[obsname]
    
    #important: always call set_histogram_function first, then get_coeff!
    def set_histogram_function(self, obsname, procname, histo):
        xmin, xmax, nbins = histo.get_xmin_xmax_nbins()
        if obsname not in self.observables:
            self.observables[obsname] = xmin, xmax, nbins
            self.observable_to_pred[obsname] = {}
        xmin2, xmax2, nbins2 = self.observables[obsname]
        assert (xmin, xmax, nbins) == (xmin2, xmax2, nbins2), "detected inconsistent binning setting histogram for (obs, proc) = (%s, %s): new binning: %s, old binning: %s" % (obsname, procname,
           (xmin, xmax, nbins), (xmin2, xmax2, nbins2))
        self.processes.add(procname)
        if procname not in self.observable_to_pred[obsname]: self.observable_to_pred[obsname][procname] = {}
        self.observable_to_pred[obsname][procname]['histogram'] = histo
        if 'coefficient-function' not in self.observable_to_pred[obsname][procname]: self.observable_to_pred[obsname][procname]['coefficient-function'] = Function()
        
    def get_coeff(self, obsname, procname):
        return self.observable_to_pred[obsname][procname]['coefficient-function']
        
    def get_processes(self, obsname):
        """
        Get the list of process names contributing to the given channel ``obsname``.
        """
        return self.observable_to_pred[obsname].keys()
        
        
    def get_signal_processes(self):
        """
        Get a set of all signal names in this model.
        """
        return set(self.signal_processes)
        
    def get_histogram_function(self, obsname, procname):
        if procname not in self.observable_to_pred[obsname]: return None
        return self.observable_to_pred[obsname][procname]['histogram']
    
    
    def fill_histogram_zerobins(self, epsilon = 0.001):
        """
        Fill all histogram bins with at least epsilon * (average bin content of histogram h).
        This can be used to prevent a prediction of exactly 0 in a bin which will lead to
        a zero likelihood / infinite negative log-likelihood.
    
        The special case ``epsilon = None`` will leave the prediction at the value of 0.0
        but will make sure that the uncertainty in each bin is at least the one corresponding 1 MC event.
        The average weight for the MC sample is estimated from the uncertainties given in the histogram.
        This makes only sense if enabling the Barlow-Beeston treatment of MC statistical uncertainties.
        """
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

    def scale_predictions(self, factor, procname = '*', obsname = '*'):
        """
        Scale the templates of the given observable ``obsname`` and process ``procname`` by ``factor``. This will
        only scale the prediction, not the data.
        """
        found_match = False
        for o in self.observable_to_pred:
            if obsname != '*' and o!=obsname: continue
            for p in self.observable_to_pred[o]:
                if procname != '*' and procname != p: continue
                found_match = True
                self.observable_to_pred[o][p]['histogram'].scale(factor)
        if not found_match: raise RuntimeError, 'did not find obname, procname = %s, %s' % (obsname, procname)

    def rebin(self, obsname, rebin_factor):
        """
        Rebin the observable ``obsname`` by the given rebinning factor.
        """
        xmin, xmax, nbins_old = self.observables[obsname]
        assert nbins_old % rebin_factor == 0
        self.observables[obsname] = (xmin, xmax, nbins_old / rebin_factor)
        pred = self.observable_to_pred[obsname]
        for p in pred:
            pred[p]['histogram'].rebin(rebin_factor)
        if obsname in self.data_histos:
            self.data_histos[obsname] = self.data_histos[obsname].rebin(rebin_factor, self.data_histos[obsname].get_name())
    
    def set_signal_processes(self, procs):
        """
        Define which processes should be considered as signal processes
    
        Any process not defined explicitly as signal is considered to be background.
        
        This method assumes that you want to treat each signal process in turn.
        For more control in situations where the signal has multiple histograms to be used simultaneously,
        use :meth:`Model.set_signal_process_groups`.
    
        ``procs`` is a list / set of glob patterns (or a single pattern).
        """
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
        
    
    def set_signal_process_groups(self, groups):
        """
        Define the signal process groups
   
        ``groups`` is a dictionary (id) --> (list of processes), see :ref:`what_is_signal`.
        """
        for spid in groups:
            for p in groups[spid]: assert type(p)==str and p in self.processes, "unknown process '%s'" % p
        self.signal_processes = set()
        for spid in groups:
            self.signal_processes.update(groups[spid])
        self.signal_process_groups = copy.deepcopy(groups)

    
    def add_asymmetric_lognormal_uncertainty(self, u_name, rel_uncertainty_minus, rel_uncertainty_plus, procname, obsname='*'):
        """
        Add a rate-only uncertainty for the given combination of (process, observable)
    
        Adds a new parameter with name ``u_name`` to the ``distribution`` (unless it already exists) with
        a Gaussian prior around 0.0 with width 1.0 and adds a factor exp(u_name * rel_uncertainty) for the 
        (process, observable) combination specified by ``procname`` and ``obsname``. In effect, this is a log-normal uncertainty
        for the coefficient of the template.
        
        Note that ``rel_uncertainty_plus`` and ``rel_uncertainty_minus`` usually have the same sign, unless the change in acceptance
        goes in the same direction for both directions of the underlying uncertainty parameters.
        """
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
    

    def add_lognormal_uncertainty(self, u_name, rel_uncertainty, procname, obsname='*'):
        """
        shortcut for add_asymmetric_lognormal_uncertainty for the symmetric case.
        """
        self.add_asymmetric_lognormal_uncertainty(u_name, rel_uncertainty, rel_uncertainty, procname, obsname)

        
    
    def get_parameters(self, signal_processes):
        """
        Get the list of parameters the model predictions depends on. In general, this depends on
        which processes are considered as signal, therefore this has to be specified in the
        ``signal_processes`` parameter.
        
        Parameters:
        
        * ``signal_processes`` - a list/set of process names to consider signal
        """
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
        # also include mean values of rvobs:
        rvobs_means = self.rvobs_distribution.get_means().values()
        for m in rvobs_means:
            if type(m)==float: continue
            assert type(m)==str
            result.add(m)
        return sorted(list(result))
    
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
    def get_cfg(self, signal_processes, signal_prior_cfg, options):
        result = {}
        if options.getboolean('model', 'use_llvm'): result['type'] = 'llvm_model'
        for sp in signal_processes:
            assert sp in self.signal_processes
        for o in self.observable_to_pred:
            result[o] = {}
            for proc in self.observable_to_pred[o]:
                if proc in self.signal_processes and proc not in signal_processes: continue
                cf = self.observable_to_pred[o][proc]['coefficient-function']
                if proc in signal_processes:
                    cf = cf.copy()
                    cf.add_factor('id', parameter = 'beta_signal')
                result[o][proc] = {'histogram': self.observable_to_pred[o][proc]['histogram'].get_cfg(), 'coefficient-function': cf.get_cfg()}
        parameters = self.get_parameters(signal_processes)
        if 'beta_signal' in parameters:            
            # optimize away product distribution in case the signal prior is flat / fixed:
            if signal_prior_cfg['type'] in ('flat_distribution', 'fixed_distribution'):
                dist = self.distribution.copy()
                if signal_prior_cfg['type'] == 'flat_distribution': range, mean, width = signal_prior_cfg['beta_signal']['range'], signal_prior_cfg['beta_signal']['fix-sample-value'], float('inf')
                else:
                    mean, width = signal_prior_cfg['beta_signal'], 0.0
                    range = mean, mean
                dist.set_distribution('beta_signal', 'gauss', mean = mean, width = width, range = range)
                result['parameter-distribution'] = dist.get_cfg(parameters)
            else:
                bkg_parameters = set(parameters)
                bkg_parameters.discard('beta_signal')
                result['parameter-distribution'] = {'type': 'product_distribution', 'distributions': (self.distribution.get_cfg(bkg_parameters), signal_prior_cfg)}
        else:
            result['parameter-distribution'] = self.distribution.get_cfg(parameters)
        #rv observables:
        rvobservables = self.rvobs_distribution.get_parameters()
        if len(rvobservables) > 0:
            result['rvobs-distribution'] = self.rvobs_distribution.get_cfg(rvobservables)
        if self.bb_uncertainties: result['bb_uncertainties'] = True
        if self.additional_nll_term: result['additional_nll_term'] = self.additional_nll_term.get_cfg()
        return result


# a function multiplying 1d functions.
# TODO: Replace by FunctionBase ...
class Function(utils.Copyable):
    """
    Instances of this class are used as coefficients for templates in a model. It is limited to simple expressions
    of a constant factor ``c``, multiplied by exponential factors (which are in effect the log-normal rate uncertainties, if the
    corresponding nuisance parameter is Gaussian), and directly using parameters as factors:
    
    ..math::
    
       c\\times \\prod_p exp(\\lambda * p) \\times \\prod_p p
       
    where ``p`` are parameters.
    """
        
    def __init__(self):
        self.value = 1.0
        self.factors = {} # map par_name -> theta cfg dictionary
    

    def add_factor(self, typ, **pars):
        """
        supported types are: 'exp', 'id', 'constant'
        
       ``pars`` depends on the value of ``typ``:
       
        * for typ=='exp': 'parameter', either 'lmbda' or 'lambda_plus' and 'lambda_minus', for the exponential term
        * for typ=='id': 'parameter', the parameter name to use for the direct term
        * for typ=='constant': 'value', the floating point value to multiple the constant term to
        """
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
        if self.value != 1.0: result['factors'].append(self.value)
        #print result
        # optimize by using several exponentials together:
        if optimize:
            result['factors'] = []
            if self.value != 1.0: result['factors'].append(self.value)
            parameters = []
            lambdas_plus = []
            lambdas_minus = []
            for p in self.factors:
                if type(self.factors[p])!=dict or self.factors[p]['type'] != 'exp_function':
                    result['factors'].append(self.factors[p])
                    continue
                parameters.append(p)
                lambdas_plus.append(self.factors[p]['lambda_plus'])
                lambdas_minus.append(self.factors[p]['lambda_minus'])
            if len(parameters) > 0:
                 result['factors'].append({'type': 'exp_function', 'parameters': parameters, 'lambdas_plus': lambdas_plus, 'lambdas_minus': lambdas_minus})
        try:
            # this can raise in case result['factors'][0] is not a dictionary:
            if len(result['factors']) == 1 and result['factors'][0]['type']=='exp_function' and self.value == 1.0:
                result = result['factors'][0]
        except Exception: pass
        return result
    
    def get_parameters(self):
        return self.factors.keys()

class Histogram(object, utils.Copyable):
    """
    This class stores the x range and 1D data, and optionally bin-by-bin uncertainties. Its main
    use is in :class:`HistogramFunction` and as data histograms in :class:`Model`.
    
    Histograms are immutable: methods such as `scale` return the new Histogram instead of modifying
    the present instance.
    """
    
    def __init__(self, xmin, xmax, values, uncertainties = None, name = None, x_low = None):
        """
        Constructor.
        
        If not None, x_low should be an array of lower bin borders in case of non-equidistant binning. It overrides
        the value of xmin.
        """
        self.xmin = xmin
        self.xmax = xmax
        self.values = values
        self.name = name
        # in case of non-equidistant binning:
        self.x_low = None
        if uncertainties is not None:
            assert len(values) == len(uncertainties)
            self.uncertainties = uncertainties
        else:
            self.uncertainties = None
        if x_low is not None:
            assert len(x_low) == len(self.values)
            assert max(x_low) < self.xmax
            self.xmin = min(x_low)
            self.x_low = x_low
            
            
    def is_compatible(self, other):
        """
        Returns whether or not this histogram has the same range and binning as other.
        """
        res = self.xmin, self.xmax, len(self.values) == other.xmin, other.xmax, len(other.values)
        if not res: return False
        if self.x_low is not None:
            if other.x_low is None: return False
            for xl1, xl2 in zip(self.x_low, other.x_low):
                if xl1 != xl2: return False
        return True
        
    def get_uncertainties(self):
        """
        Return the array of uncertainties, or ``None``, if no uncertainties were defined.
        """
        return self.uncertainties
    
    def get_values(self):
        """
        Return the array of values.
        """
        return self.values
    
    def get_value_sum(self):
        """
        Calculate the sum of values.
        """
        return sum(self.values)
    
    def get_value_sum_uncertainty(self):
        """
        Get the uncertainty on the sum of values. Returns ``None`` if no uncertainties were defined.
        """
        if self.uncertainties is None: return None
        return math.sqrt(sum([x**2 for x in self.uncertainties]))
        
    def get_xmin(self):
        """
        Get the lower end of the histogram range
        """
        return self.xmin
    
    def get_xmax(self):
        """
        Get the upper end of the histogram range
        """
        return self.xmax
    
    def get_nbins(self):
        """
        Get the number of bins. Same as ``len(self.get_values())``
        """
        return len(self.values)
    
    def get_x_low(self, ibin):
        """
        Get the low end of the bin number ``ibin``. Bins are numberes from ``0`` to ``get_nbins-1``.
        """
        if self.x_low is not None: return self.x_low[ibin]
        else: return self.xmin + (self.xmax - self.xmin) / len(self.values) * ibin        
    
    def get_cfg(self, include_error = True):
        """
        Get a dictionary representing a HistogramFunction configuration for theta.
        """
        result = {'type': 'direct_data_histo', 'range': [self.xmin, self.xmax], 'nbins': len(self.values), 'data': self.values}
        if self.uncertainties is not None and include_error: result['uncertainties'] = self.uncertainties
        return result
    
    def get_name(self):
        """
        Get the name of the histogram
        """
        return self.name
    
    def strip_uncertainties(self):
        """
        Get a copy of the Histogram without uncertainties.
        """
        h = Histogram(self.xmin, self.xmax, self.values, None, self.name, x_low = self.x_low)
        return h
    
    def scale(self, factor, new_name = None):
        """
        Scale the histogram by the given factor ``factor``; optionally setting a new hitsogram name.
        """
        uncs = None if self.uncertainties is None else array.array('d', [v * abs(factor) for v in self.uncertainties])
        h = Histogram(self.xmin, self.xmax, array.array('d', [v * factor for v in self.values]), uncs, new_name, x_low = self.x_low)
        return h
        
    # calculate self + coeff * other_h and return the result as new Histogram (does not modify self).
    def add(self, coeff, other_h):
        """
        Calculate ``self + coeff * other_h`` and returns the resulting histogram. Uncertainties are propagated.
        """
        assert self.is_compatible(other_h)
        result = Histogram(self.xmin, self.xmax, self.values[:], x_low = self.x_low)
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
        """
        Return a rebinned histogram. Assumes that the uncertainties in different bins are uncorrelated.
        """
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


class ExponentialHistogramFunction:
    def __init__(self, lambda0, c, parameter, normalize_to, binborders):
        self.lambda0 = lambda0
        self.c = c
        self.parameter = parameter
        self.normalize_to = normalize_to
        self.binborders = binborders

    def rebin(self):
        raise RuntimeError, "rebinning not supported for ExponentialHistogramFunction"
        
    def remove_parameter(self, par):
        raise RuntimeError, "removing a parameter not supported for ExponentialHistogramFunction"
    
    def evaluate(self, parameters):
        p = parameters[self.parameter]
        lmbda = self.lambda0 + p * self.c
        nbins = len(self.binborders) - 1
        values = []
        s = 0.0
        for i in range(nbins):
            bincontent = abs(math.exp(lmbda * self.binborders[i+1]) - math.exp(lmbda * self.binborders[i]))
            values.append(bincontent)
            s += bincontent
        if s==0.0: s = 1.0
        return Histogram(self.binborders[0], self.binborders[-1], values, x_low = self.binborders[:-1]).scale(self.normalize_to / s)
            
    
    def rename_parameter(self, par_old, par_new):
        if self.parameter==par_old: self.parameter = par_new
        
    def get_cfg(self):
        result = {'type': 'exponential_hf', 'parameter' : self.parameter, 'lambda0' : self.lambda0, 'c' : self.c, 'normalize_to': self.normalize_to, 'binborders' : self.binborders }
        return result
        
    def get_parameters(self):
        return set([self.parameter])
        
    def get_xmin_xmax_nbins(self):
        return self.binborders[0], self.binborders[-1], len(self.binborders)-1
        
        


# for morphing histos:
class HistogramFunction:
    """
    A parameter-dependent Histogram. This is used in the :class:`Model` as building block.
    
    The parameter dependence currently modeled is a template morphing which has one "nominal" histogram
    and 2*n_syst "alternate" histograms where n_syst is the number of systematic uncertainties / the number of nuisance parameters
    included in the description of this histogram.
    
    It corresponds to the ``cubiclinear_histomorph`` plugin in theta; see details about the kind of morphing there.
    """
    def __init__(self, typ = 'cubiclinear'):
        assert typ in ('cubiclinear',)
        self.typ = typ
        self.parameters = set()
        self.nominal_histo = None
        self.factors = {} # parameter -> factor
        self.normalize_to_nominal = False
        self.syst_histos = {} # map par_name -> (plus histo, minus_histo)
        self.histrb = None # xmin, xmax nbins
    
    def get_xmin_xmax_nbins(self): return self.histrb
    def get_nominal_histo(self): return self.nominal_histo
    def get_plus_histo(self, par): return self.syst_histos[par][0]
    def get_minus_histo(self, par): return self.syst_histos[par][1]
    def get_factor(self, par): return self.factors[par]
    
    def rebin(self, factor):
        assert self.histrb[2] % factor == 0
        self.nominal_histo = self.nominal_histo.rebin(factor)
        for p in self.syst_histos:
            self.syst_histos[p] = self.syst_histos[p][0].rebin(factor), self.syst_histos[p][1].rebin(factor)
        self.histrb = self.histrb[0], self.histrb[1], self.histrb[2] / factor
    
    def remove_parameter(self, par_name):
        assert par_name in self.parameters
        self.parameters.discard(par_name)
        del self.factors[par_name]
        del self.syst_histos[par_name]
        
    def rename_parameter(self, current_name, new_name):
        """
        Rename a parameter. ``current_name`` does not need to exist in the list of parameters this HistogramFunction depends on.
        """
        if current_name not in self.parameters: return
        assert current_name != new_name
        self.parameters.remove(current_name)
        self.parameters.add(new_name)
        self.factors[new_name] = self.factors[current_name]
        del self.factors[current_name]
        self.syst_histos[new_name] = self.syst_histos[current_name]
        del self.syst_histos[current_name]
        
    def evaluate(self, par_values):
        """
        Return the Histogram evaluated at the given parameter values.
        
        ``par_values`` is a dictionary where the key is the parameter name (a string) and the value is the value of this parameter (a float).
        """
        if self.typ == 'cubiclinear':
            result = self.nominal_histo.copy()
            for p in self.parameters:
                delta = par_values[p] * self.factors[p]
                if abs(delta) > 1:
                    if delta > 0: h = self.nominal_histo.add(-1.0, self.syst_histos[p][0].strip_uncertainties())
                    else: h = self.nominal_histo.add(-1.0, self.syst_histos[p][1].strip_uncertainties())
                    result = result.add(-abs(delta), h)
                else:
                    hsum = self.syst_histos[p][0].strip_uncertainties().add(1.0, self.syst_histos[p][1].strip_uncertainties())
                    hsum = hsum.add(-2.0, self.nominal_histo)
                    hdiff = self.syst_histos[p][0].strip_uncertainties().add(-1.0, self.syst_histos[p][1])
                    diff = hdiff.scale(0.5 * delta)
                    diff = diff.add(delta**2 - 0.5 * abs(delta)**3, hsum)
                    result = result.add(1.0, diff.strip_uncertainties())
            for i in range(len(result[2])):
                result[2][i] = max(0.0, result[2][i])
            if self.normalize_to_nominal: result = result.scale(self.nominal_histo.get_value_sum() / result.get_value_sum())
            return result
        raise RuntimeError, "unknown typ '%s'" % self.typ
    
    def set_nominal_histo(self, nominal_histo, reset_binning = False):
        """
        Set the nominal Histogram
        
        ``nominal_histo`` is the :class:`Histogram` instance to use as the new nominal histogram.
        If ``reset_binning`` is ``True``, the internal binning is reset to the binning of ``nominal_histo``.
        Otherwise, ``nominal_histo`` must have a binning which is consistent with the binning sued so far.
        """
        if reset_binning:
            self.histrb = None
            assert len(self.syst_histos) == 0
        histrb = nominal_histo[0], nominal_histo[1], len(nominal_histo[2])
        if self.histrb is None: self.histrb = histrb
        assert histrb == self.histrb, "histogram range / binning inconsistent!"
        self.nominal_histo = nominal_histo
    
    def set_syst_histos(self, par_name, plus_histo, minus_histo, factor = 1.0):
        """
        Set the shifted Histograms for 
        """
        self.parameters.add(par_name)
        self.factors[par_name] = factor
        histrb_plus = plus_histo[0], plus_histo[1], len(plus_histo[2])
        histrb_minus = minus_histo[0], minus_histo[1], len(minus_histo[2])
        assert histrb_plus == histrb_minus, "histogram range / binning inconsistent between plus / minus histo"
        if self.histrb is None: self.histrb = histrb_plus
        assert histrb_plus == self.histrb, "histogram range / binning inconsistent for syst=%s" % par_name
        self.syst_histos[par_name] = (plus_histo, minus_histo)
    
    def scale(self, factor):
        self.nominal_histo = self.nominal_histo.scale(factor)
        for s in self.syst_histos:
            self.syst_histos[s] = self.syst_histos[s][0].scale(factor), self.syst_histos[s][1].scale(factor)

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
    

class GaussDistribution:
    """
    A multivariate Gauss distribution with known mean and covariance.
    """
    def __init__(self, parameters, mu, covariance, ranges = None):
        n = len(parameters)
        assert n== len(mu) and n==len(covariance)
        self.parameters = list(parameters)
        self.mean = mu
        self.covariance = covariance
        if ranges is None:
            self.ranges = [("-inf", "inf") for i in range(len(mu))]
        else: self.ranges = ranges
        
    def get_parameters(self):
        return self.parameters[:]
        
    def get_cfg(self, rpars):
        p_to_i = {}
        for par in rpars: p_to_i[par] = self.parameters.index(par)
        mu = [self.mean[p_to_i[par]] for par in rpars]
        cov = [[self.covariance[p_to_i[row_par]][p_to_i[col_par]] for col_par in rpars] for row_par in rpars]
        ranges = [self.ranges[p_to_i[p]] for p in rpars]
        return {'type': 'gauss', 'parameters': rpars, 'mean': mu, 'covariance': cov, 'ranges' : ranges}
    
    def __str__(self):
        return 'GaussDistribution(' + str(self.parameters) + ", " + str(self.mean) + ", " + str(self.covariance) + ")"
        
    def rename_parameter(self, old_name, new_name):
        if old_name not in self.parameters: return
        self.parameters[self.parameters.index(old_name)] = new_name

# product of 1d distributions
class Distribution(utils.Copyable):
    def __init__(self):
        # map par_name -> dictionary with parameters 'mean', 'width', 'range', 'typ'.
        self.distributions = {} 
    
    # supported type: so far only 'gauss'
    # note that width can be infinity or 0.0 to get flat and delta distributions, resp. In this case, the
    # type does not matter
    def set_distribution(self, par_name, typ, mean, width, range):
        assert typ == 'gauss'
        assert range[0] <= range[1]
        if mean is not None and type(mean) != str: # otherwise, the mean is a parameter ...
            mean = float(mean)
            assert range[0] <= mean and mean <= range[1] and width >= 0.0
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
        result = {'type': 'igauss', 'parameters': [], 'mu':[], 'sigma':[], 'ranges': []}
        for p in sorted(parameters):
            if p not in self.distributions: continue
            result['parameters'].append(p)
            result['mu'].append(self.distributions[p]['mean'])
            result['sigma'].append(self.distributions[p]['width'])
            result['ranges'].append(self.distributions[p]['range'])
        return result



def get_fixed_dist(template_dist):
    """
    return a Distribution object in which all parameters are fixed to their default values, using the Distribution template_dist.
    """
    result = Distribution()
    for p in template_dist.get_parameters():
        val = template_dist.get_distribution(p)['mean']
        result.set_distribution(p, 'gauss', val, 0.0, [val, val])
    return result


def get_fixed_dist_at_values(par_values):
    """
    return a Distribution object in which all parameters are fixed to the value given in par_values; par_values
    is a dictionary (parameter name) -> (parameter value)
    """
    result = Distribution()
    for p in par_values:
        val = par_values[p]
        result.set_distribution(p, 'gauss', val, 0.0, [val, val])
    return result


