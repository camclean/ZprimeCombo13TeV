# -*- coding: utf-8 -*-
import config, sqlite3, os.path, glob, re
import numpy.linalg, array, math

import theta_interface, plotutil

import scipy.stats

import Model

inf = float("inf")
cl_1sigma = 0.68268949213708585
cl_2sigma = 0.95449973610364158

# returns a dictionary float x --> signal process group id
#
# As default, the x value is the first number encountered in the signal process group id string. 
#
# options can contain a dictionary 'spid_to_xvalue' in which case that number is used.
def get_x_to_sp(spgids, **options):
    x_to_sp = {}
    next_x = 0
    for sp in spgids:
        if 'spid_to_xvalue' in options and sp in options['spid_to_xvalue']: x = options['spid_to_xvalue'][sp]
        else: x = extract_number(sp)
        if x is None:
            warning("cannot find x-value for signal process id '%s', using %d (try passing the option spid_to_xvalue = {'%s': XXX})" % (sp, next_x, sp))
            x = next_x
            next_x += 1
        x_to_sp[x] = sp
    assert len(spgids) == len(x_to_sp)
    return x_to_sp


def make_tuple(*args): return tuple(args)

# options: booleans load_root_ plugins  (default: True)
# use_llvm (default: False)
def get_common_toplevel_settings(**options):
    cfg_options = {'plugin_files': ['$THETA_DIR/lib/core-plugins.so']}
    if options.get('use_llvm', False):
        cfg_options['plugin_files'].append('$THETA_DIR/lib/llvm-plugins.so')
    if options.get('load_root_plugins', True):
        cfg_options['plugin_files'].append('$THETA_DIR/lib/root.so')
    toplevel_settings = {'options': cfg_options}
    return toplevel_settings

def reldiff(d1, d2):
    return abs(d1 - d2) / max(abs(d1), abs(d2))

# returns a default minimizer specification which should be pretty robust
def minimizer(need_error = True, always_mcmc = False, minimizer_insane = False):
    n_iterations = 1000
    if minimizer_insane:
        n_iterations = 100000
        always_mcmc = True
    minimizers = []
    #try, in this order: migrad, mcmc+migrad, simplex, mcmc+simplex.
    if not always_mcmc: minimizers.append({'type': 'root_minuit'})
    minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min0', 'iterations': n_iterations, 'after_minimizer': {'type': 'root_minuit'}})
    if not always_mcmc: minimizers.append({'type': 'root_minuit', 'method': 'simplex'})
    minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min1', 'iterations': n_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
    #minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min2', 'iterations': n_iterations * 10, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
    result = {'type': 'minimizer_chain', 'minimizers': minimizers}
    if need_error: result['last_minimizer'] = {'type': 'root_minuit'}
    return result


def mul_list(l, c):
    for i in range(len(l)):
        l[i] *= c

# in theta, names must begin with a letter and consist only of A-Za-z0-9_- and not contain '__'
def transform_name_to_theta(name):
    result = ''
    for c in name:
        if c >= 'a' and c <= 'z' or c >= 'A' and c <= 'Z' or c >='0' and c <='9' or c in ('-', '_'): result += c
        else: result += '_'
    result = result.replace('__', '_')
    if result[0] >= '0' and result[0] <= '9' or result[0]=='-': result = 'tn_' + result
    return result

# returns a Distribution object, given the model, signal process and nuisance_prior specification ('shape:X;rate:Y'...)
def nuisance_prior_distribution(model, spec):
    if spec.__class__ == Model.Distribution: result = Model.Distribution.merge(spec, spec)  # note: merging copies ...
    else:
        result = Model.Distribution()
        min_p0, max_p1 = len(spec), 0
        shape_spec, rate_spec = None, None
        if 'shape:' in spec:
            p0 = spec.find('shape:')
            p1 = spec.find(';', p0)
            if p1 == -1: p1 = len(spec)
            min_p0 = min(min_p0, p0)
            max_p1 = max(max_p1, p1)
            shape_spec = spec[p0 + len('shape:'):p1]
        if 'rate:' in spec:
            p0 = spec.find('rate:')
            p1 = spec.find(';', p0)
            if p1 == -1: p1 = len(spec)
            min_p0 = min(min_p0, p0)
            max_p1 = max(max_p1, p1)
            rate_spec = spec[p0 + len('rate:'):p1]
        if (min_p0, max_p1) != (0, len(spec)): raise RuntimeError, 'nuisance_prior_distribution: could not parse specification %s (%d, %d)' % (spec, min_p0, max_p1)
        rc_pars, sc_pars = model.get_rate_shape_parameters()
        if shape_spec is not None:
            if shape_spec not in ('fix', 'free'): raise RuntimeError, 'unknown shape specification "%s"' % shape_spec
            for par in sc_pars:
                dist_pars = model.distribution.get_distribution(par)
                if shape_spec == 'fix':
                    #print 'fixing %s' % par
                    dist_pars['width'] = 0.0
                    result.set_distribution(par, **dist_pars)
                elif shape_spec == 'free':
                    #print 'freeing %s' % par
                    dist_pars['width'] = float("inf")
                    result.set_distribution(par, **dist_pars)
        if rate_spec is not None:
            if rate_spec not in ('fix', 'free'): raise RuntimeError, 'unknown rate specification "%s"' % rate_spec
            for par in rc_pars:
                dist_pars = model.distribution.get_distribution(par)
                if rate_spec == 'fix':
                    #print 'fixing %s' % par
                    dist_pars['width'] = 0.0
                    result.set_distribution(par, **dist_pars)
                elif rate_spec == 'free':
                    #print 'freeing %s' % par
                    dist_pars['width'] = float("inf")
                    result.set_distribution(par, **dist_pars)
    result = Model.Distribution.merge(model.distribution, result)
    return result

# get mean and rms estimate from list l:
def get_mean_width(l):
   n = len(l) * 1.0
   assert n > 0
   mean = sum(l) / n
   if n == 1: width = float('inf')
   else: width = math.sqrt(sum([(x - mean)**2 for x in l]) / (n-1))
   return mean, width

# get truncated mean / width similar to MarkovChainMC method from combine:
def get_trunc_mean_width(l):
   n = len(l)
   l_sorted = sorted(l)
   assert n > 0
   median = l_sorted[n/2]
   width = l_sorted[3*n/4] - l_sorted[n/4]
   l2 = [x for x in l if x >= median-width and x <= median+width]
   return get_mean_width(l2)

   
def p_to_Z(p_value):
   return -scipy.stats.norm.ppf(p_value)
   
def Z_to_p(z_value):
    return scipy.stats.norm.sf(z_value)
   
   
# returns the theta config dictionary, given a signal_prior specification ('flat' / 'fix:X')
def signal_prior_dict(spec):
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
            signal_prior_dict = {'type': 'flat_distribution', 'beta_signal': {'range': [xmin, xmax], 'fix-sample-value': value}}
        elif spec.startswith('fix:'):
            v = float(spec[4:])
            signal_prior_dict = theta_interface.delta_distribution(beta_signal = v)
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if type(spec) != dict: raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a dictionary!"
        signal_prior_dict = spec
    return signal_prior_dict

# returns a tuple (data_source dict, beta_signal distribution dict), given the input_spec ('toys:X' / 'toys-asimov:X' / 'data')
# the name is always 'source'.
#
# possible options: 'toydata_seed'
def data_source_dict(model, input_spec, **options):
    if input_spec == 'data':
        source_dict = {'type': 'histo_source', 'name': 'source'}
        for o in model.observables:
            source_dict[o] = model.data_histos[o].get_cfg(False)
        if len(model.data_rvobsvalues) > 0:
            source_dict['rvobs-values'] = dict(model.data_rvobsvalues)
        return (source_dict, theta_interface.delta_distribution(beta_signal = 0.0))
    elif input_spec.startswith('toys:') or input_spec.startswith('toys-asimov:'):
        beta_signal_value = float(input_spec[input_spec.find(':') + 1:])
        seed = 1
        if 'toydata_seed' in options: seed = int(options['toydata_seed'])
        result = {'type': 'model_source', 'name': 'source', 'model': '@model', 'rnd_gen': {'seed': seed}}
        if input_spec.startswith('toys-asimov:'):
            result['dice_poisson'] = False
            result['dice_template_uncertainties'] = False
        return (result, theta_interface.delta_distribution(beta_signal = beta_signal_value))
    else: raise RuntimeError, 'data_source dict: not implemented for input = "%s"' % input_spec

# given (x,y) points of a function find the x value where the maximum is taken, suitable
# especially for binned data.
# xs must be monotonically increasing.
# Uses quadratic interpolation with the three points around the maximum to make the result more accurate.
# If the interpolated maximum is found outside the x range, the returned x value will be cut off at "half a binwidth", i.e. at
#         xs[0] - 0.5 * (xs[1] - xs[0])       and        xs[-1] + 0.5 * (xs[-1] - xs[-2])
# Only works well for non-noisy data.
def argmax_xy(xs, ys):
    assert len(xs) == len(ys)
    assert len(xs) >= 3
    maxy = -float("inf")
    maxi = 0
    for i in range(len(xs)):
        if ys[i] > maxy:
            maxy = ys[i]
            maxi = i
    if maxi==0: imin, imax = 0, 3
    elif maxi==len(xs)-1: imin, imax = len(xs)-3, len(xs)
    else: imin, imax = maxi-1, maxi+2
    xleft, x0, xright = xs[imin:imax]
    yleft, y0, yright = ys[imin:imax]
    a,b,c = numpy.linalg.solve([[xleft**2, xleft, 1], [x0**2, x0, 1], [xright**2, xright, 1]], [yleft, y0, yright])
    if a>=0: return xs[maxi]
    return -b/(2*a)
    
# same as argmax_xy, but use an instance of poltutil.plotdata as argument
def argmax(pd): return argmax_xy(pd.x, pd.y)

def extract_number(s):
    r = re.compile('(\d+)')
    m = r.search(s)
    if m is None: return None
    return float(m.group(1))


def get_p(n, n0):
    p = n*1.0 / n0
    p_error = max(math.sqrt(p*(1-p) / n0), 1.0 / n0)
    return p, p_error


def info(s):
    if not config.suppress_info:
        print "[INFO] ", s

def warning(s):
    print "[WARN] ", s

# return a list of result rows for the given query on the .db filename.
def sql_singlefile(filename, query, return_colnames = False):
    if not os.path.exists(filename): raise RuntimeError, "sql: the file %s does not exist!" % filename
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    try:  c.execute(query)
    except Exception, ex:
        print "exception executing %s on file %s: %s" % (query, filename, str(ex))
        raise ex
    result = c.fetchall()
    if return_colnames:
        desc = c.description
    c.close()
    conn.close()
    if return_colnames: return result, desc
    return result

def sql(filename_pattern, query):
    result = []
    for f in glob.glob(filename_pattern):
        result.extend(sql_singlefile(f, query))
    return result

def plotdata_from_histoColumn(h):
    a = array.array('d')
    a.fromstring(h)
    pdata = plotutil.plotdata()
    xmin = a[0]
    xmax = a[1]
    nbins = len(a) - 4
    binwidth = (xmax - xmin) / nbins
    pdata.x = [xmin + i*binwidth for i in range(nbins)]
    pdata.y = a[3:-1]
    pdata.color = '#000000'
    return pdata

