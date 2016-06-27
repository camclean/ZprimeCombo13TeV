from scipy.stats import *
from math import *

# The probability to observe >= n events for a Poisson
# distribution with mean mu ("ge" = "greater or equal")
def poisson_p_ge(n, mu): return poisson.sf(n-0.5, mu)

# The probability to observe <= n events for a Poisson
# distribution with mean mu ("le" = "less or equal")
def poisson_p_le(n, mu): return poisson.cdf(n, mu)

# the probability to observe exactly n events for a Poisson
# distribution with mean mu ("eq" = "equal")
def poisson_p_eq(n, mu):
    if mu==0.0: return 0
    return poisson.pmf(n, mu)


# convert a p-value to a Z-value using the "one-tailed normal distribution" convention
def p_to_Z(p): return -norm.ppf(p)

# plot the value in data as histogram in the given filename.
# If xmin and xmax is are None, they are set to the minimum and xmaxmum value of the data.
def plot_histogram(data, fname = 'histogram.pdf', nbins = 100, xmin = None, xmax = None, xlabel = '$x$', **plotargs):
    if xmin is None: xmin = min(data)
    if xmax is None: xmax = max(data)
    pd = plotdata(color = '#ff0000')
    pd.histogram(data, xmin, xmax, nbins)
    plot(pd, xlabel, '$N$', fname, **plotargs)

# plot x,y data and save in the given filename.
def plot_xy(x, y, fname = 'xy.pdf', xlabel = '$x$', ylabel = '$y$', **plotargs):
    assert len(x)==len(y), "x and y must be of same length!"
    pd = plotdata(as_function = True, color = '#ff0000')
    pd.x = x
    pd.y = y
    plot(pd, xlabel, ylabel, fname, **plotargs)

# plot x,y,ye data, where ye is the list of errors for y.
def plot_xye(x, y, ye, fname = 'xey.pdf', xlabel = '$x$', ylabel = '$y$', **plotargs):
    assert len(x)==len(y), "x and y must be of same length!"
    assert len(y)==len(ye), "y and ye must be of same length!"
    pd = plotdata(as_function = True, color = '#ff0000')
    pd.x = x
    pd.y = y
    pd.yerrors = ye
    plot(pd, xlabel, ylabel, fname, **plotargs)

# Get likelihood ratio test statistic values for an ensemble of toys generated according to 
# the background-only hypothesis mu = 0. This test statistic is defined as
#   t = log(    max_mu,theta  L(mu, theta|data) /   max_theta L(mu=0, theta|data)    )
# where in the numerator, the maximum is taken over all parameter values while in the denominator,
# the signal strength mu is fixed to 0.
#
# Note: this routine assumes that there is only one signal process.
#
# The return value is a list of floats, t-values.
def get_bkg_t(model, n = 1000):
    res = deltanll(model, 'toys:0.0', n)
    # in case of only one signal, directly return the dnll values:
    return res[res.keys()[0]]['dnll']
    
# get the likelihood ratio test statistic value for data; see get_bkg_lr for more information.
# The result is a float, the observed t-value.
def get_data_t(model):
    res = deltanll(model, 'data', 1)
    return res[res.keys()[0]]['dnll'][0]

# small wrapper for the CLs limit calculation in theta. Also assumes that there is only one signal.
def get_observed_cls_limit(model, **cls_args):
    exp, obs = cls_limits(model, cls_options = {'expected_bands': 0, 'reltol_limit': 0.02})
    return obs.y[0], obs.yerrors[0]

# get the posterior (x,y) data from theta. Note that the posterior values are NOT normalized!
def get_posterior(model, nscan = 100, mu_min = 0.0, mu_max = 3.0, **args):
    res = bayesian_posteriors(model, 'data', 1, {'beta_signal': (nscan, mu_min, mu_max)})
    h = res[res.keys()[0]]['beta_signal'][0]
    return [h.get_x_low(ibin) for ibin in range(h.get_nbins())], h.get_values()
    
# build the shape model from the lecture, i.e. a model with 30 bins with exponentially
# falling background with 1000. events and gaussian signal at the given mass.
#
# b_rate_unc is the (relative) uncertainty on the background rate
# if b_shape unc is set to True, it implements the shape uncerainty as discussed in the lecture.
def build_shape_model(signal_mass = 500., b_rate_unc = 0.0, b_shape_unc = False, seed = 50):
    #seed = 50  makes a nice 3.1sigma effect for mass = 500.
    nbins = 30
    model = Model()
    hf_b = HistogramFunction()
    xmin, xmax = 0.0, 1000.0
    h_b = exp_histo(xmin, xmax, nbins, 1000.0)
    hf_b.set_nominal_histo(h_b)
    if b_shape_unc:
        bplus = exp_histo(xmin, xmax, nbins, 1100., 0.002*1.2, 0.2 * 1.2)
        bminus = exp_histo(xmin, xmax, nbins, 900., 0.002*0.7, 0.2 * 0.6)
        hf_b.set_syst_histos('bshape', bplus, bminus)
        model.distribution.set_distribution('bshape', typ='gauss', mean = 0.0, width = 1.0, range = [-inf, inf])
    model.set_histogram_function('obs', 'b', hf_b)
    if b_rate_unc > 0.0:
        model.add_lognormal_uncertainty('brate', math.log(1 + b_rate_unc), 'b')
    # build the signal template:
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo(gauss_histo(xmin, xmax, nbins, 100.0, signal_mass, 150.))
    model.set_histogram_function('obs', 's%d' % int(signal_mass), hf_s)
    model.set_signal_processes('s*')
    # generate ranomd data according to background-only:
    scipy.random.seed(seed)
    data = [float(x) for x in scipy.random.poisson(h_b.get_values())]
    model.set_data_histogram('obs', Histogram(xmin, xmax, data))
    return model
    
# helper routine for build_modelS1: return a Histogram with nevents events following
# the form
#   exp(-lmbda*c) + c
def exp_histo(xmin, xmax, nbins, nevents, lmbda = 0.002, c = 0.2):
    bs = (xmax - xmin) / nbins
    data = [math.exp(-lmbda*x)+c for x in [xmin + (i+0.5)*bs for i in range(nbins)]]
    result = Histogram(xmin, xmax, data)
    return result.scale(nevents / result.get_value_sum())

# helper routine for build_modelS1: return a Histogram with nevents in form of a normal
# distribution around mu with standard deviation sigma
def gauss_histo(xmin, xmax, nbins, nevents, mu, sigma):
    bs = (xmax - xmin) / nbins
    data = [math.exp(-((x-mu)/sigma)**2 / 2.) for x in [xmin + (i+0.5)*bs for i in range(nbins)]]
    result = Histogram(xmin, xmax, data)
    return result.scale(nevents / result.get_value_sum())
    
    
def mle_print(model, input, n, **args):
    # see online documentation for mle http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/theta-auto/likelihood.html#maximum-likelihood-estimate
    # and explanation of return value structure: http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/theta-auto/examples.html#maximum-likelihood-fits
    options = Options()
    options.set('minimizer', 'strategy', 'newton_vanilla')
    result = mle(model, input, n, options = options, **args)
    spgs = args.get('signal_process_groups', model.signal_process_groups)
    for sp in spgs:
        print "MLE result for signal process group '%s' (%s)" % (sp, spgs[sp])
        n = None
        for p in model.get_parameters(spgs[sp]):
            if n is None: n = len(result[sp][p])
            print "%20s =" % p,
            for i in range(min([n, 10])):
                print  " %5.2f +- %5.2f" % (result[sp][p][i][0], result[sp][p][i][1]),
            print ""
    return result

            

            