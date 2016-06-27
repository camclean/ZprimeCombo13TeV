import config, utils, os.path, datetime, math, bisect

from theta_interface import *
import Report
import bayesian

from utils import *

from scipy import optimize


# get x value at y, using linear interpolation.
# y must be between x1 and x2
def get_x(x1, x2, y1, y2, y):
    if y==y1: return x1
    if y==y2: return x2
    # assert y is between y1 and y2:
    assert (y1 - y) * (y2 - y) < 0
    f = (y - y1) / (y2 - y1)
    return x1 + f * (x2 - x1)


class Belt:
    def __init__(self):
        self.truth_to_ts = {}
        self.sorted = False
        
    def sort(self):
        if self.sorted: return
        for t in self.truth_to_ts:
            self.truth_to_ts[t] = sorted(self.truth_to_ts[t])
        self.sorted = True
    
    def add_point(self, truth_value, ts_value):
        if truth_value not in self.truth_to_ts:
            self.truth_to_ts[truth_value] = []
        self.truth_to_ts[truth_value].append(ts_value)
        self.sorted = False
        # note: sorting later is much faster than keeping the list sorted by something like this:
        #bisect.insort_left(self.truth_to_ts[truth_value], ts_value)
    
    def merge(self, other_belt):
        self.sorted = False
        for t in other_belt.truth_to_ts:
            if t in self.truth_to_ts: self.truth_to_ts[t].extend(other_belt.truth_to_ts[t])
            else: self.truth_to_ts[t] = other_belt.truth_to_ts[:]
    
    # returns a tuple (cls upper limit, uncertainty)
    # uncertainty is from finite number of toys only (!)
    def cls(self, ts_value, cl = 0.95):
        self.sort()
        index_bkg = bisect.bisect_left(self.truth_to_ts[0.0], ts_value)
        # 'omp' = one minus p
        n_bkg = len(self.truth_to_ts[0.0])
        omp_bkg = index_bkg * 1.0 / n_bkg
        index_bkg_uncertainty = math.sqrt(omp_bkg * (1 - omp_bkg) * n_bkg)
        omp_bkg_plus = (index_bkg + index_bkg_uncertainty) * 1.0 / n_bkg
        omp_bkg_minus = (index_bkg - index_bkg_uncertainty) * 1.0 / n_bkg
        #print '1 - p_bkg = ', omp_bkg, omp_bkg_plus, omp_bkg_minus
        cls_values = []
        cls_uncertainties = []
        truth_values = sorted(self.truth_to_ts.keys())
        for t in truth_values:
            index = bisect.bisect_left(self.truth_to_ts[t], ts_value)
            n = len(self.truth_to_ts[t])
            omp = index * 1.0 / n
            index_uncertainty = math.sqrt(omp * (1 - omp) * n)
            omp_plus = (index + index_uncertainty) * 1.0 / n
            omp_minus = (index - index_uncertainty) * 1.0 / n
            #omp_uncertainties.append(omp_uncertainty)
            # the relative uncertainties add ...
            cls = omp / omp_bkg
            cls_plus = omp_plus / omp_bkg_minus
            cls_minus = omp_minus / omp_bkg_plus
            cls_values.append(cls)
            cls_uncertainties.append(max(cls_plus - cls, cls - cls_minus))
        if cls_values[0] - cls_uncertainties[0] <= 1 - cl or cls_values[-1] + cls_uncertainties[-1] > 1 - cl:
            print 'belt inversion: did not find intersection (ts=%f, max truth=%f); probably bet_signal scan range was not large enough; returning inf +- inf as limit!' % (ts_value, truth_values[-1])
            return float("inf"), float("inf")
        i = 0
        while i < len(cls_values) and cls_values[i] > 1 - cl: i+=1
        assert i>0 and i<len(cls_values)
        limit = get_x(truth_values[i-1], truth_values[i], cls_values[i-1], cls_values[i], 1 - cl)
        
        values_plus = [cls_values[i] + cls_uncertainties[i] for i in range(len(cls_values))]
        i = 0
        while i < len(truth_values) and values_plus[i] > 1 - cl: i+=1
        assert i>0 and i<len(truth_values)
        limit_plus = get_x(truth_values[i-1], truth_values[i], values_plus[i-1], values_plus[i], 1 - cl)
        
        values_minus = [cls_values[i] - cls_uncertainties[i] for i in range(len(cls_values))]
        i = 0
        while i < len(truth_values) and values_minus[i] > 1 - cl: i+=1
        assert i>0 and i<len(truth_values)
        limit_minus = get_x(truth_values[i-1], truth_values[i], values_minus[i-1], values_minus[i], 1 - cl)
        limit_uncertainty = max(abs(limit - limit_plus), abs(limit - limit_minus))
        #print 'limit for ts=%f: %f +- %f' % (ts_value, limit, limit_uncertainty)
        return limit, limit_uncertainty



# get the chi2 value for the given vector. cov_inv is the inverse of the covariance
# matrix
def chi2(vec, cov_inv):
    n = len(vec)
    s = 0
    for i in range(n):
        for j in range(n):
            s += vec[i] * cov_inv[i][j] * vec[j]
    return s
    

# return A, lambda, chi2 parameters
def exp_interpolation(xs, ys, y_covariance):
    lmbda0 = math.log(ys[-1] / ys[0]) / (xs[-1] - xs[0])
    #print lmbda0
    n = len(xs)
    i = n / 2
    A0 = math.exp(xs[i] * lmbda0) * ys[i]
    fitfunc = lambda p, x: A0 * p[0] * math.exp(p[1] * lmbda0 * x)
    cov_inv = numpy.linalg.inv(y_covariance)
    # the cost function to minimize:
    costfunc = lambda p: chi2([fitfunc(p, xs[i]) - ys[i] for i in range(n)], cov_inv)
    p0 = [1.0, 1.0]
    #print 'cost before: ', costfunc(p0)
    xopt, fopt, grad, bopt, fcalls, gradcalls, warnflags = optimize.fmin_bfgs(costfunc, p0, disp=0, full_output=1)
    #print xopt
    #print res
    #print 'cost after: ', costfunc(res)
    #print bopt
    print grad
    print bopt
    assert bopt[0][0] >= 0.0
    bopt = [[A0**2 * bopt[0][0], A0*lmbda0 * bopt[0][1]], [A0*lmbda0*bopt[1][0], lmbda0**2 * bopt[1][1]] ]
    result = {'A': xopt[0] * A0, 'lambda': xopt[1] * lmbda0, 'chi2': fopt, 'bopt': bopt}
    return result
    
    
def test_interpolation():
    A, lmbda = 2.0, -0.01
    n = 10
    y_width = 0.1
    xs = numpy.linspace(0, 10, num=n)
    ys = numpy.array(map(lambda x: A*math.exp(lmbda*x) * max(0.001, numpy.random.normal(1, y_width)), xs))
    y_cov = numpy.ndarray((n,n))
    for i in range(n):
        for j in range(n):
            y_cov[i][j] = (ys[i] * y_width)**2 if i==j else 0.0
    res = exp_interpolation(xs, ys, y_cov)
    #print res['chi2']
    return res

def test_i2():
    results = [test_interpolation() for i in range(100)]
    print 'chi2 mean width: ', get_mean_width([res['chi2'] for res in results])
    print 'A mean width: ', get_mean_width([res['A'] for res in results])
    print 'lambda mean width: ', get_mean_width([res['lambda'] for res in results])
    # calculate pulls for A and lambda:
    print 'A uncertainty: ', get_mean_width([math.sqrt(res['bopt'][0][0] * 2) for res in results])
    print 'A pull: ', get_mean_width([(res['A'] - 2.0) / math.sqrt(res['bopt'][0][0] * 2) for res in results])
    #print 'lambda pull: ', get_mean_width([res['lambda'] for res in results])
    #print 'bopt[0]: ', results[0]['bopt']
    #pd = plotutil.plotdata()
    #pd.histogram(values, 0.0, 20.0, 20)
    #pd.color = '#ff0000'
    #plot((pd,), "$\\chi^2$", "N", "chi2.pdf")

#test_i2()
#test_interpolation()


def submatrix(M, i0, i1):
    result = []
    for i in range(i0, i1):
        result.append([])
        for j in range(i0, i1):
            result[-1].append(M[i][j])
    return result
    
# make an exponential interpolation to a falling spectrum to get an estimate for the x value at y0.
# Does not necessarily use all x / y values, but usually only values next to y0.
#
# returns a dictionary with these values:
# * x0: estimate of x0
# * x0_error: error estimate for x0
# * chi2: chi2 value for the fit (will be zero for two points included ...)
# * A, lambda: the parameters of the fit function (f(x) = A * exp(lambda * x))
def exp_interpolation_smart(xs, ys, y_covariance, y0):
    assert len(xs) == len(ys) and len(ys) > 1
    i_intersect = 0
    n = len(xs)
    while i_intersect < n and ys[i_intersect] > y0: i_intersect += 1
    i_first = max(0, i_intersect - 1)
    i_last1 = min(n, i_intersect + 1)
    if i_first + 1 == i_last1:
        if i_first == 0: i_last += 1
        else: i_first -= 1
    res0 = exp_interpolation(xs[i_first:i_last1], ys[i_first:i_last1], submatrix(y_covariance, i_first, i_last1))
    
    
    
    



# belts is a dictionary (signal process name) -> (Belt instance) as returned by prepare_belts
def plot_belts(belts, cl = 0.95):
    plotsdir = os.path.join(config.workdir, 'plots')
    config.report.new_section('plot_belts')
    colors = ['#edd400', '#f57900', '#c17d11', '#73d216', '#3465a4', '#75507b', '#d3d7cf', '#555753']
    # for each sp, make a plot different ts values with x axis = truth, y-axis = CLs value
    for sp in belts:
        name = sp
        config.report.add_html('<h2>signal process "%s"</h2>' % sp)
        belt = belts[sp]
        belt.sort()
        truth_values = sorted(belt.truth_to_ts.keys())
        # as "resonable" values for the ts value, use the values at truth=0.0 which would enter band construction
        at = lambda frac, l: l[int(frac * len(l))]
        ts0 = belt.truth_to_ts[truth_values[0]]
        ts_values = (at(0.025, ts0), at(0.16, ts0),  at(0.5, ts0), at(0.84, ts0),at(0.975, ts0))
        # for each ts value, plot CLs versus truth:
        # returns the fraction of x within the (sorted) list l
        frac = lambda x, l: bisect.bisect_left(l, x) * 1.0 / len(l)
        #binom_error  = lambda p, n: math.sqrt(p * (1-p) / n)
        binom_error  = lambda p, n: math.sqrt(p * (1-p) / n + 1.0 / (4 * n**2))
        frac_error = lambda x,l: (frac(x,l), binom_error(frac(x,l), len(l)))
        icol = -1
        pds = []
        for ts in ts_values:
            icol = (icol + 1) % len(colors)
            pd = plotutil.plotdata()
            pd.x = truth_values[:]
            pd.y = []
            pd.yerrors = []
            pd.lw = 1.0
            omp0, omp0_error = frac_error(ts, ts0)
            for truth in pd.x:
                 omp, omp_error = frac_error(ts, belt.truth_to_ts[truth])
                 pd.y.append(omp / omp0)
                 pd.yerrors.append(math.sqrt(omp_error **2 / omp0**2 + omp**2 / omp0**4 * omp0_error**2))
            pd.color = colors[icol]
            pd.as_function = True
            pd.draw_line = False
            pds.append(pd)
            # try to fit with a power law (x - mu)^{-k}, but only use  few points to the left, not all;
            # for now, use (apart from the point immidiately to the left), 2 more points to the left,
            # and all points to the right
            i_intersect = 0
            while i_intersect < len(pd.y) and pd.y[i_intersect] > 1-cl: i_intersect += 1
            if i_intersect == len(pd.y):
                print 'did not find index for cl level; beta range too small!'
                continue
            i_first = max(0, i_intersect - 3)
            i_first_plot = max(0, i_intersect - 1)
            i_first_plot = i_first
            i_last1 = min(len(pd.y), i_intersect + 3)
            lmbda = math.log(pd.y[i_last1 - 1] / pd.y[i_first]) / (pd.x[i_last1 - 1] - pd.x[i_first])
            A = math.exp(-pd.x[i_intersect] * lmbda) * pd.y[i_intersect]
            fitfunc = lambda p, x: A * p[0] * math.exp(p[1] * lmbda * x)
            errfunc = lambda p, x, y, yerrors: [(fitfunc(p,x[i]) - y[i]) / yerrors[i] for i in range(len(y))]
            p0 = [1.0, 1.0]
            args = (pd.x[i_first:i_last1],pd.y[i_first:i_last1], pd.yerrors[i_first:i_last1])
            p1, err = optimize.leastsq(errfunc, p0[:], args = args)
            print p1, err
            print errfunc(p1, *args)
            #p1, cov_p1 = optimize.curve_fit(fitfunc, pd.x, pd.y, p0, pd.yerrors)
            #print out
            #print p1, success
            pd2 = plotutil.plotdata()
            pd2.x = numpy.linspace(pd.x[i_first_plot], pd.x[i_last1-1], num=100)
            pd2.y  = map(lambda x: fitfunc(p1, x), pd2.x)
            pd2.color = pd.color
            pd2.as_function = True
            pds.append(pd2)
            
        plot(pds, 'truth value', 'CLs value', os.path.join(plotsdir, 'plot_belts-%s.png' % name))
        config.report.add_html('<p><img src="plots/plot_belts-%s.png" /></p>' % name)


'''
        ts_start, ts_stop = ts_values[0] - 0.5 * abs(ts_values[0] - ts_values[-1]), ts_values[-1] + 0.5 * abs(ts_values[0] - ts_values[-1])
        pd = plotutil.plotdata()
        pd.x = numpy.linspace(ts_start, ts_stop, num=100)
        pd.y = []
        for i in range(len(pd.x)-1):
            pd.y.append(len([x for x in ts0 if x >= pd.x[i] and x < pd.x[i+1]]))
        i = len(pd.x)-1
        pd.y.append(len([x for x in ts0 if x >= pd.x[i]]))
        pd.color = '#000000'
        plot((pd,), 'ts value', '$N_{toys}$', os.path.join(plotsdir, 'plot_belts-ts0-%s.png' % name))
        config.report.add_html('<p><img src="plots/plot_belts-ts0-%s.png" /></p>' % name)
        # returns the fraction of x within the (sorted) list l
        frac = lambda x, l: bisect.bisect_left(l, x) * 1.0 / len(l)
        binom_error  = lambda p, n: math.sqrt(p * (1-p) / n)
        frac_error = lambda x,l: (frac(x,l), binom_error(frac(x,l), len(l)))
        icol = 0
        pds = []
        for ts in ts_values:
            pd = plotutil.plotdata()
            pd.x = truth_values[:]
            pd.y = []
            pd.yerrors = []
            pd.lw = 1.0
            omp0, omp0_error = frac_error(ts, ts0)
            for truth in pd.x:
                 omp, omp_error = frac_error(ts, belt.truth_to_ts[truth])
                 pd.y.append(omp / omp0)
                 pd.yerrors.append(math.sqrt(omp_error **2 / omp0**2 + omp**2 / omp0**4 * omp0_error**2))
            pd.color = colors[icol]
            pd.as_function = True
            pd.draw_line = False
            pds.append(pd)
            # try to fit:
            fitfunc = lambda p, x: 1 / (1 + math.exp(p[0]**2 * (x - p[1])))
            errfunc = lambda p, x, y: [fitfunc(p,x[i]) - y[i] for i in range(len(y))]
            p0 = [1., 5.]
            p1, success = optimize.leastsq(errfunc, p0[:], args=(pd.x,pd.y))
            print p1, success
            pd2 = plotutil.plotdata()
            pd2.x = numpy.linspace(min(pd.x), max(pd.x), num=100)
            pd2.y  = map(lambda x: fitfunc(p1, x), pd2.x)
            pd2.color = pd.color
            pd2.as_function = True
            pds.append(pd2)
            icol = (icol + 1) % len(colors)

        plot(pds, 'ts value', 'CLs value', os.path.join(plotsdir, 'plot_belts-%s.png' % name))
        config.report.add_html('<p><img src="plots/plot_belts-%s.png" /></p>' % name)
'''

# ts is either 'lr' for likelihood ratio or 'mle' for maximum likelihood estimate.
#
# for 'lr', there are different variants: (note: the names are taken from
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit )
# * 'LEP' will do a 'LEP-like' likelihood ratio in which no minimization is done at all; the likelihood ratio is evaluated using the
#    most probable values for the nuisance parameters and beta_signal=1, beta_signal=0.
#    this can be achieved by
#    nuisance_prior = 'shape:fix;rate:fix'
#    signal_prior = 'fix:1'
# * 'Tevatron' minimizes the likelihood ratio w.r.t. all nuisance parameters and use beta_signal=1, beta_signal=0
#    This can be achieved by using
#     nuisance_prior = ''
#     signal_prior = 'fix:1'
# * 'ATLAS' method is similar to 'Tevatron' but will also let beta_signal free in the numerator of the likelihood ratio (and use beta_signal=0 for
#    the denominator). This is achieved by
#     nuisance_prior = ''
#     signal_prior = 'flat'
#
# Additional variants exist, for example the default is similar to 'ATLAS', however, shape-changing parameters are fixed.
#
# options:
# * 'beta_signal_ranges': a dictionary (signal process) -> (min_beta, max_beta) to use for the scan; default is to run Bayesian limits on background only and use 
#    [0, median  + 5 * width]  as range. Note that a min_beta of 0.0 must always be included if you evaluate the belts with cls.
#
# the signal_prior setting is the likelihood constraint added for minimization while calculating the likelihood ratio
# or mle. For the likelihood ratio, it is only used to calculate the maximised likelihood value for the signal+background
# hypothesis. (For background only, beta_signal = 0 is always fixed to zero.)
#
# returns: a dictionary (signal process group id) -> Belt instance (see above)
def prepare_belts(model, ts = 'lr', signal_prior='flat', nuisance_prior = 'shape:fix', beta_scan_n = 20, n_toys_per_beta = 10000, signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    beta_signal_ranges = {}
    # 1. find a good range for the beta_signal-scan by running Bayesian limits on background only
    if 'beta_signal_ranges' not in options:
        quant = bayesian.bayesian_quantiles(model, input = 'toys:0', n = 30, signal_prior = 'flat', nuisance_prior = '', quantile = 0.95, write_report = False, signal_processes = signal_processes)
        for sp in quant:
            limits = sorted(quant[sp])
            n = len(limits)
            assert n > 20, 'too many failures in calculating Bayesian quantiles'
            median, upper1s, upper2s = limits[int(0.5 * n)], limits[int(0.84 * n)], limits[int(0.975 * n)]
            beta_max = median + max(5 * (upper1s - median), 2.5 * (upper2s - median))
            utils.info("prepare_belts: For signal process '%s', using beta_max=%f" % (sp, beta_max))
            beta_signal_ranges[sp] = (0.0, beta_max)
    else:
        beta_signal_ranges = options['beta_signal_ranges']
        for sp in signal_processes:
            sp_id = ''.join(sp)
            assert sp_id in beta_signal_ranges, 'the option beta_signal_ranges does not provide a beta range for signal processes group "%s"!' % sp_id
    # 2. use these beta_signal_max for scanning:
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': beta_scan_n * n_toys_per_beta, 'model': '@model', 'producers': ('@ts_producer',),
        'output_database': sqlite_database(), 'log-report': False}
    minimizer = {'type': 'root_minuit'}
    if ts == 'mle':
        ts_producer = {'type': 'mle', 'name': 'mle', 'minimizer': minimizer, 'parameter': 'beta_signal',
           'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior")}
        ts_colname = 'mle__beta_signal'
    elif ts == 'lr':
        ts_producer = {'type': 'deltanll_hypotest', 'name': 'lr', 'minimizer': minimizer,
        'signal-plus-background-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
        'background-only-distribution': product_distribution(delta_distribution(beta_signal = 0.0), "@nuisance_prior")}
        ts_colname = 'lr__nll_diff'
    else: raise RuntimeError, 'unknown ts "%s"' % ts
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer}    
    toplevel_settings.update(get_common_toplevel_settings(**options))
    seed = 1
    if 'toydata_seed' in options: seed = int(options['toydata_seed'])
    main['data_source'] = {'type': 'model_source', 'name': 'source', 'model': '@model', 'rnd_gen': {'seed': seed}}
    cfg_names_to_run = []
    for sp in signal_processes:
        sp_id = ''.join(sp)
        beta_signal_range = beta_signal_ranges[sp_id]
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        toplevel_settings['model-distribution-signal'] = equidistant_deltas('beta_signal', beta_signal_range, beta_scan_n)
        name = write_cfg(model, sp, 'prepare_belts', 'scan%g' % beta_max, additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    
    # 3. build belts
    result = {}
    cachedir = os.path.join(config.workdir, 'cache')
    for name in cfg_names_to_run:
        method, sp_id, dummy = name.split('-',2)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select source__beta_signal, %s from products' % ts_colname)
        #dt_start = datetime.datetime.now()
        belt = Belt()
        for beta, ts in data:
            belt.add_point(beta, ts)
        belt.sort()
        #dt_end = datetime.datetime.now()
        #print dt_start, dt_end
        result[sp_id] = belt
    return result


# calculates CLs limits on the specified input ensemble.
# 
# options are passed to prepare_belts, including 'nuisance_prior', 'signal_prior' and 'ts' which define the 
# test statistic. To control belt construction range and statistics, use 'beta_scan_n', 'n_toy_per_beta' and 'beta_signal_ranges'
#
# returns: nothing so far, use the html report for results.
def cls_limits(model, input = 'toys:0', n = 1000, cl = 0.95, ts = 'lr', signal_prior = 'flat', nuisance_prior = 'shape:fix', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    # 1. make the belts:
    options['ts'] = ts
    options['signal_prior'] = signal_prior
    options['nuisance_prior'] = nuisance_prior
    belts = prepare_belts(model, **options)
    # 2. calculate the test statistic values for the ensemble 'input'
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ('@ts_producer',),
        'output_database': sqlite_database(), 'log-report': False}
    minimizer = {'type': 'root_minuit'}
    if ts == 'mle':
        ts_producer = {'type': 'mle', 'name': 'mle', 'minimizer': minimizer, 'parameter': 'beta_signal',
           'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior")}
        ts_colname = 'mle__beta_signal'
    elif ts == 'lr':
        ts_producer = {'type': 'deltanll_hypotest', 'name': 'lr', 'minimizer': minimizer,
        'signal-plus-background-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
        'background-only-distribution': product_distribution(delta_distribution(beta_signal = 0.0), "@nuisance_prior")}
        ts_colname = 'lr__nll_diff'
    else: raise RuntimeError, 'unknown ts "%s"' % ts
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer}
    toplevel_settings.update(get_common_toplevel_settings(**options))
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input, **options)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'cls_limits', '', additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    
    cachedir = os.path.join(config.workdir, 'cache')
    
    result = {}
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    if input=='data': header = '%f %% C.L. upper limit' % (cl * 100)
    else: header = '%f %% C.L. upper limit (median; central 68%%; central 95%%)' % (cl * 100)
    result_table.add_column('limit', header)
    for name in cfg_names_to_run:
        method, sp, dummy = name.split('-',2)
        result_table.set_column('process', sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select %s from products' % ts_colname)
        #TODO: make belt evaluation more efficient !!!
        #result[sp] = map(lambda ts: belts[sp].cls(ts, cl)[0], [row[0] for row in data])
        ts_values = sorted([row[0] for row in data])
        n = len(ts_values)
        if n == 0:
           result_table.set_column('%f quantile' % q, 'N/A')
           continue
        if input == 'data':
            limit, limit_uncertainty = belts[sp].cls(ts_values[n/2], cl)
            result_table.set_column('limit', '%.5g +- %.3g' % (limit, limit_uncertainty))
        else:
            limit_median, limit1low, limit1high, limit2low, limit2high = map(lambda ts: belts[sp].cls(ts, cl)[0], [ts_values[n / 2], ts_values[int(0.16 * n)],
                ts_values[int(0.84 * n)], ts_values[int(0.025 * n)], ts_values[int(0.975 * n)]])
            result_table.set_column('limit', '%.3g  (%.3g, %.3g; %.3g, %.3g)' % (limit_median, limit1low, limit1high, limit2low, limit2high))
        result_table.add_row()
    config.report.new_section("CLs limits on ensemble '%s'" % input)
    #config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_prior)))
    config.report.add_html(result_table.html())
    return result

