#!/usr/bin/python
# coding=utf8

# This script produces the result pdf plots from the *.db files.
# go to the very end of this file for some more annotations of which
# plot is produced where.

import matplotlib

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.text
import matplotlib.lines
import matplotlib.patches
import numpy

import scipy.special
import scipy.optimize as optimize
from scipy import interpolate

import sqlite3, math, os, os.path

# save the median Z_est value for the signal plus background model here:
median_z_est = 0.0

# return a list of result rows for the given query on the .db filename.
def sql(filename, query):
    if not os.path.exists(filename): raise RuntimeError, "sql: the file %s does not exist!" % filename
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    c.execute(query)
    result = c.fetchall()
    c.close()
    conn.close()
    return result
                            
def add_xlabel(axes, text, *args, **kwargs):
    label = axes.set_xlabel(text, size='large', ha='right', *args, **kwargs)
    label.set_position((1.0, 0.03))
    return label

def add_ylabel(axes, text, *args, **kwargs):
    label = axes.set_ylabel(text, size='large', va='top', *args, **kwargs)
    label.set_position((-0.03, 1.0))
    return label

def add_lumititle(ax):
    title = ax.set_title('$L = 200\,\mathrm{pb}^{-1}$', ha='right')
    title.set_position((1.0, 1.0))

#add secondary title:
def add_stitle(ax, title):
    return ax.text(0.0, 1.03, title, transform = ax.transAxes, ha='left', va='bottom')

class plotdata:
    def __init__(self):
        self.x = None
        self.y = None
        self.legend = None
        self.yerrors = None
        self.xerrors = None
        self.fill_xrange = None
        self.fill_color = None
        self.fill_hatch = None
        self.color = None
        self.marker = 'None'
        self.lw = 2
        self.extra_plot_args = {}
        # an array of bands; a band is a three-tuple (y1, y2, color). y1 and y2 are 
        # arrays of y values.
        # bands to draw first should come first 
        self.bands = None

    def set_from_nphisto(self, histo):
        self.x = histo[1]
        self.y = histo[0]
        #shift x values:
        self.xerrors = [0.5 * (x2 - x1) for x1, x2 in zip(self.x[:-1], self.x[1:])]
        self.x = [0.5 * (x1 + x2) for x1, x2 in zip(self.x[:-1], self.x[1:])]
        self.yerrors = map(math.sqrt, self.y)
        

#histos is a list of tuples (label, histogram)
def plot(histos, xlabel, ylabel, outname, log = False, logx = False, legend_args = {}, ax_modifier=None, norm=False, stitle=None, xmax=None, ymin=None, ymax=None):
    if norm: map(lambda h: h[1].Scale(1.0/h[1].Integral()), histos)
    cm = 1.0/2.54
    fsize = 15*cm, 12*cm
    fp = fm.FontProperties(size=10)
    fig = plt.figure(figsize = fsize)
    rect = 0.15, 0.15, 0.8, 0.8
    ax = fig.add_axes(rect)
    if log:  ax.set_yscale('log')
    if logx: ax.set_xscale('log')
    #else: ax.set_ylim(ymin=0)
    add_xlabel(ax, xlabel, fontproperties=fp)
    add_ylabel(ax, ylabel, fontproperties=fp)
    if type(stitle)==type(""):
        add_stitle(ax, stitle)
    if xmax!=None:
        ax.set_xlim(xmax=xmax)
    if ymin!=None:
        ax.set_ylim(ymin=ymin)
    
    if ymax!=None:
        ax.set_ylim(ymax=ymax)
    draw_legend = False
    for histo in histos:
        if histo.legend: draw_legend = True
        if histo.bands is not None:
            for band in histo.bands:
                ax.fill_between(histo.x, band[0], band[1], lw=0, facecolor=band[2])
        if histo.yerrors is not None: ax.errorbar(histo.x, histo.y, histo.yerrors, xerr = histo.xerrors, elinewidth=1, lw=0, label=histo.legend, color=histo.color)
        else: ax.plot(histo.x, histo.y, lw=histo.lw, label=histo.legend, color=histo.color, marker=histo.marker, *histo.extra_plot_args)
        if histo.fill_xrange is not None:
            where = [(x > histo.fill_xrange[0] and x < histo.fill_xrange[1]) for x in histo.x]
            color = histo.fill_color
            #ax.fill_between(histo.x, histo.y, y2 = [ax.get_ylim()[0]]*len(histo.y), where = where, lw=0, facecolor=color)
            poly_xy = [(x,y) for (x,y) in zip(histo.x, histo.y) if x >= histo.fill_xrange[0] and x <= histo.fill_xrange[1]]
            poly_x = [x for (x,y) in poly_xy]
            poly_y = [y for (x,y) in poly_xy]
            poly_x = [poly_x[0]] + poly_x + [poly_x[-1]]
            poly_y = [ax.get_ylim()[0]] + poly_y + [ax.get_ylim()[0]]
            ax.fill(poly_x, poly_y, facecolor=color, hatch=histo.fill_hatch, lw=0)
    if draw_legend: ax.legend(prop=fp,**legend_args)
    if ax.get_legend() is not None:
        map(lambda line: line.set_lw(1.5), ax.get_legend().get_lines())
    
    if ymin!=None:
        ax.set_ylim(ymin=ymin)
    
    if ax_modifier!=None: ax_modifier(ax)
    fig.savefig(outname)
    del fig

def p_to_Z(p):
    return math.sqrt(2) * scipy.special.erfinv(1 - 2*p)

def deltanll():
    z_values = sql("results/deltanll_hypo.db", "select hypotest__nll_sb - hypotest__nll_b from products")
    z_values = map(lambda x: math.sqrt(2*abs(x[0])), z_values)
    z_plot = plotdata()
    z_plot.set_from_nphisto(numpy.histogram(z_values, bins = 50, range = (0.0, 5.0)))
    z_plot.color = '#00aa00'
    z_plot.xerrors = None
    plot((z_plot,), '$Z_{\mathrm{est}}$', '$N_{\mathrm{PE}}$', 'results/deltanll_hypo.pdf')
    global median_z_est
    median_z_est = numpy.median(z_values)
    print "deltanll: median Z_est = %f" % median_z_est

def deltanll2():
    z_values_bonly = sql("results/deltanll_hypo_bonly.db", "select hypotest__nll_sb - hypotest__nll_b from products")
    z_values_bonly = map(lambda x: math.sqrt(2*abs(x[0])), z_values_bonly)
    z_plot_bonly = plotdata()
    z_plot_bonly.legend = "B only hypothesis"
    z_plot_bonly.color = '#0000aa'
    z_plot_bonly.set_from_nphisto(numpy.histogram(z_values_bonly, bins = 100, range = (0.0, 5.0)))
    
    z_values_sb = sql("results/deltanll_hypo.db", "select hypotest__nll_sb - hypotest__nll_b from products")
    z_values_sb = map(lambda x: math.sqrt(2*abs(x[0])), z_values_sb)
    z_plot_sb = plotdata()
    z_plot_sb.legend = "S + B hypothesis"
    z_plot_sb.color = '#00aa00'
    z_plot_sb.set_from_nphisto(numpy.histogram(z_values_sb, bins = 100, range = (0.0, 5.0)))
    
    z_plot_bonly.fill_xrange = (median_z_est, 999.0)
    
    # determine the expected p value.
    # a. all "background only" events:
    N = len(z_values_bonly)
    # b. count "background only" events above the median_z_est for "signal + background": 
#    global median_z_est
    Np = len(filter(lambda x: x > median_z_est, z_values_bonly)) * 1.0
    # c. p-value is now very easy; convert it to a Z value:
    p = Np / N
    print "deltanll2: expected Z = %f" % p_to_Z(p)
    
    median_line = lambda ax: (ax.add_line(matplotlib.lines.Line2D([median_z_est, median_z_est], [1.0, 10000.], lw=1.0, color=z_plot_sb.color, ls='-', drawstyle='default')),
         ax.add_artist(matplotlib.text.Text(median_z_est, 10000 * 1.1, "median of S + B", horizontalalignment = "center", color=z_plot_sb.color, size='smaller')),
         ax.add_artist(matplotlib.text.Text(3.5, 80, "$N_{\mathrm{PE}} \cdot \hat p$", color=z_plot_bonly.color, size='smaller'))
          )
    plot((z_plot_bonly,z_plot_sb), '$Z_{\mathrm{est}}$', '$N_{\mathrm{PE}}$ per bin', 'results/deltanll_hypo2.pdf', log = True, ymin = 1.0, ax_modifier=median_line)
    

def bayesfactor():
    bf = sql('results/bayesfactor.db', 'select bayesfactor__nl_posterior_sb - bayesfactor__nl_posterior_b from products;')
    bf = map(lambda x: math.exp(x[0]), bf)
    median = numpy.median(bf)
    print "bayesfactor: median bayes factor: %g" % median
    bf_plot = plotdata()
    bf_plot.color = '#0000aa'
    bf_plot.set_from_nphisto(numpy.histogram(bf, bins = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100]))
    plot((bf_plot,), '$B_{01}$', '$N_{PE}$', 'results/bf.pdf', logx = True, ax_modifier = (lambda ax: ax.set_xlim(0.99e-4)))
    
def profile_likelihood_intervals():
    data = sql('results/profile_likelihood_intervals.db', 'select intervals__maxl, intervals__lower06800, intervals__upper06800, intervals__lower09500, intervals__upper09500 from products')
    coverage68 = len([1 for x in data if (10.0 > x[1] and x[2] > 10.0)]) * 1.0 / len(data)
    coverage95 = len([1 for x in data if (10.0 > x[3] and x[4] > 10.0)]) * 1.0 / len(data)
    print "profile_likelihood_intervals: coverage 68: %f, 95: %f" % (coverage68, coverage95)
    mu_s_median = numpy.median([x[0] for x in data])
    print "profile_likelihood_intervals: median mu_s: ", mu_s_median
    one_sigma_lower_median = numpy.median([x[1] for x in data])
    one_sigma_upper_median = numpy.median([x[2] for x in data])
    print "profile_likliehood_intervals: median interval: ", one_sigma_lower_median, one_sigma_upper_median
    
def neyman(data, nbins, rangex, xlabel, ylabel, filename, median_ts):
    mu_s_binborders = numpy.linspace(rangex[0], rangex[1], nbins + 1)
    binned_mu_s = [[] for i in mu_s_binborders[1:]]
    for (x,y) in data:
        ix = int((x - rangex[0]) / (rangex[1] - rangex[0]) * nbins)
        if ix<0 or ix>=nbins: continue
        binned_mu_s[ix] += [y]
    map(lambda ylist: ylist.sort(), binned_mu_s)
    pdata = plotdata()
    pdata.color = '#000000'
    pdata.marker = '+'
    pdata.lw = 0.5
    pdata.x = mu_s_binborders
    pdata.x = [0.5 * (x1 + x2) for x1, x2 in zip(pdata.x[:-1], pdata.x[1:])]
    pdata.y = [ylist[len(ylist)/2] for ylist in binned_mu_s]
    yhigh_1s = [ylist[int(len(ylist) * 0.84)] for ylist in binned_mu_s]
    ylow_1s = [ylist[int(len(ylist) * 0.16)] for ylist in binned_mu_s]
    yhigh_2s = [ylist[int(len(ylist) * 0.975)] for ylist in binned_mu_s]
    ylow_2s = [ylist[int(len(ylist) * 0.025)] for ylist in binned_mu_s]
    
    ymid_in = scipy.interpolate.interp1d(pdata.x, pdata.y)
    x = optimize.fsolve(lambda x: ymid_in(x) - median_ts, 10.0)
    
    yhigh_in = scipy.interpolate.interp1d(pdata.x, yhigh_1s)
    lowx = optimize.fsolve(lambda x: yhigh_in(x) - median_ts, 10.0)
    
    ylow_in = scipy.interpolate.interp1d(pdata.x, ylow_1s)
    highx = optimize.fsolve(lambda x: ylow_in(x) - median_ts, 10.0)
    
    print "68%% interval for %s: (%f, %f), 0%% interval: %f" % (filename, lowx, highx, x)
    
    pdata.bands = [(ylow_2s, yhigh_2s, '#e4e400'), (ylow_1s, yhigh_1s, '#00e400')]
    band1_l = matplotlib.patches.Rectangle((0,0), 1, 1, fc='#00e400')
    band2_l = matplotlib.patches.Rectangle((0,0), 1, 1, fc='#e4e400')
    
    interval_drawer = lambda ax: (ax.add_line(matplotlib.lines.Line2D([0.001, highx], [median_ts, median_ts], lw=1.0, color='#ff0000', ls='-')), 
      ax.add_line(matplotlib.lines.Line2D([lowx, lowx], [0.001, median_ts], lw=1.0, color='#ff0000', ls=':')),
      ax.add_line(matplotlib.lines.Line2D([highx, highx], [0.001, median_ts], lw=1.0, color='#ff0000', ls=':')),
      ax.add_line(matplotlib.lines.Line2D([x, x], [0.001, median_ts], lw=1.0, color='#ff0000', ls=':')),
      ax.add_artist(matplotlib.text.Text(0, median_ts, ' median '+ylabel, size='small', color='#ff0000')),
      ax.legend((band1_l, band2_l), ('68% central belt', '95% central belt'), loc='upper left'))
    plot((pdata,), xlabel, ylabel, filename, ax_modifier=interval_drawer)


def neyman_z_est():
    data = sql('results/neyman_z_est.db', 'select data_source__mu_s, hypotest__nll_b - hypotest__nll_sb from products')
    data = map(lambda x: (x[0], math.sqrt(2.0*abs(x[1]))) , data)
    nbins = 30
    rangex = (0.0, 30.0)
    ts = [d[1] for d in data if abs(d[0] - 10.0)< 0.2]
    median = numpy.median(ts)
    print "neyman_z_est: median TS is ", median
    neyman(data, nbins, rangex, '$\mu_s$', r'$Z_{\mathrm{est}}$', 'results/neyman_z_est.pdf', median)

def neyman_mle():
    data = sql('results/neyman_mle.db', 'select data_source__mu_s, mle__mu_s from products')
    nbins = 30
    rangex = (0.0, 30.0)
    ts = [d[1] for d in data if abs(d[0] - 10.0)< 0.2]
    median = numpy.median(ts)
    print "neyman_mle: median TS is ", median
    neyman(data, nbins, rangex, '$\mu_s$', '$\hat{\mu}_s$', 'results/neyman_mle.pdf', median)

def neyman_mle_syst():
    data = sql('results/neyman_mle_syst.db', 'select data_source__mu_s, mle__mu_s from products')
    nbins = 30
    rangex = (0.0, 30.0)
    ts = [d[1] for d in data if abs(d[0] - 10.0)< 0.2]
    median = numpy.median(ts)
    print "neyman_mle_syst: median TS is ", median
    neyman(data, nbins, rangex, '$\mu_s$', '$\hat{\mu}_s$', 'results/neyman_mle_syst.pdf', median)
    
def neyman_mle_syst_interp():
    data = sql('results/neyman_mle_syst_interp.db', 'select data_source__mu_s, mle__mu_s from products')
    nbins = 30
    rangex = (0.0, 30.0)
    ts = [d[1] for d in data if abs(d[0] - 10.0)< 0.2]
    neyman(data, nbins, rangex, '$\mu_s$', '$\hat{\mu}_s$', 'results/neyman_mle_syst_interp.pdf', numpy.median(ts))
    
def cls_z_est():
    data = sql('results/neyman_z_est.db', 'select data_source__mu_s, hypotest__nll_b - hypotest__nll_sb from products')
    data0 = sql('results/deltanll_hypo_bonly.db', 'select data_source__mu_s, hypotest__nll_b - hypotest__nll_sb from products')
    data = data + data0
    data = map(lambda x: (x[0], math.sqrt(2.0*abs(x[1]))) , data)
    nbins = 30
    rangex = (0.0, 30.0)
    mu_s_binborders = numpy.linspace(rangex[0], rangex[1], nbins + 1)
    binned_mu_s = [[] for i in mu_s_binborders[1:]]
    for (x,y) in data:
        ix = int((x - rangex[0]) / (rangex[1] - rangex[0]) * nbins)
        if ix<0 or ix>=nbins: continue
        binned_mu_s[ix] += [y]
    map(lambda ylist: ylist.sort(), binned_mu_s)
    pdata = plotdata()
    pdata.color = '#000000'
    pdata.x = mu_s_binborders
    pdata.x = [0.5 * (x1 + x2) for x1, x2 in zip(pdata.x[:-1], pdata.x[1:])]
    p_b = lambda ts_crit: len([1 for ts in binned_mu_s[0] if ts > ts_crit]) * 1.0 / len(binned_mu_s[0])
    cls = lambda ts_list, ts_crit, p_b__ts_crit: (1.0 - len([1 for ts in ts_list if ts > ts_crit])*1.0 / len(ts_list)) / (1.0 - p_b__ts_crit)
    
    def calc_upper(measured_ts = median_z_est, cl = 0.95):
        p_b__ts_crit = p_b(measured_ts)
        y = map(lambda ts_list: cls(ts_list, measured_ts, p_b__ts_crit), binned_mu_s)
        # use linear interpolation to find out mu_s for which CLs = 0.05:
        cls_inter = scipy.interpolate.interp1d(pdata.x, y)
        x0 = scipy.optimize.fsolve(lambda x: cls_inter(x) - (1 - cl), 7)
        return x0
    print "cls_z_est: expected 95%% upper limit in case of signal is:", calc_upper(median_z_est)
    print "cls_z_est: expected 95%% upper limit in case of no signal is:", calc_upper(0.01)
    
def posterior():
    data = sql('results/posterior.db', 'select post__posterior_mu_s from products limit 0, 1')
    import array
    a = array.array('d')
    a.fromstring(data[0][0])
    pdata = plotdata()
    xmin = a[0]
    xmax = a[1]
    nbins = len(a) - 4
    binwidth = (xmax - xmin) / nbins
    pdata.x = numpy.linspace(xmin, xmax, nbins + 1)
    #delete overflow, underflow and range:
    a = a[3:-1]
    nbins = len(a)
    pdata.set_from_nphisto((a, pdata.x))
    pdata.y = numpy.array(pdata.y) * (1.0 / sum(pdata.y) / binwidth)
    pdata.yerrors = None
    pdata.xerrors = None
    pdata.color = '#0000ff'
    
    cum = pdata.y.cumsum()
    cum *= binwidth
    f = interpolate.interp1d(pdata.x, cum)
    xlow = optimize.fsolve(lambda x: f(x) - 0.16, 5.0)
    xhigh = optimize.fsolve(lambda x: f(x) - 0.84, 15.0)
    print "posterior: 68%% credible interval is %f -- %f" %(xlow, xhigh)
    
    fy = interpolate.interp1d(pdata.x, pdata.y)
    interval_drawer = lambda ax: (ax.add_line(matplotlib.lines.Line2D([xlow, xlow], [0., fy(xlow)], lw=1.0, color='#ff0000', ls='-')),
      ax.add_line(matplotlib.lines.Line2D([xhigh, xhigh], [0., fy(xhigh)], lw=1.0, color='#ff0000', ls='-')),
      ax.add_artist(matplotlib.text.Text(0.5*(xlow + xhigh), 0.62 * ax.get_ylim()[1], '68%', color='#ff0000', verticalalignment='center', horizontalalignment='center'))
      )
    
    p = [1,10,5]
    fitrange_bins = [0,nbins-1]
    while pdata.y[fitrange_bins[0]] < 0.06: fitrange_bins[0] += 1
    while pdata.y[fitrange_bins[1]] < 0.06: fitrange_bins[1] -= 1

	# fit a normal distribution around the maximum to determine the position of the maximum:    
    func = lambda p, x: p[0] * math.exp(-(x - p[1])**2 / p[2])
    fitfunc = lambda p, x: numpy.array(map(lambda xx: func(p, xx), x))
    errfunc = lambda p, x, y: (fitfunc(p, x) - y)[fitrange_bins[0]:fitrange_bins[1]]
    p1, success = optimize.leastsq(errfunc, p[:], args=(pdata.x, pdata.y))
    print "posterior: found argmax is at mu_s = %f " % p1[1]
    
    pdata_fit = plotdata()
    pdata_fit.x = pdata.x[fitrange_bins[0]:fitrange_bins[1]]
    pdata_fit.y = fitfunc(p1, pdata_fit.x)
    pdata_fit.color = '#ff0000'
    
    #to draw the fit, uncomment the following line:
    #plot((pdata,pdata_fit), "$\mu_s$", "$p(\mu_s)$", "results/posterior.pdf", xmax = 30.0, ax_modifier = interval_drawer)
    plot((pdata,), "$\mu_s$", "$p(\mu_s)$", "results/posterior.pdf", xmax = 30.0, ax_modifier = interval_drawer)


def cls_illustration():

    gauss = lambda mean, width, norm, x: 1.0 / math.sqrt(2*math.pi*width) * math.exp(-(x-mean)**2 / (2*width**2))
    ts_b = lambda x: gauss(10, 3, 1, x)
    pd_b = plotdata()
    pd_b.x = numpy.linspace(0, 30, 301)
    pd_b.y = map(ts_b, pd_b.x)
    pd_b.color = '#0000aa'
    pd_b.legend = "B hypothesis"
    pd_b.fill_xrange = (0, 17.)
    pd_b.fill_color = pd_b.color
    
    ts_sb = lambda x: gauss(18, 4, 1, x)
    pd_sb = plotdata()
    pd_sb.x = pd_b.x[:]
    pd_sb.y = map(ts_sb, pd_sb.x)
    pd_sb.color = '#00aa00'
    pd_sb.legend = "S + B hypothesis"
    pd_sb.fill_xrange = (17., 99.)
    pd_sb.fill_color = pd_sb.color
    
    draw_text = lambda ax: (ax.add_artist(matplotlib.text.Text(5, 0.185, '$1-p_{b}$', color=pd_b.color)),
       ax.add_artist(matplotlib.text.Text(20, 0.185, '$p_{s+b}$', color=pd_sb.color)),
       ax.add_artist(matplotlib.text.Text(17, 0.21, '$\hat{T}$', horizontalalignment='center')),
       ax.add_line(matplotlib.lines.Line2D([17, 17], [0, 0.206], lw=1, color='black'))
       )
    
    plot((pd_b,pd_sb), "$TS$", "$d(TS)$", "results/cls_illustration.pdf", ax_modifier=draw_text)
    


# Figure 2: Distribution of Z_est for S+B.
deltanll()

#Figure 3: Distribution of Z_est for S+B (as in Figure 2) and high stat. for Background only;
# calculate Z values via tail of B only.
deltanll2()

# Figure 4: expected Bayes factors for S+B
bayesfactor()

#(no figure): profile likelihood method
profile_likelihood_intervals()

# Figure 5: Neyman construction using Z_est as test statistic:
neyman_z_est()

# Figure 6: Neyman construction using maxmimum likelihood estimate of mu_s as test statistic:
neyman_mle()

# Figure 7: posterior for mu_s
posterior()

# Figure 8: illustration of the definitions used in the CLs method
cls_illustration()

#(no figure): CLs method based on Z_est as test statistic
cls_z_est()

# Figue 9: Neyman construction using maximum likelihood estimate as test statistic, including an uncertainty on tau
neyman_mle_syst()

# Figue 11: Neyman construction using maximum likelihood estimate as test statistic, including an uncertainty on the signal shape
neyman_mle_syst_interp()

