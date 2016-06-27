# -*- coding: utf-8 -*-
# coding=utf8

import matplotlib
# note: some matplotlib backends are broken and cause segfaults in fig.save.
# try commenting out and in the use of Cairo in case of problems ...
#try:
#matplotlib.use('Cairo')
#except ImportError: pass

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.text
import matplotlib.lines
import matplotlib.patches

from scipy import interpolate
import numpy as np

import math, StringIO

import theta_auto

def add_xlabel(axes, text, *args, **kwargs):
    label = axes.set_xlabel(text, size='large', ha='right', *args, **kwargs)
    label.set_position((1.0, 0.03))
    return label

def add_ylabel(axes, text, *args, **kwargs):
    label = axes.set_ylabel(text, size='large', va='top', *args, **kwargs)
    label.set_position((-0.03, 1.0))
    return label

class plotdata:
    """
    Class holding (x,y) data and formatting information for plotting (1D-)histograms and functions.
    
    :ivar x: array of x values. For Histograms, these are the lower bin borders.
    :ivar xmax: The maximum x value for Histograms.
    :ivar y: array of y values
    :ivar yerrors: array of errors in y or ``None`` if there are no y errors
    :ivar xerrors: array of x errors or ``None`` if there are no x errors
    :ivar bands: array of bands to draw. A single band is a thre-tuple ``(ymin, ymax, color)`` where ``ymin`` and ``ymax`` are arrays of y values and color is a matplotlib color directive, e.g. "#rrggbb".
    
    :ivar as_function: ``True`` if the data should be displayed as function. Otherwise, it will be drawn as histogram
    :ivar color: The color used for drawing the line
    :ivar fmt: matplotlib format string for the line. Default is "-", i.e., a solid line.
    :ivar lw: line width for the line
    :ivar legend: The legend string or ``None`` in case no legend should be drawn
    :ivar fill_color: color to be used for filling. The default ``None`` does not draw a filled area
    :ivar fill_xrange: a tuple ``(xmin, xmax)`` which should be filled with fill_color. Currently only implemented for histograms.
    :ivar yerrors_mode: a string specifying how errors are displayed. The default "bars" will draw error bars. Currently, the only other option is "area" which draws a shaded area instead, which is useful if the errors are so small (or points so dense in x-direction) that the error bars would overlap.
    :ivar capsize: cap size to use for y error bars (default: 1.5)
    """
    
    def __init__(self, color = '#000000', legend = None, as_function = False, lw = 2, legend_order = 0):
        self.x = []
        self.y = []
        self.xmax = None # required only for one-binned histograms ...
        self.legend = legend
        self.legend_order = legend_order
        self.yerrors = None
        # capsize for the error bars:
        self.capsize = 1.5
        # how to draw the yerrors; valid are "bars" for error bars (using color, capsize and lw), 'bars0' to also show the marker for y values of 0,
        # or 'area' to draw only a shaded area.
        self.yerrors_mode = 'bars'
        self.yerrors_fill_alpha = 0.5
        self.xerrors = None
        # filling will be done only is fill_color is not None
        self.fill_color = None
        # in case filling is done, the x value range to fill. The default "None"
        # is to fill everything (i.e., whole x-range). Otherwise, specify tuple (xmin, xmax) here.
        self.fill_xrange = None
        self.fill_to_y = 0.0
        self.color = color
        self.marker = 'None'
        self.markersize = 1.0
        self.lw = lw
        self.fmt = '-'
        # an array of bands; a band is a three-tuple (y1, y2, color). y1 and y2 are 
        # arrays of y values.
        # bands to draw first should come first 
        self.bands = None
        self.band_lw = 0
        self.bands_fill = True
        self.as_function = as_function
        self.draw_histo = True
        self.draw_line = True
        
    
    def histogram(self, values, xmin, xmax, nbins, errors = False, include_uoflow = False):
        """
        create new data by making a histogram from xmin to xmax with nbins from the given values
        which should be a iterable yielding floats.
        If errors is True, yerrors is set to sqrt(n) in each bin.
        if include_uoflow is True, values under (over) the range are inserted in the first (last) bin.
        """
        xmin, xmax = float(xmin), float(xmax)
        self.xmax = xmax
        self.x = [xmin + (xmax - xmin) / nbins * i for i in range(nbins)]
        self.y = [0.0] * nbins
        for v in values:
            ibin = int((v - xmin) / (xmax - xmin) * nbins)
            if not include_uoflow:
                if ibin < 0 or ibin >= nbins: continue
            else:
                if ibin < 0: ibin = 0
                if ibin >= nbins: ibin = nbins-1
            self.y[ibin] += 1
        if errors: self.yerrors = map(math.sqrt, self.y)
        
    # add the values in values to the current histogram, i.e. increments y[ibin] by 1
    # for the corresponding bin for each v in values.
    def histogram_add(self, values, include_uoflow = False):
        xmin, xmax, nbins = self.x[0], self.xmax, len(self.x)
        if len(self.y) == 0: self.y = [0.0] * nbins
        assert len(self.y) == nbins
        for v in values:
            ibin = int((v - xmin) / (xmax - xmin) * nbins)
            if not include_uoflow:
                if ibin < 0 or ibin >= nbins: continue
            else:
                if ibin < 0: ibin = 0
                if ibin >= nbins: ibin = nbins-1
            self.y[ibin] += 1
        

    def scale_y(self, factor):
        """
        Scale all y values by the given factor. Also scales the y-values in the bands; y errors are not scaled.
        """
        self.y = [y*factor for y in self.y]
        if self.bands is None: return
        for band in self.bands:
            band[0][:] = [y * factor for y in band[0]]
            band[1][:] = [y * factor for y in band[1]]

    def scale_x(self, factor):
        """ Scale all x-values by the given factor, including xmax.
        """
        self.x = [x*factor for x in self.x]
        self.xmax *= factor

    # set data according to the "histo triple" h = (xmin, xmax, data) or Histogram instance
    def histo_triple(self, h):
        binwidth = (h[1] - h[0]) / len(h[2])
        self.x = [h[0] + i * binwidth for i in range(len(h[2]))]
        self.y = h[2][:]
        self.xmax = h[1]
        
    def set_histogram(self, histo):
        """
        Set the variables x, xmax, y, yerrors from ``histo``. ``histo`` should be an instance of :class:`theta_auto.Histogram`.
        """
        self.xmax = histo.get_xmax()
        self.x = [histo.get_x_low(i) for i in range(histo.get_nbins())]
        self.y = histo.get_values()
        self.yerrors = histo.get_uncertainties()
        self.xerrors = None
    
    
    def histo(self):
        """
        Return current x,y data as :class:`theta_auto.Histogram` instance
        """
        ye = self.yerrors[:] if self.yerrors is not None else None
        # note: copy all data!
        return theta_auto.Histogram(self.x[0], self.xmax, self.y[:], ye, x_low = self.x[:])
    
    
    # replace x, y and bands by a smoothed version, obtained by cubic interpolation
    # evaluated n times more points than original
    #
    # s is the amount of smoothing (default None is the scipy default)
    #
    # relunc is the relative uncertainty for each y value to assume for s
    #
    # in case more control is needed, make the smoothing outside and re-set the x,y values.
    #
    # should only be used in cases self.as_function = True
    def smooth(self, n = 3, s = None, relunc = 0.05, miny_factor = 0.0):
        oldx = self.x[:]
        # assume a 5% uncertainty for the smoothing
        y_average = sum(self.y) / len(self.y)
        tck = interpolate.splrep(oldx, self.y, w = [1 / (relunc * max(self.y[i], y_average * miny_factor)) for i in range(len(oldx))], s = s)
        self.x = list(np.linspace(min(self.x), max(self.x), n * len(self.x)))
        self.y = interpolate.splev(self.x, tck)
        if self.bands is None: return
        for band in self.bands:
            tck = interpolate.splrep(oldx, band[0], w = [1 / (relunc * band[0][i]) for i in range(len(oldx))], s = s)
            band[0][:] = interpolate.splev(self.x, tck)
            tck = interpolate.splrep(oldx, band[1], w = [1 / (relunc * band[1][i]) for i in range(len(oldx))], s = s)
            band[1][:] = interpolate.splev(self.x, tck)
        
    def write_txt(self, ofile):
        """
        Write the content in a text file.
        
        :param ofile: The output file, either as string or an open file handle
        
        One line is written per (x,y) point. The first line is a comment line starting with "#" explaining the fields. The general data layout is::
        
          x y yerror band0ymin band0ymax band1ymin band1ymax ...
          
        where in general some entries are missing if not available.
        """
        if type(ofile)==str: ofile = open(ofile, 'w')
        ofile.write('# x; y')
        if self.yerrors is not None: ofile.write('; yerror')
        if self.bands is not None:
            for k in range(len(self.bands)):
                ofile.write('; band %d low; band %d high' % (k, k))
        ofile.write("\n")
        for i in range(len(self.x)):
            ofile.write("%10.5g %10.5g " % (self.x[i], self.y[i]))
            if self.yerrors is not None:
                ofile.write("%10.5g " % self.yerrors[i])
            if self.bands is not None:
                for k in range(len(self.bands)):
                    ofile.write("%10.5g %10.5g" % (self.bands[k][0][i], self.bands[k][1][i]))
            ofile.write("\n")
            
    # printing support, similar to write_txt, but in a string. Note that this might cxhange in the future; it's mainly for
    # quick inspection in the script by the user. Use write_txt to be safe.
    def __str__(self):
        sio = StringIO.StringIO()
        self.write_txt(sio)
        return sio.getvalue()

    def __repr__(self): return str(self)
    
    def read_txt(self, infile):
        """
        Read data from a file produced by :meth:`write_txt`. This replaces the instance variables x, y, yerrors and bands.
        """
        values = []
        for line in file(infile):
            if len(line) == 0: continue
            if line[0] == '#': continue
            line_values = map(lambda s: float(s), line.split())
            # check that the number of values in the lines agree:
            if len(values) > 0:
                if len(values[0]) != len(line_values): raise RuntimeError, "number of values given is inconsistent!"
            values.append(line_values)
        n_values = len(values[0])
        have_yerrors = (n_values % 2 == 1)
        # read x, y values:
        self.x = [row[0] for row in values]
        self.y = [row[1] for row in values]
        if have_yerrors:
            self.yerrors = [row[2] for row in values]
        else:
            self.yerrors = None
        self.xerrors = None
        # read bands:
        n_bands = (n_values - 2) / 2
        self.bands = []
        colors = ['#ffff00', '#00ff00']
        yerror_offset = 0
        if have_yerrors: yerror_offset = 1
        for i in range(n_bands):
            band = ([row[2+2*i +  yerror_offset] for row in values], [row[3+2*i  + yerror_offset] for row in values], colors[i % len(colors)])
            self.bands.append(band)
        
def plot(histos, xlabel, ylabel, outname = None, logy = False, logx = False, ax_modifier=None, title_ul=None, title_ur = None,
 extra_legend_items = [], xmin = None, xmax=None, ymin=None, ymax=None, legend_args = {}, fig = None, figsize_cm = (15, 12), fontsize = 10, axes_creation_args = {}):
    """
    Plot the given :class:`plotutil.plotdata` objects. Many drawing options are controlled by those instances; see documentation there.
    
    :param histos: A list of :class:`plotutil.plotdata` instances or a single instance
    :param xlabel: The label for the x axis. Latex is allowed in $$-signs
    :param ylabel: The label for the y axis. Latex is allowed in $$-signs
    :param outname: name of the output file; the file extension will be used to guess the file type (by matplotlib); typical choices are ".pdf" and ".png".
    :param logy: use log-scale in y direction
    :param logx: use log scale in x direction
    :param ax_modifier: function called with the matplotlib.Axes object as argument. Allows to perform additional manipulation in case you need to "extend" this method
    :param title_ul: Title for the upper left corner. Can be either a string or a list of strings. A list of strings will be displayed as multiple lines.
    :param title_ur: Title for the upper right corner. Unlike ``title_ul``, only a single string is allowed
    :param extra_legend_items: allows to specify extra items for the legend. It is a list of two-tuples ``(handle, legend)`` where ``handle`` is a matplotlib object to use to draw the legend and ``legend`` is the string to be drawn
    :param xmin: The minimum x value to draw
    :param xmax: The maximum x value to draw
    :param ymin: The minimum y value to draw
    :param ymax: The maximum y value to draw
    :param figsize_cm: The figure size in cm as tuple ``(width, height)``
    :param fontsize: The font size in points
    """
    #extra_legend_items = extra_legend_items[:]
    legend_items = []
    cm = 1.0/2.54
    fsize = figsize_cm[0]*cm, figsize_cm[1]*cm
    fp = fm.FontProperties(size = fontsize)
    if fig is None:
        fig = plt.figure(figsize = fsize)
    axes_creation_args = dict(axes_creation_args) # make copy to prevent next line from modifying argument.
    rect = axes_creation_args.pop('rect', (0.15, 0.15, 0.8, 0.75))
    ax = fig.add_axes(rect, **axes_creation_args)
    if logy: ax.set_yscale('log')
    if logx: ax.set_xscale('log')
    add_xlabel(ax, xlabel, fontproperties=fp)
    add_ylabel(ax, ylabel, fontproperties=fp)
    if title_ul is not None:
        if type(title_ul) == type([]):
            yoffset = 1.02
            for s in title_ul:
                ax.text(0.0 if yoffset>1 else 0.02, yoffset, s, transform = ax.transAxes, ha='left', va = 'bottom' if yoffset > 1 else 'top')
                yoffset -= fontsize * 1.0 / (72 * fsize[1]) * 1.5
        else:
            ax.text(0.0, 1.02, title_ul, transform = ax.transAxes, ha='left', va='bottom')
    if title_ur is not None: ax.text(1.0, 1.02, title_ur, transform = ax.transAxes, ha='right', va='bottom')
    draw_legend = False
    if isinstance(histos, plotdata): histos = [histos]
    for histo in histos:
        legend_added = False
        assert len(histo.x)==len(histo.y), "number of x,y coordinates not the same for '%s'" % histo.legend
        if histo.legend: draw_legend = True
        # allow empty "dummy" plots which have legend but no content:
        if len(histo.x)==0:
            if histo.fill_color is not None:
                legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw), histo.legend))
            else:
                legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, lw=histo.lw), histo.legend))
            continue
        if histo.bands is not None:
            for band in histo.bands:
                if histo.bands_fill:
                    if len(band) > 3: c = band[3]
                    else: c = band[2]
                    ax.fill_between(histo.x, band[0], band[1], lw=histo.band_lw, facecolor=band[2], color=c)
                else:
                    xs = histo.x + [x for x in reversed(histo.x)]
                    ys = band[0] + [y for y in reversed(band[1])]
                    xs.append(xs[0])
                    ys.append(ys[0])
                    ax.plot(xs, ys, lw=histo.band_lw, color=band[2])
        if not histo.as_function:
            if histo.yerrors is None or histo.draw_histo:
               new_x = [histo.x[0]]
               for x in histo.x[1:]: new_x += [x]*2
               new_x += [histo.xmax]
               new_y = []
               for y in histo.y: new_y += [y]*2
               if logy and ymin is not None:
                    for i in range(len(new_y)): new_y[i] = max(new_y[i], ymin)
               if histo.fill_color is not None:
                    if histo.fill_xrange is None:
                        ax.fill_between(new_x, new_y, [0] * len(new_y), lw=histo.lw, color=histo.color, facecolor = histo.fill_color)
                        if histo.legend is not None: legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw), histo.legend))
                    else:
                        x_clipped = [histo.fill_xrange[0]]
                        y_clipped = []
                        for x,y in zip(new_x, new_y):
                            if x >= histo.fill_xrange[0]:
                                if len(y_clipped)==0: y_clipped.append(y)
                                if x <= histo.fill_xrange[1]:
                                    x_clipped.append(x)
                                    y_clipped.append(y)
                                else:
                                    x_clipped.append(histo.fill_xrange[1])
                                    y_clipped.append(y_clipped[-1])
                                    break
                        ax.fill_between(x_clipped, y_clipped, [histo.fill_to_y] * len(y_clipped), lw=0, color=None, facecolor = histo.fill_color)
                        ax.plot(new_x, new_y, histo.fmt, lw=histo.lw, color=histo.color)
                        if histo.legend is not None: legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw), histo.legend))
                        #legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, lw=histo.lw), histo.legend))
               else:
                   ax.plot(new_x, new_y, histo.fmt, lw=histo.lw, color=histo.color)
                   if histo.legend is not None: legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, lw=histo.lw, ls = histo.fmt), histo.legend))
               legend_added = True
            # if histo.yerrors is set, draw with errorbars, shifted by 1/2 binwidth ...
            if histo.yerrors is not None:
                if histo.yerrors_mode.startswith('bars'):
                    if histo.xmax is None: raise RuntimeError, "need xmax in histogram for y errors"
                    low_x = histo.x + [histo.xmax]
                    x_centers = [0.5 * (low_x[i] + low_x[i+1]) for i in range(len(histo.x))]
                    ys = histo.y
                    yerrors = histo.yerrors
                    if histo.yerrors_mode != 'bars0':
                        x_new, y_new, ye_new = [], [], []
                        for x,y,ye in zip(x_centers, ys, yerrors):
                            if y != 0:
                                x_new.append(x)
                                y_new.append(y)
                                ye_new.append(ye)
                        x_centers, ys, yerrors = x_new, y_new, ye_new
                        if len(x_centers) == 0: continue # do not plot empty list at all
                    ax.plot(x_centers, ys, marker = histo.marker, markersize=histo.markersize, ls='None', mew = 0.0, mfc = histo.color)
                    ax.errorbar(x_centers, ys, yerrors, ecolor = histo.color, capsize = histo.capsize, lw = histo.lw, fmt = None)
                    if not legend_added:
                        if histo.legend is not None: legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, marker=histo.marker, markevery=(1,10), mew=0.0, markersize=histo.markersize, lw=histo.lw), histo.legend))
                        legend_added = True
                else:
                    new_x = [histo.x[0]]
                    for x in histo.x[1:]: new_x += [x]*2
                    new_x += [histo.xmax]
                    new_y_low, new_y_high = [], []
                    for y, ye in zip(histo.y, histo.yerrors):
                        new_y_low += [y - ye]*2
                        new_y_high += [y + ye]*2
                    ax.fill_between(new_x, new_y_high, new_y_low, lw=histo.lw, color = histo.color, facecolor = histo.fill_color, alpha = histo.yerrors_fill_alpha)
                    if not legend_added:
                        if histo.legend is not None: legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw, alpha = histo.yerrors_fill_alpha), histo.legend))
                        legend_added = True
        else: # i.e. as_function:
            if histo.yerrors is not None:
                if histo.yerrors_mode == 'bars':
                    lw = histo.lw
                    if histo.draw_line is False: lw = 0
                    larg = {}
                    if histo.legend is not None: larg = {'label': histo.legend}
                    ax.errorbar(histo.x, histo.y, histo.yerrors, elinewidth = histo.lw, lw=lw, color=histo.color, marker=histo.marker, markersize=histo.markersize, **larg)
                elif histo.yerrors_mode == 'area':
                    y_high = [y + yerror for (y, yerror) in zip(histo.y, histo.yerrors)]
                    y_low = [y - yerror for (y, yerror) in zip(histo.y, histo.yerrors)]
                    ax.fill_between(histo.x, y_high, y_low, lw = histo.lw, color = histo.color, facecolor = histo.fill_color, alpha = histo.yerrors_fill_alpha)
                else:
                    raise RuntimeError, "yerrors_mode='%s' for as_function=True not supported!" % histo.yerrors_mode
            else:
                if histo.fill_color is not None:
                    if histo.fill_xrange is not None:
                        x_draw, y_draw = [], []
                        for x,y in zip(histo.x, histo.y):
                            if x >= histo.fill_xrange[0] and x <= histo.fill_xrange[1]:
                                x_draw.append(x)
                                y_draw.append(y)
                    else:
                        x_draw, y_draw = pd.x, pd.y
                    ax.fill_between(x_draw, y_draw, [histo.fill_to_y] * len(y_draw), lw=histo.lw, label=histo.legend, color=histo.color, facecolor = histo.fill_color)
                else:
                    ax.plot(histo.x, histo.y, histo.fmt, lw=histo.lw, color=histo.color, marker=histo.marker, markersize=histo.markersize)
                    if histo.legend is not None: legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, ls = histo.fmt, lw=histo.lw), histo.legend))
    legend_items = [(h,l) for (o,h,l) in sorted(legend_items, cmp = lambda x,y: cmp(x[0], y[0]))]
    if draw_legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = handles[:]
        labels = labels[:]
        for h, l in  legend_items + extra_legend_items:
            labels.append(l)
            if type(h)==str:
                h = matplotlib.patches.Rectangle((0, 0), 1, 1, fc=h)
            handles.append(h)
        ax.legend(handles, labels, prop = fp, numpoints = 1, **legend_args)
    #if ax.get_legend() is not None:
    #    map(lambda line: line.set_lw(1.5), ax.get_legend().get_lines())

    if ymin!=None:
        ax.set_ylim(ymin=ymin)
    if ymax!=None:
        ax.set_ylim(ymax=ymax)
    if xmin!=None:
        ax.set_xlim(xmin=xmin)
    if xmax is None:
        if len(histos) > 0:
            xmax = max([max(pd.x) for pd in histos])
            xmax = max([xmax] + [pd.xmax for pd in histos if pd.xmax is not None])
    if xmax is not None:
        ax.set_xlim(xmax=xmax)
    
    if ax_modifier!=None: ax_modifier(ax)
    if outname is not None: fig.savefig(outname)
    del fig
    
def make_stack(pdatas):
    for i in range(len(pdatas)):
        for j in range(i+1, len(pdatas)):
            pdatas[i].y = map(lambda x: x[0] + x[1], zip(pdatas[i].y, pdatas[j].y))
            if pdatas[i].yerrors is not None and pdatas[j].yerrors is not None:
                pdatas[i].yerrors = map(lambda x: math.sqrt(x[0]**2 + x[1]**2), zip(pdatas[i].yerrors, pdatas[j].yerrors))


def scatter_ax_m(x, y, xycol = None, s = 8):
   if xycol is not None:
       if type(xycol)==str:
           color = xycol
       else:
           color = [xycol(xv, yv) for xv, yv in zip(x,y)]
   else:
       color = '#000000'
   return lambda ax: (ax.scatter(x,y, color = color, s = s))


