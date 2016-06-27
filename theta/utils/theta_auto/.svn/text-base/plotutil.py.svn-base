# -*- coding: utf-8 -*-
# coding=utf8

import matplotlib
# note: some matplotlib backends are broken and cause segfaults in fig.save.
# try commenting out and in the use of Cairo in case of problems ...
#matplotlib.use('Cairo')

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.text
import matplotlib.lines
import matplotlib.patches

from scipy import interpolate
import numpy as np

import math

def add_xlabel(axes, text, *args, **kwargs):
    label = axes.set_xlabel(text, size='large', ha='right', *args, **kwargs)
    label.set_position((1.0, 0.03))
    return label

def add_ylabel(axes, text, *args, **kwargs):
    label = axes.set_ylabel(text, size='large', va='top', *args, **kwargs)
    label.set_position((-0.03, 1.0))
    return label

# plotdata represents the data of a single curve in a plot, including drawing options, legend, etc.
class plotdata:
    def __init__(self, color = '#000000', legend = None, as_function = False, lw = 2):
        self.x = []
        self.y = []
        self.legend = legend
        self.yerrors = None
        self.xerrors = None
        self.fill_color = None
        self.color = color
        self.marker = 'None'
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
        
    # make a histogram of the given values
    def histogram(self, values, xmin, xmax, nbins, errors = False):
        xmin, xmax = float(xmin), float(xmax)
        self.x = [xmin + (xmax - xmin) / nbins * i for i in range(nbins)]
        self.y = [0.0] * nbins
        if errors: self.yerrors = [0.0] * nbins
        for v in values:
            ibin = int((v - xmin) / (xmax - xmin) * nbins)
            if ibin < 0 or ibin >= nbins: continue
            self.y[ibin] += 1
            if errors: self.yerrors[ibin] += 1
        if errors: self.yerrors = map(math.sqrt, self.yerrors)

    # scale all y values by factor
    def scale_y(self, factor):
        self.y = [y*factor for y in self.y]
        if self.bands is None: return
        for band in self.bands:
            band[0][:] = [y * factor for y in band[0]]
            band[1][:] = [y * factor for y in band[1]]

    # set data according to the "histo triple" h = (xmin, xmax, data)
    def histo_triple(self, h):
        binwidth = (h[1] - h[0]) / len(h[2])
        self.x = [h[0] + i * binwidth for i in range(len(h[2]))]
        self.y = h[2][:]
    
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
    def smooth(self, n = 3, s = None, relunc = 0.05):
        oldx = self.x[:]
        # assume a 5% uncertainty for the smoothing
        tck = interpolate.splrep(oldx, self.y, w = [1 / (relunc * self.y[i]) for i in range(len(oldx))], s = s)
        self.x = list(np.linspace(min(self.x), max(self.x), n * len(self.x)))
        self.y = interpolate.splev(self.x, tck)
        if self.bands is None: return
        for band in self.bands:
            tck = interpolate.splrep(oldx, band[0], w = [1 / (relunc * band[0][i]) for i in range(len(oldx))], s = s)
            band[0][:] = interpolate.splev(self.x, tck)
            tck = interpolate.splrep(oldx, band[1], w = [1 / (relunc * band[1][i]) for i in range(len(oldx))], s = s)
            band[1][:] = interpolate.splev(self.x, tck)
        
    # ofile is a string (filename) or a handle to an open file
    def write_txt(self, ofile):
        if type(ofile)==str: ofile = open(ofile, 'w')
        ofile.write('# x; y')
        if self.bands is not None:
            for k in range(len(self.bands)):
                ofile.write('; band %d low; band %d high' % (k, k))
        ofile.write("\n")
        for i in range(len(self.x)):
            ofile.write("%10.5g %10.5g " % (self.x[i], self.y[i]))
            if self.bands is not None:
                for k in range(len(self.bands)):
                    ofile.write("%10.5g %10.5g" % (self.bands[k][0][i], self.bands[k][1][i]))
            ofile.write("\n")
    
    # infile is the filename
    def read_txt(self, infile):
        # values is a list of lines in the file; each line is a list of floats
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
        assert n_values % 2 == 0, "invalid number of values (has to be even)"
        # read x, y values:
        self.x = [row[0] for row in values]
        self.y = [row[1] for row in values]
        # read bands:
        n_bands = (n_values - 2) / 2
        self.bands = []
        colors = ['#ffff00', '#00ff00']
        for i in range(n_bands):
            band = ([row[2+2*i] for row in values], [row[3+2*i] for row in values], colors[i % len(colors)])
            self.bands.append(band)
        
        

## \brief Make a plot and write it to an output file
#
#
# histos is a list / tuple of plotdata instances, xlabel and ylabel are lables for the axes, outname is the output filename.
#
# logx and logy control whether the x/y-scale should be logarithmic.
#
# If set, ax_modifier should be a function / callable. It will be called with the matplotlib Axes instance as only argument.
# This can be used for arbitrary mainulation, such as adding other objects to the plot, etc.
#
# title_ul and title_ur are strings shown on the upper left and upper right of the plot, resp.
#
# extra_legend_items is a list of extra items to add to the legend, as tuple (handle, lable) where handle is a matplotlib Artist
# and label the legend string. As special case, we also allow handle to be a string in which case it is assumed to encode a color.
# For other cases, see matplotlib legend guide for details.
#
# xmin, xmax, ymin, ymax control the region shown in the plot. the default is to determine the region automatically (by matplotblib)
def plot(histos, xlabel, ylabel, outname = None, logy = False, logx = False, ax_modifier=None, title_ul=None, title_ur = None,
extra_legend_items = [], xmin = None, xmax=None, ymin=None, ymax=None, legend_args = {}, fig = None):
    cm = 1.0/2.54
    fsize = 15*cm, 12*cm
    fp = fm.FontProperties(size=10)
    if fig is None:
        fig = plt.figure(figsize = fsize)
    ax = fig.add_axes((0.15, 0.15, 0.8, 0.75))
    if logy:  ax.set_yscale('log')
    if logx: ax.set_xscale('log')
    add_xlabel(ax, xlabel, fontproperties=fp)
    add_ylabel(ax, ylabel, fontproperties=fp)
    if title_ul is not None:
        if type(title_ul) == type([]):
            yoffset = 1.02
            for s in title_ul:
                ax.text(0.0 if yoffset>1 else 0.02, yoffset, s, transform = ax.transAxes, ha='left', va = 'bottom' if yoffset > 1 else 'top')
                yoffset -= 0.05
        else:
            ax.text(0.0, 1.02, title_ul, transform = ax.transAxes, ha='left', va='bottom')
    if title_ur is not None: ax.text(1.0, 1.02, title_ur, transform = ax.transAxes, ha='right', va='bottom')
    draw_legend = False
    for histo in histos:
        assert len(histo.x)==len(histo.y), "number of x,y coordinates not the same for '%s'" % histo.legend
        if histo.legend: draw_legend = True
        # allow empty "dummy" plots which have legend but no content:
        if len(histo.x)==0: continue
        if histo.bands is not None:
            for band in histo.bands:
                if histo.bands_fill:
                    ax.fill_between(histo.x, band[0], band[1], lw=histo.band_lw, facecolor=band[2], color=band[2])
                else:
                    xs = histo.x + [x for x in reversed(histo.x)]
                    ys = band[0] + [y for y in reversed(band[1])]
                    xs.append(xs[0])
                    ys.append(ys[0])
                    ax.plot(xs, ys, lw=histo.band_lw, color=band[2])
        if not histo.as_function:
            # histo.x is assumed to contain the lower bin edges in this case ...
            if len(histo.x) >= 2:  x_binwidth = histo.x[1] - histo.x[0]
            else: x_binwidth = 1.0
            # if histo.yerrors is set, draw with errorbars, shifted by 1/2 binwidth ...
            if histo.yerrors is not None:
               new_x = [x + 0.5 * x_binwidth for x in histo.x]
               ax.errorbar(new_x, histo.y, histo.yerrors, label=histo.legend, ecolor = histo.color, marker='o', ms = 2, capsize = 0, fmt = None)
            if histo.yerrors is None or histo.draw_histo:
               new_x = [histo.x[0]]
               for x in histo.x[1:]: new_x += [x]*2
               new_x += [histo.x[-1] + x_binwidth]
               new_y = []
               for y in histo.y: new_y += [y]*2
               if logy and ymin is not None:
                    for i in range(len(new_y)): new_y[i] = max(new_y[i], ymin)
               if histo.fill_color is not None:
                   ax.fill_between(new_x, new_y, [0] * len(new_y), lw=histo.lw, label=histo.legend, color=histo.color, facecolor = histo.fill_color)
               else:
                   ax.plot(new_x, new_y, histo.fmt, lw=histo.lw, label=histo.legend, color=histo.color)
               #ax.bar(map(lambda x: x - 0.5  * x_binwidth, histo.x), histo.y, width=x_binwidth, color=histo.color, label=histo.legend, ec=histo.color)
        else:
            if histo.yerrors is not None:
                lw = histo.lw
                if histo.draw_line is False: lw = 0
                ax.errorbar(histo.x, histo.y, histo.yerrors, elinewidth = histo.lw, lw=lw, label=histo.legend, color=histo.color, marker=histo.marker)
            else:
                ax.plot(histo.x, histo.y, histo.fmt, lw=histo.lw, label=histo.legend, color=histo.color, marker=histo.marker)
    if draw_legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = handles[:]
        labels = labels[:]
        for h, l in extra_legend_items:
            labels.append(l)
            if type(h)==str:
                h = matplotlib.patches.Rectangle((0, 0), 1, 1, fc=h)
            handles.append(h)
        ax.legend(handles, labels, prop = fp, **legend_args)
    if ax.get_legend() is not None:
        map(lambda line: line.set_lw(1.5), ax.get_legend().get_lines())

    if ymin!=None:
        ax.set_ylim(ymin=ymin)
    if ymax!=None:
        ax.set_ylim(ymax=ymax)
    if xmin!=None:
        ax.set_xlim(xmin=xmin)
    if xmax!=None:
        ax.set_xlim(xmax=xmax)
    
    if ax_modifier!=None: ax_modifier(ax)
    if outname is not None: fig.savefig(outname)
    del fig
    
def make_stack(pdatas):
    for i in range(len(pdatas)):
        for j in range(i+1, len(pdatas)):
            pdatas[i].y = map(lambda x: x[0] + x[1], zip(pdatas[i].y, pdatas[j].y))

