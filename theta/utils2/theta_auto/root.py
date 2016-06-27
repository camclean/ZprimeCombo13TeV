# -*- coding: utf-8 -*-
import ROOT
ROOT.gROOT.SetBatch(True)
import re, fnmatch, math, copy
import os, os.path
import array
import utils

from theta_auto.Model import *


class rootfile:
    # these are caching dictionaries indexing by filename
    tfiles = {}
    
    def __init__(self, filename):
        assert os.path.isfile(filename), "File %s not found (cwd: %s)" % (filename, os.getcwd())
        self.filename = filename
        if self.filename not in rootfile.tfiles: rootfile.tfiles[self.filename] = ROOT.TFile(filename, "read")
        self.tfile = rootfile.tfiles[self.filename]
        
    @staticmethod
    def th1_to_histo(th1, include_uncertainties, include_uoflow = False):
        xmin, xmax, nbins = th1.GetXaxis().GetXmin(), th1.GetXaxis().GetXmax(), th1.GetNbinsX()
        if include_uoflow:
            binwidth = (xmax - xmin) / nbins
            xmin -= binwidth
            xmax += binwidth
            nbins += 2
            values = array.array('d', [th1.GetBinContent(i) for i in range(nbins)])
            uncertainties = array.array('d', [th1.GetBinError(i) for i in range(nbins)]) if include_uncertainties else None
            x_low  = [xmin] + [th1.GetBinLowEdge(i) for i in range(1, nbins)]
        else:
            values = array.array('d', [th1.GetBinContent(i) for i in range(1, nbins+1)])
            uncertainties = array.array('d', [th1.GetBinError(i) for i in range(1, nbins+1)]) if include_uncertainties else None
            x_low  = [th1.GetBinLowEdge(i) for i in range(1, nbins+1)]
        h = Histogram(xmin, xmax, values, uncertainties, th1.GetName(), x_low = x_low)
        return h
  
    # get all templates as dictionary (histogram name) --> Histogram instance
    # only checks type of histogram, not naming convention
    def get_all_templates(self, include_uncertainties, warn = True, include_uoflow = False):
        result = {}
        l = self.tfile.GetListOfKeys()
        for key in l:
            clas = key.GetClassName()
            if clas == 'TDirectoryFile': continue
            if clas not in ('TH1F', 'TH1D') and warn:
                print "WARNING: ignoring key %s in input file %s because it is of ROOT class %s, not TH1F / TH1D" % (key.GetName(), self.filename, clas)
                continue
            th1 = key.ReadObj()
            result[str(key.GetName())] = rootfile.th1_to_histo(th1, include_uncertainties, include_uoflow = include_uoflow)
        return result
        
    def get_filename(self): return self.filename
        
    def get_histogram(self, hname, include_uncertainties, fail_with_exception = False):
        h = self.tfile.Get(hname)
        if not h.Class().InheritsFrom("TH1"):
            if fail_with_exception: raise RuntimeError, "histogram '%s' in root file '%s' not found!" % (hname, self.tfile.GetName())
            else: return None
        return rootfile.th1_to_histo(h, include_uncertainties)



# flatten a nested dictionary with string indices into a flat dictionary. The new key names
# are given by <key1>__<key2>...
#
# note that the nesting is only flattened if the key is a str and the value is a dict.
def _flatten_nested_dict(d, result, current_key = ''):
    for k in d:
        if type(k)==str and type(d[k]) == dict:
            if current_key == '':   new_current_key = k
            else: new_current_key = current_key + '__' + k
            _flatten_nested_dict(d[k], result, new_current_key)
        else:
            if current_key == '':   new_current_key = k
            else: new_current_key = current_key + '__' + k
            result[new_current_key] = d[k]


# in the dictionary d, lists are converted to dictionaries
def _lists_to_dict(d):
    for k in d:
        if type(d[k])==list:
            l = d[k]
            d[k] = dict([('%d' % i, l[i]) for i in range(len(l))])
        elif type(d[k])==dict:
            _lists_to_dict(d[k])


def write_histograms_to_rootfile(histograms_, rootfilename):
    """
    Parameters:
    
     * ``histograms`` - a dictionary with strings as key name and :class:`Histogram` s as value; nested dictionaries are allowed
     * ``rootfilename`` - the filename of the root file to create. Will be overwritten if it already exists.
    
    Note that the name of the TH1Ds in the output root file is constructed via the key names in the dictionary: the
    key name is given by concatenating all key names required to retrive the histogram in the ``histograms`` parameter,
    separated by "__". For example if histograms['channel1']['proc1'] is a Histogram, its name in the root output file
    will be ``channel1__proc1``.
    """
    outfile = ROOT.TFile(rootfilename, "recreate")
    histograms = copy.deepcopy(histograms_)
    flattened = {}
    _lists_to_dict(histograms)
    _flatten_nested_dict(histograms, flattened)
    for name in flattened:
        h = flattened[name]
        if not isinstance(h, Histogram): raise RuntimeError, "a non-Histogram value encountered in histograms parameter (name: %s; value: %s); this is not allowed." % (name, str(h))
        root_h = ROOT.TH1D(name, name, h.get_nbins(), h.get_xmin(), h.get_xmax())
        h_values = h.get_values()
        for i in range(h.get_nbins()):
            root_h.SetBinContent(i+1, h_values[i])
        root_h.Write()
    outfile.Close()
    


def build_model_from_rootfile(filenames, histogram_filter = lambda s: True, root_hname_to_convention = lambda s: s, transform_histo = lambda h: h, include_mc_uncertainties = False):
    """
    Build a multi-channel model based on template morphing from histograms in a root file
    
    This root file is expected to contain all the templates of the model adhering to a certain naming scheme:
      ``<observable>__<process>``     for the "nominal" templates (=not affect by any uncertainty) and
      ``<observable>__<process>__<uncertainty>__(plus,minus)``  for the "shifted" templates to be used for template morphing.

    
    ``<observable>``, ``<process>``, and ``<uncertainty>`` are names you can choose at will as long as it does not contain '__'. You are encouraged
    to choose sensible names as these names are used in the output a lot.

    For example, if you want to make a combined statistical evaluation of a muon+jets and an electron+jets ttbar cross section measurement,
    you can name the observables "mu" and "ele"; the processes might be "ttbar", "w", "nonw", the uncertainties "jes", "q2". Provided
    all uncertainties affect all template shapes, you would supply 6 nominal and 24 "uncertainty" templates:
    The 6 nominal would be: mu__ttbar, mu__w, mu__nonw, ele__ttbar, ele__w, ele__nonw
    Some of the 24 "uncertainty" histograms would be: mu__ttbar__jes__plus, mu__ttbar__jes__minus, ..., ele__nonw__q2__minus
    
    All templates of one observable must have the same range and binning. All templates should be normalized
    to the same luminosity (although normalization can be changed from the analysis python script later, this is generally not recommended, unless
    scaling everything to a different lumi).

    It is possible to omit some of the systematic templates completely. In this case, it is assumed
    that the presence of that uncertainty has no influence on this process in this observable.

    Observed data has the special process name "DATA" (all capitals!), so for each observable, there should be exactly one ``<observable>_DATA``
    histogram, if you have data at all. If you do not have data, just omit this; the methods will be limited to calculating the expected
    result.

    To identify which process should be considered as signal, call :meth:`Model.set_signal_processes` after constructing the Model.
        

    The model built is based on the given templates where the systematic uncertainties are fully correlated across different
    observables and processes, i.e., the same parameter is used to interpolate between the nominal and shifted templates
    if the name of the uncertainty is the same. Two different systematic uncertainties (=with different names) are assumed to be uncorrelated.
    Each parameter has a Gaussian prior with width 1.0 and mean 0.0 and has the same name as the uncertainty. You can use
    the functions in Distribution (e.g., via model.distribution) to override this prior. This is useful if the "plus" and "minus" templates
    are not the +-1sigma deviations, but, say, the +-2sigma in which case you can use a prior with width 0.5.

    Parameters:
    
    * ``filenames`` is either a single string or a list of strings speficiying the root file names to read the histograms from.
    * ``histogram_filter`` is a function which -- given a histogram name as in the root file --
      returns either ``True`` to keep histogram or ``False`` to ignore the histogram. The default is to keep all histograms.
      This is useful if you want to consider only a subset of channels or uncertainties.
    * ``root_hname_to_convention`` is a function which get the "original" histogram name (as in the root file) to histogram names as expected by the
       naming convention as described above. The default is to not modify the names.
    * ``transform_histo`` is a function which takes one parameter, the :class:`Histogram` instance as read from the root file (but the name already transformed
       using ``root_hname_to_convention``). This method should return a :class:`Histogram` instance which should be used. This is useful e.g. for re-binning or scaling
       Histograms "on the fly", without having to re-create the root input file. Note that the name of the returned Histogram must be the same as the input Histogram(!)
    * ``include_mc_uncertainties`` is a boolean which specifies whether or not to use the Histogram uncertainties as Monte-Carlo statistical uncertainties and include
      their treatment in the statistical methods using the "barlow-Beeston light" method (see also :ref:`model_intro`).
      
    """
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
                try:
                   hf.set_syst_histos('%s' % u, histos['%s__%s__%s__plus' % (o, p, u)], histos['%s__%s__%s__minus' % (o, p, u)])
                except Exception, e:
                   print "Exception while setting syst histos for channel %s, process %s" % (o, p)
                   raise e
            result.set_histogram_function(o, p, hf)
    for u in uncertainties:
        result.distribution.set_distribution('%s' % u, 'gauss', mean = 0.0, width = 1.0, range = (-float("inf"), float("inf")))
    result.bb_uncertainties = include_mc_uncertainties
    return result

