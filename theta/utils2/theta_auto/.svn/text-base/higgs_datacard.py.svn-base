# -*- coding: utf-8 -*-
import math
import utils
from Model import *
from root import *

_debug = False


inf = float("inf")
# default nuisance parameter ranges for shape / rate uncertainties:
default_shape_range = [-inf, inf] # combine uses: [-4.0, 4.0]
default_rate_range = [-inf, inf]  # combine uses: [-7.0, 7.0]

def is_int(s):
    try:
        int(s)
        return True
    except ValueError: return False

# line is either a string or a tuple (string, int)
def get_cmds(line):
    if type(line) == tuple:
        line = line[0]
    cmds = [c.strip() for c in line.split()]
    return [c for c in cmds if c!='']


# in theta, names must begin with a letter and consist only of A-Za-z0-9_- and not contain '__'
def transform_name_to_theta(name):
    result = ''
    for c in name:
        if c >= 'a' and c <= 'z' or c >= 'A' and c <= 'Z' or c >='0' and c <='9' or c in ('-', '_'): result += c
        else: result += '_'
    result = result.replace('__', '_')
    if result[0] >= '0' and result[0] <= '9' or result[0]=='-': result = 'tn_' + result
    return result
    
# add entries to the dictionary d. For example,
# d = {}
# add_entry(d, 'a', 'b', 'c', 1.0)
# will make
# d['a']['b']['c'] == 1.0
# valid and True
def add_entry(d, *l):
    if len(l) == 2:
       d[l[0]] = l[1]
       return
    if l[0] not in d: d[l[0]] = {}
    add_entry(d[l[0]], *l[1:])

# general remark on shape uncertainties:
#
# From https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit#Datacard_for_Shape_analyses:
#
# The block of lines defining the mapping (first block in the datacard) contains one or more rows in the form
#
#   * shapes process channel file histogram [ histogram_with_systematics ] 
#
# In this line
#
#   * process is any one the process names, or * for all processes, or data_obs for the observed data
#   * channel is any one the process names, or * for all channels
#   * file, histogram and histogram_with_systematics identify the names of the files and of the histograms within the file, after doing some replacements (if any are found):
#          o $PROCESS is replaced with the process name (or "data_obs" for the observed data)
#          o $CHANNEL is replaced with the channel name
#          o $SYSTEMATIC is replaced with the name of the systematic + (Up, Down)
#          o $MASS is replaced with the higgs mass value which is passed as option in the command line used to run the limit tool 
#
#  [end citation]
#
# In theta, the same implementation is use with some variations:
# * $MASS in NOT implemented.
# * We allow to use some variants for $DIRECTION:
#    - $DIRECTION_updown
#    - $DIRECTION_plusminus
#   If one of these if given, 'Up' or 'Down' are not appended while replacing '$SYSTEMATIC' as decribed above. Instead,
#   '$DIRECTION_...' is replaced with the according direction strip ('up'/'down' or 'plus'/'minus', resp.)
# * So far, each uncertainty specification line MUST be 'shape'. Some higgs datacards use 'shapeN2' which I don't know what it does ...
#   'shape' will interpolate between the nominal and shifted templates using cubiclinear_histomorph using the templates as given
#   in the input root file. After the interpolation, the interpolated histogram is re-normalized to the integral of the nominal one
#   (which has to be identical to the 'rate' specification in the datacard). To treat the rate uncertainty, a seperate log-normal
#   uncertainty is added as part of the coefficient-function. This effectively separates the rate and shape part of such uncertainties.
# * For each channel/process pair to find a shape for, theta will go through the list of 'shapes' lines and pick the first matching one.
#   Make sure that the ordering of these 'shapes' lines is Ok (i.e., put more specific ones first). This will be changed in a future version
#   for better compatibility with "combine" (TODO).



# Exception class used to indicate that something was not found (e.g. file or histogram wityin the file)
class NotFoundException(RuntimeError): pass

# Exception class to indicate that some real inconsistency has been found; the program should probably be aborted:
class InconsistentDataException(RuntimeError): pass



# This method replaces the one-bin data-histogram of the model predictions for the given
# observable and process by a morphing histogram function using the root histogram from filename.
# uncs is a dictionary from uncertainties to factors.
# Models the rate uncertainty via a lognormal term in the coefficient function.
#
# In case len(uncs) is zero, hname_with_systematics can be an empty string (but should be a str, not None)
#
# uncs is a dictionary (uncertainty) -> factor
#
# names for obs, proc, uncs are according to the datacard, not the theta-converted names!
#
# searchpath is the list of paths to look for the root file references form the datacard; it is processed in the given order and the first file found
# is used (usually, it contains '.' and the datacard path). Note that the filename must be unique across all search paths!
#
# rhandling governs how the "rate" part of the shape uncertainty is handled:
#   * "renormalize-lognormal" will perform the template morphing between the three histograms nominal, plus and minus after they have all been scaled to the "nominal" histogram yield.
#     After the template morphing, the resulting interpolated histogram is scaled to the nominal yield. The rate part is then handeled as
#     a log-normal factor in the same parameter as used for the morphing.
#   * "morph" will perform the template morphing between the nominal, plus and minus histograms without any (additional) rescaling.
add_shapes_rootfiles = {}
def add_shapes(model, obs, proc, uncs, filename, hname, hname_with_systematics, include_uncertainties, searchpaths = ['.'], variables = {}, rhandling = 'renormalize-lognormal'):
    assert rhandling in ('renormalize-lognormal', 'morph')
    if filename not in add_shapes_rootfiles:
        path = None
        for s in searchpaths:
            if os.path.isfile(os.path.join(s, filename)):
                path = s
                break
        if path is None: raise NotFoundException, "did not find file '%s' in the paths %s" % (filename, str(searchpaths))
        add_shapes_rootfiles[filename] = rootfile(os.path.join(path, filename))
    rf = add_shapes_rootfiles[filename]
    theta_obs = transform_name_to_theta(obs)
    theta_proc = transform_name_to_theta(proc)
    hname = hname.replace('$CHANNEL', obs)
    hname_with_systematics = hname_with_systematics.replace('$CHANNEL', obs)
    for varname, value in variables.iteritems():
        hname = hname.replace('$%s' % varname, value)
        hname_with_systematics = hname_with_systematics.replace('$%s' % varname, value)
    if proc == 'DATA':
        hname_tmp = hname.replace('$PROCESS', 'DATA')
        histo = rf.get_histogram(hname_tmp, include_uncertainties = False)
        if histo is None:
            hname_tmp = hname.replace('$PROCESS', 'data_obs')
            histo = rf.get_histogram(hname_tmp, include_uncertainties = False)
        if histo is None:
            if _debug: print "note: did not find data histogram in %s" % rf.get_filename()
            raise NotFoundException, "did not find histo"
        model.set_data_histogram(theta_obs, histo, reset_binning = True)
        return
    hf = model.get_histogram_function(theta_obs, theta_proc)
    assert hf is not None, "model has no process '%s' in channel '%s'" % (theta_proc, theta_obs)
    assert len(hf.get_parameters())==0, "model has non-trivial shape uncertainty already"
    old_nominal_histogram = hf.get_nominal_histo()
    assert len(old_nominal_histogram[2])==1, "expected a counting-only histogram with only one bin"
    hname = hname.replace('$PROCESS', proc)
    hname_with_systematics = hname_with_systematics.replace('$PROCESS', proc)
    nominal_histogram = rf.get_histogram(hname, include_uncertainties = include_uncertainties)
    if nominal_histogram is None:
        if _debug: print "note: did not find histogram %s in %s" % (hname, rf.get_filename())
        raise NotFoundException, "did not find histo"
    if _debug:
        print "norm(%s) = %.3f" % (hname, nominal_histogram.get_value_sum())
    # check that histogram in rootfile matches definition in datacard (allow deviations up to 1% / 1e-4 absolute):
    nominal_is_zero = False
    if old_nominal_histogram.get_value_sum() > 0.0 or nominal_histogram.get_value_sum() > 0.0:
        if old_nominal_histogram.get_value_sum() != -1.0 and utils.reldiff(old_nominal_histogram.get_value_sum(), nominal_histogram.get_value_sum()) > 0.01 and abs(old_nominal_histogram.get_value_sum() - nominal_histogram.get_value_sum()) > 1e-4:
            raise InconsistentDataException("add_shapes: histogram normalisation given in datacard and from root file differ by more than 1%% "
                         "(and absolute difference is > 1e-4) for channel %s, process %s (histogram name '%s')" % (obs, proc, hname))
    else:
        print "WARNING: channel '%s' process '%s': yield is <=0. Process will ALWAYS have 0 contribution; please delete it from the datacard." % (obs, proc)
        nominal_is_zero = True
    # even for nominal_is_zero, make sure to set the histogram to ensure that the binning is correct:
    hf.set_nominal_histo(nominal_histogram, reset_binning = True)
    model.reset_binning(theta_obs, nominal_histogram[0], nominal_histogram[1], len(nominal_histogram[2]))
    if len(uncs) == 0: return
    if nominal_is_zero: return
    for u in uncs:
        theta_unc = transform_name_to_theta(u)
        if '$DIRECTION_' in hname_with_systematics:
            hname_plus = hname_with_systematics.replace('$SYSTEMATIC', u)
            hname_minus = hname_plus
            hname_plus = hname_plus.replace('$DIRECTION_plusminus', 'plus')
            hname_minus = hname_plus.replace('$DIRECTION_plusminus', 'minus')
            hname_plus = hname_plus.replace('$DIRECTION_updown', 'up')
            hname_minus = hname_plus.replace('$DIRECTION_updown', 'down')
        else:
            hname_plus = hname_with_systematics.replace('$SYSTEMATIC', u + 'Up')
            hname_minus = hname_with_systematics.replace('$SYSTEMATIC', u + 'Down')
        histo_plus = rf.get_histogram(hname_plus, include_uncertainties = include_uncertainties)
        if histo_plus is None:
            if _debug: print "note: did not find histogram %s in %s" % (hname_plus, rf.get_filename())
            raise NotFoundException, "did not find histo"
        histo_minus = rf.get_histogram(hname_minus, include_uncertainties = include_uncertainties)
        if histo_minus is None:
            if _debug: print "note: did not find histogram %s in %s" % (hname_minus, rf.get_filename())
            raise NotFoundException, "did not find histo"
        if _debug:
            print "norm(%s) = %.3f" % (hname_plus, histo_plus.get_value_sum())
            print "norm(%s) = %.3f" % (hname_minus, histo_minus.get_value_sum())
            
        if rhandling == 'renormalize-lognormal':
            # make the rate uncertainty part of the coefficient function, i.e., normalize plus and minus histograms
            # to nominal and add a lognormal uncertainty to the coefficient function:
            lambda_plus = math.log(histo_plus.get_value_sum() / nominal_histogram.get_value_sum()) * uncs[u]
            lambda_minus = -math.log(histo_minus.get_value_sum() / nominal_histogram.get_value_sum()) * uncs[u]
            model.get_coeff(theta_obs, theta_proc).add_factor('exp', parameter = u, lambda_plus = lambda_plus, lambda_minus = lambda_minus)
            f_plus = nominal_histogram.get_value_sum() / histo_plus.get_value_sum()
            histo_plus = histo_plus.scale(f_plus)
            f_minus = nominal_histogram.get_value_sum() / histo_minus.get_value_sum()
            histo_minus = histo_minus.scale(f_minus)
            hf.set_syst_histos(u, histo_plus, histo_minus, uncs[u])
            hf.normalize_to_nominal = True
        else:
            hf.set_syst_histos(u, histo_plus, histo_minus, uncs[u])

            
class ParametricShapeBuilder:
    def __init__(self, model):
        self.model = model
        self.obs, self.proc =  None, None
    
    # c = None mean c = lambda0
    # free = True means to give the prior an infinite width
    def exponential(self, parameter, lambda0, lambda0_delta, binning = None, binborders = None):
        free = False
        if lambda0_delta == inf:
            free = True
            lambda0_delta = 1.0
        if parameter in self.model.distribution.get_parameters():
            # double-check that the parameter prior is set to the value we need:
            pars = self.model.get_distribution_parameters(parameter)
            assert pars['mean'] != 0.0, "incompatible distribution for parameter %s: exponential needs mean=0.0" % paramater
            if free:
                assert pars['width'] == inf, "incompatible distribution for parameter %s: exponential with infinite lambda0_delta requires width=inf" % paramater
            else:
                assert pars['width'] == 1.0, "incompatible distribution for parameter %s: exponential with finite lambda0_delta requires width=1.0" % paramater
        else:
            self.model.distribution.set_distribution(parameter, 'gauss', mean = 0.0, width = inf if free else 1.0, range = [-inf, inf])
        hf_old = self.model.get_histogram_function(self.obs, self.proc)
        norm = hf_old.get_nominal_histo().get_value_sum()
        if _debug: print "exponential: found old norm = ", norm
        # find out binning:
        if binning is not None and binborders is not None:
            raise RuntimeError, "both binning and binborders specifies, please specify only one of those!"
        if binning is not None:
            xmin, xmax, nbins = binning
            dx = (xmax - xmin) / nbins
            binborders = [xmin + dx * i for i in range(nbins+1)]
        elif binborders is None:
            xmin, xmax, nbins = self.model.observables[self.obs]
            print "NOTE: pshape exponential: assuming equidistant binning with (xmin, xmax, nbins) = (%.2f, %.2f, %d)" % (xmin, xmax, nbins)
            dx = (xmax - xmin) / nbins
            binborders = [xmin + dx * i for i in range(nbins+1)]
        # check that binborders make sense:
        for i in range(1, len(binborders)):
            assert binborders[i-1] < binborders[i], "binborders must increase monotonically!"
        hf = ExponentialHistogramFunction(lambda0, lambda0_delta, parameter, norm, binborders)
        self.model.set_histogram_function(self.obs, self.proc, hf)
        
        
    def apply_(self, theta_obs, theta_proc, command):
        if _debug: print "PSB: apply_ %s, %s, %s" % (theta_obs, theta_proc, command)
        self.obs = theta_obs
        self.proc = theta_proc
        exec('self.' + command)

def build_model(fname, filter_cp = lambda chan, p: True, filter_uncertainty = lambda unc: True, include_mc_uncertainties = False, variables = {}, rmorph_method = 'renormalize-lognormal'):
    """
    Build a Model from a text-based datacard as used in LHC Higgs analyses

    See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit

    Note that not the complete set of features is supported, in particular no unbinned fits.
    Supported uncertainties are: lnN (symmetric and asymmetric), lnU, gmN, shape

    The 'shape' uncertainty uses a slightly different interpolation: the Higgs tool uses a quadratic interpolation with linear extrapolation
    whereas theta uses a cubic interpolation and linear extrapolation. It is expected that this has negligible impact
    on the final result, but it might play a role in extreme cases.
    
    Parameters:
    
    * ``fname`` is the filename of the datacard to process
    * ``filter_cp`` is a function which, for a given pair of a channel name and process name (as given in the model configuration file), returns ``True`` if this channel/process should be kept and ``False`` otherwise. The default is to keep all channel/process combinations.
    * ``filter_uncertainty`` is a filter function for the uncertainties. The default is to keep all uncertainties
    * ``include_mc_uncertainties`` if ``True`` use the histogram uncertainties of shapes given in root files for Barlow-Beeston light treatment of MC stat. uncertainties
    * ``variables`` is a dictionary for replacing strings in the datacards. For example, use ``variables = {'MASS': '125'}`` to replace each appearance of '$MASS' in the datacard with '125'. Both key and value should be strings.
    * ``rmorph_method`` controls how the rate part of a shape uncertainty is handled: "renormalize-lognormal" will re-scale the plus and minus histogram to the nominal one,
      perform the morphing on those histograms, re-scale the morphed histogram to the nominal one and add an exponential (=log-normal) rate factor using the same parameter as is used
      for the interpolation. Instead "morph" will simply interpolate between the nominal, plus and minus histograms as they are.
    """
    model = Model()
    lines = [l.strip() for l in file(fname)]
    lines = [(lines[i], i+1) for i in range(len(lines)) if not lines[i].startswith('#') and lines[i]!='' and not lines[i].startswith('--')]
    
    cmds = get_cmds(lines[0])
    while cmds[0] != 'imax':
        print 'WARNING: ignoring line %d ("%s") at beginning of file as first token is "%s", not "imax", although not marked as comment' % (lines[0][1], lines[0][0], cmds[0])
        lines = lines[1:]
        cmds = get_cmds(lines[0])
    assert cmds[0]=='imax', "Line %d: Expected imax statement as first statement in the file" % lines[0][1]
    imax = cmds[1]
    if imax !='*': imax = int(imax)
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='jmax', "Line %d: Expected 'jmax' statement directly after 'imax' statement" % lines[0][1]
    #jmax = int(cmds[1])
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='kmax', "Line %d: Expected 'kmax' statement directly after 'jmax' statement" % lines[0][1]
    if cmds[1] == '*': kmax = -1
    else: kmax = int(cmds[1])
    lines = lines[1:]

    shape_lines = []
    shape_observables = set()
    cmds = get_cmds(lines[0])
    while cmds[0].lower() == 'shapes':
        assert len(cmds) in (5,6)
        if len(cmds) == 5: cmds.append('')
        shape_lines.append(cmds[1:]) # (process, channel, file, histogram, histogram_with_systematics)
        obs = cmds[2]
        shape_observables.add(obs)
        lines =lines[1:]
        cmds = get_cmds(lines[0])
        
    pshape_lines = []
    while cmds[0].lower() == 'pshape':
        assert len(cmds) >= 4 # pshape channel proc command
        pshape_lines.append([cmds[1], cmds[2], ' '.join(cmds[3:])]) # (process, channel, command)
        if _debug: print "Line %d: found pshape line %s" % (lines[0][1], str(pshape_lines[-1]))
        lines =lines[1:]
        cmds = get_cmds(lines[0])

    assert cmds[0].lower() in ('bin', 'observation'), "Line %d: Expected 'bin' or 'observation' statement" % lines[0][1]
    if cmds[0].lower() == 'bin':
        # prepend a 'c' so we can use numbers as channel names:
        channel_labels = cmds[1:]
        if imax=='*': imax = len(channel_labels)
        assert len(channel_labels) == imax, "Line %d: Number of processes from 'imax' and number of labels given in 'bin' line (%s) mismatch" % (lines[0][1], str(channel_labels))
        lines = lines[1:]
        cmds = get_cmds(lines[0])
    else:
        channel_labels = [ '%d' % i for i in range(1, imax + 1)]
    assert cmds[0].lower()=='observation', "Line %d: Expected 'observation' statement directly after fist 'bin' statement" % lines[0][1]
    observed_flt = [float(o) for o in cmds[1:]]
    observed_int = map(lambda f: int(f), observed_flt)
    if observed_flt != observed_int: raise RuntimeError, "Line %d: non-integer events given in 'observed' statement!" % lines[0][1]
    if imax=='*': imax = len(observed_int)
    assert len(observed_int) == imax, "Line %d: Number of processes from 'imax' and number of bins given in 'observed' mismatch: imax=%d, given in observed: %d" % (lines[0][1], imax, len(observed))
    for i in range(len(channel_labels)):
        theta_obs = transform_name_to_theta(channel_labels[i])
        model.set_data_histogram(theta_obs, Histogram(0.0, 1.0, [observed_flt[i]]))
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0] == 'bin', "Line %d: Expected 'bin' statement"% lines[0][1]
    # save the channel 'headers', to be used for parsing the next line:
    channels_for_table = cmds[1:]
    for c in channels_for_table:
        if c not in channel_labels: raise RuntimeError, "Line % d: unknown channel '%s'" % (lines[0][1], c)
    lines = lines[1:]
    n_cols = len(channels_for_table)

    cmds = get_cmds(lines[0])
    assert cmds[0]=='process'
    processes1 = cmds[1:]
    if len(processes1) != n_cols:
        raise RuntimeError, "Line %d: 'bin' statement and 'process' statement have different number of elements" % lines[0][1]
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='process', "Line %d: Expected second 'process' line directly after first" % lines[0][1]
    processes2 = cmds[1:]
    if n_cols != len(processes2):
        raise RuntimeError, "Line %d: 'process' statements have different number of elements" % lines[0][1]
    lines = lines[1:]
    
    # get process names and numeric process ids:
    if(all(map(is_int, processes1))):
        process_ids_for_table = [int(s) for s in processes1]
        processes_for_table = processes2
    else:
        if not all(map(is_int, processes2)): raise RuntimeError("just before line %d: one of these 'process' lines should contain only numbers!" % lines[0][1])
        process_ids_for_table = [int(s) for s in processes2]
        processes_for_table = processes1
    
    # build a list of columns to keep ( = not filtered by filter_cp):
    column_indices = []
    for i in range(n_cols):
        if filter_cp(channels_for_table[i], processes_for_table[i]): column_indices.append(i)

    # check process label / id consistency:
    p_l2i = {}
    p_i2l = {}
    for i in range(n_cols):
        p_l2i[processes_for_table[i]] = process_ids_for_table[i]
        p_i2l[process_ids_for_table[i]] = processes_for_table[i]
    # go through again to make check, also save signal processes:
    signal_processes = set()
    for i in range(n_cols):
        if p_l2i[processes_for_table[i]] != process_ids_for_table[i] or p_i2l[process_ids_for_table[i]] != processes_for_table[i]:
            raise RuntimeError, "Line %d: mapping process id <-> process label (defined via the two 'process' lines) is not one-to-one as expected!" % lines[0][1]
        if p_l2i[processes_for_table[i]] <= 0:
            signal_processes.add(processes_for_table[i])

    cmds = get_cmds(lines[0])
    assert cmds[0]=='rate', "Line %d: Expected 'rate' statement after the two 'process' statements" % lines[0][1]
    if n_cols != len(cmds)-1:
        raise RuntimeError, "Line %d: 'rate' statement does specify the wrong number of elements" % lines[0][1]
    for i in column_indices:
        theta_obs, theta_proc = transform_name_to_theta(channels_for_table[i]), transform_name_to_theta(processes_for_table[i])
        n_exp = float(cmds[i+1])
        #print o,p,n_exp
        hf = HistogramFunction()
        hf.set_nominal_histo(Histogram(0.0, 1.0, [n_exp]))
        #print "setting prediction for (theta) channel '%s', (theta) process '%s'" % (theta_obs, theta_proc)
        model.set_histogram_function(theta_obs, theta_proc, hf)
        assert model.get_histogram_function(theta_obs, theta_proc) is not None
    lines = lines[1:]
    
    kmax  = len(lines)

    if kmax != len(lines):
        raise RuntimeError, "Line %d--end: wrong number of lines for systematics (expected kmax=%d, got %d)" % (lines[0][1], kmax, len(lines))
    
    # save uncertainty names to avoid duplicates:
    uncertainty_names = set()
    # shape systematics is a dictionary (channel) --> (process) --> (parameter) --> (factor)
    # factors of 0 are omitted.
    shape_systematics = {}
    for i in range(kmax):
        if _debug: print "processing line %d" % lines[i][1]
        cmds = get_cmds(lines[i])
        assert len(cmds) >= len(processes_for_table) + 2, "Line %d: wrong number of entries for uncertainty '%s'" % (lines[i][1], cmds[0])
        if not filter_uncertainty(cmds[0]): continue
        if cmds[0] in uncertainty_names:
            raise RuntimeError, "Uncertainty '%s' specified more than once; this is not supported." % cmds[0]
        uncertainty_names.add(cmds[0])
        uncertainty = transform_name_to_theta(cmds[0])
        if cmds[1] == 'gmN':
            values = cmds[3:]
            n_affected = 0
            k = float(cmds[2])
            for icol in column_indices:
                if values[icol]=='-': continue
                val = float(values[icol])
                if val==0.0: continue
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                # add the same parameter (+the factor in the table) as coefficient:
                model.get_coeff(obsname, procname).add_factor('id', parameter = uncertainty)
                n_affected += 1
                n_exp = model.get_histogram_function(obsname, procname).get_nominal_histo()[2][0]
                if max(n_exp, val*k) != 0:
                     if abs(n_exp - val*k)/max(n_exp, val*k) > 0.03:
                            raise RuntimeError, "gmN uncertainty %s for process %s is inconsistent: the rate expectation should match k*theta but N_exp=%f, k*theta=%f!" % (cmds[0], procname, n_exp, val*k)
            if n_affected > 0:
                n_obs_sb = float(cmds[2]) # the number of observed events in the sideband
                n_model_sb = n_obs_sb     # the number of events in the model template in the sideband. This is pretty arbitrary, as a scale factor is used to fit this anyway.
                if n_model_sb == 0.0: n_model_sb = 1.0
                hf = HistogramFunction()
                hf.set_nominal_histo(Histogram(0.0, 1.0, [n_model_sb]))
                obs_sb = '%s_sideband' % uncertainty
                model.set_histogram_function(obs_sb, 'proc_%s_sideband' % uncertainty, hf)
                model.set_data_histogram(obs_sb, Histogram(0.0, 1.0, [n_obs_sb]))
                model.get_coeff(obs_sb, 'proc_%s_sideband' % uncertainty).add_factor('id', parameter = uncertainty)
                # as mean, use the value at the observation such that toys reproduce this value ...
                model.distribution.set_distribution(uncertainty, 'gauss', mean = n_obs_sb / n_model_sb, width = float("inf"), range = (0.0, float("inf")))
        elif cmds[1] in ('lnN', 'lnU'):
            n_affected = 0
            values = cmds[2:]
            for icol in column_indices:
                if values[icol]=='-': continue
                if '/' in values[icol]:
                    p = values[icol].find('/')
                    lambda_minus = -math.log(float(values[icol][0:p]))
                    lambda_plus = math.log(float(values[icol][p+1:]))
                else:
                    lambda_minus = math.log(float(values[icol]))
                    lambda_plus = lambda_minus
                if lambda_plus == 0.0 and lambda_minus == 0.0: continue
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                n_affected += 1
                model.get_coeff(obsname, procname).add_factor('exp', parameter = uncertainty, lambda_minus = lambda_minus, lambda_plus = lambda_plus)
            if n_affected > 0:
                if cmds[1] == 'lnN':  model.distribution.set_distribution(uncertainty, 'gauss', mean = 0.0, width = 1.0, range = default_rate_range)
                else:  model.distribution.set_distribution(uncertainty, 'gauss', mean = 0.0, width = inf, range = [-1.0, 1.0])
        elif cmds[1] == 'gmM':
            values = cmds[2:]
            values_f = set([float(s) for s in values if float(s)!=0.0])
            if len(values_f)>1: raise RunetimeError, "gmM does not support different uncertainties"
            if len(values_f)==0: continue
            n_affected = 0
            for icol in column_indices:
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                model.get_coeff(obsname, procname).add_factor('id', parameter = uncertainty)
                n_affected += 1
            if n_affected > 0:
                model.distribution.set_distribution(uncertainty, 'gamma', mean = 1.0, width = float(values[icol]), range = (0.0, float("inf")))
        elif cmds[1] in 'shape':
            factors = cmds[2:]
            n_affected = 0
            for icol in column_indices:
                if factors[icol] == '-' or float(factors[icol]) == 0.0: continue
                factor = float(factors[icol])
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                n_affected += 1
                add_entry(shape_systematics, channels_for_table[icol], processes_for_table[icol], cmds[0], factor)
            if n_affected > 0:
                model.distribution.set_distribution(uncertainty, 'gauss', mean = 0.0, width = 1.0, range = default_shape_range)
        else: raise RuntimeError, "Line %d: unknown uncertainty type %s" % (lines[0][1], cmds[1])
    # add shape systematics:
    if '*' in shape_observables: shape_observables = set(channel_labels)
    data_done = set()
    searchpaths = ['.', os.path.dirname(fname)]
    if _debug: print "adding shapes now ..."
    psb = ParametricShapeBuilder(model)
    # loop over processes and observables:
    for icol in column_indices:
        obs = channels_for_table[icol]
        if _debug: print "adding shape for channel '%s'" % obs
        if obs not in shape_observables: continue
        proc = processes_for_table[icol]
        found_matching_shapeline = False
        # try all lines in turn, until adding the shapes from that file succeeds:
        for l in shape_lines: # l = (process, channel, file, histogram, histogram_with_systematics)
            try:
                if _debug: print "   shape line: %s" % str(l)
                if l[1]!='*' and l[1]!=obs: continue
                if obs not in data_done and l[0] in ('*', 'data_obs', 'DATA'):
                    try:
                        add_shapes(model, obs, 'DATA', {}, l[2], l[3], '', include_mc_uncertainties, searchpaths = searchpaths, variables = variables, rhandling = rmorph_method)
                        data_done.add(obs)
                    except RuntimeError: pass # ignore missing data
                if l[0]!='*' and l[0]!=proc: continue
                uncs = {}
                if obs in shape_systematics: uncs = shape_systematics[obs].get(proc, {})
                if _debug: print "   adding shapes for channel %s, process %s, trying file %s, line %s" % (obs, proc, l[2], ' '.join(l))
                add_shapes(model, obs, proc, uncs, l[2], l[3], l[4], include_mc_uncertainties, searchpaths = searchpaths, variables = variables, rhandling = rmorph_method)
                found_matching_shapeline = True
                break
            except NotFoundException: pass # ignore the case that some histo has not been found for now; raise a RuntimeError later if no line matched
        if found_matching_shapeline: continue
        # now that nothing is found, also try the pshapelines
        for l in pshape_lines:
            if l[0] != proc or l[1] != obs: continue
            print "Trying to apply pshape line %s" % str(l)
            theta_obs, theta_proc = transform_name_to_theta(obs), transform_name_to_theta(proc)
            psb.apply_(theta_obs, theta_proc, l[2])
            found_matching_shapeline = True
        if not found_matching_shapeline:
            raise RuntimeError, "Did not find all the (nominal / systematics) histograms for channel '%s', process '%s'" % (obs, proc)
    model.set_signal_processes([transform_name_to_theta(proc) for proc in signal_processes])
    if include_mc_uncertainties: model.bb_uncertainties = True
    return model

