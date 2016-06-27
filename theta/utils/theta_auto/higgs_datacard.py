# -*- coding: utf-8 -*-
import math
import utils
from Model import *

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
add_shapes_rootfiles = {}
def add_shapes(model, obs, proc, uncs, filename, hname, hname_with_systematics, include_uncertainties):
    if filename not in add_shapes_rootfiles:
        add_shapes_rootfiles[filename] = rootfile(filename)
    rf = add_shapes_rootfiles[filename]
    theta_obs = transform_name_to_theta(obs)
    theta_proc = transform_name_to_theta(proc)
    hname = hname.replace('$CHANNEL', obs)
    hname_with_systematics = hname_with_systematics.replace('$CHANNEL', obs)
    if proc == 'DATA':
        hname_tmp = hname.replace('$PROCESS', 'DATA')
        histo = rf.get_histogram(hname_tmp, include_uncertainties = False)
        if histo is None:
            hname_tmp = hname.replace('$PROCESS', 'data_obs')
            histo = rf.get_histogram(hname_tmp, include_uncertainties = False)
        if histo is None: raise RuntimeError, "did not find data histogram in rootfile"
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
    if utils.reldiff(sum(old_nominal_histogram[2]), sum(nominal_histogram[2])) > 0.01 and abs(sum(old_nominal_histogram[2]) - sum(nominal_histogram[2])) > 1e-4:
        raise RuntimeError, "add_shapes: histogram normalisation given in datacard and from root file differ by more than >1% (and absolute difference is > 1e-4)"
    hf.set_nominal_histo(nominal_histogram, reset_binning = True)
    model.reset_binning(theta_obs, nominal_histogram[0], nominal_histogram[1], len(nominal_histogram[2]))
    if len(uncs) == 0: return
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
        histo_plus = rf.get_histogram(hname_plus, include_uncertainties = include_uncertainties, fail_with_exception = True)
        histo_minus = rf.get_histogram(hname_minus, include_uncertainties = include_uncertainties, fail_with_exception = True)
        # make the rate uncertainty part of the coefficient function, i.e., normalize plus and minus histograms
        # to nominal and add a lognormal uncertainty to the coefficient function:
        lambda_plus = math.log(sum(histo_plus[2]) / sum(nominal_histogram[2])) * uncs[u]
        lambda_minus = -math.log(sum(histo_minus[2]) / sum(nominal_histogram[2])) * uncs[u]
        model.get_coeff(theta_obs, theta_proc).add_factor('exp', parameter = u, lambda_plus = lambda_plus, lambda_minus = lambda_minus)
        f_plus = sum(nominal_histogram[2]) / sum(histo_plus[2])
        utils.mul_list(histo_plus[2], f_plus)
        f_minus = sum(nominal_histogram[2]) / sum(histo_minus[2])
        utils.mul_list(histo_minus[2], f_minus)
        hf.set_syst_histos(u, histo_plus, histo_minus, uncs[u])
        hf.normalize_to_nominal = True
 
## \brief Build a Model from a datacard as used in LHC Higgs analyses
# 
# See https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit
#
# Note that not the complete set of features is supported, in particular no unbinned fits.
# Supported uncertainties are: lnN (symmetric and asymmetric), gmM, gmN, shape
#
# The 'shape' uncertainty uses a slightly different interpolation: the Higgs tool uses a quadratic interpolation with linear extrapolation
# whereas theta uses a cubic interpolation and linear extrapolation. It is expected that this has negligible impact
# on the final result, but it might play a role in extreme cases (?)
#
# \param fname is the filename of the datacard to process.
#
# \param filter_channel is a function which, for each channel name (as given in the model configuration in fname), returns
# True if this channel should be kept and False otherwise. The default is to keep all channels.
def build_model(fname, filter_channel = lambda chan: True, filter_uncertainty = lambda unc: True, debug = False, include_mc_uncertainties = False):
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
        if not filter_channel(channel_labels[i][1:]): continue
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
    processes_for_table = cmds[1:]
    if len(processes_for_table) != n_cols:
        raise RuntimeError, "Line %d: 'bin' statement and 'process' statement have different number of elements" % lines[0][1]
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='process', "Line %d: Expected second 'process' line directly after first" % lines[0][1]
    process_ids_for_table = [int(s) for s in cmds[1:]]
    if n_cols != len(process_ids_for_table):
        raise RuntimeError, "Line %d: 'process' statements have different number of elements" % lines[0][1]
    lines = lines[1:]

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
    for i in range(n_cols):
        if not filter_channel(channels_for_table[i]): continue
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
    
    # shape systematics is a dictionary (uncertainty) --> (channel) --> (process) --> (factor)
    # factors of 0 are omitted.
    shape_systematics = {}
    for i in range(kmax):
        if debug: print "processing line %d" % lines[i][1]
        cmds = get_cmds(lines[i])
        assert len(cmds) >= len(processes_for_table) + 2, "Line %d: wrong number of entries for uncertainty '%s'" % (lines[i][1], cmds[0])
        if not filter_uncertainty(cmds[0]): continue
        uncertainty = transform_name_to_theta(cmds[0])
        if cmds[1] == 'gmN':
            values = cmds[3:]
            n_affected = 0
            k = float(cmds[2])
            for icol in range(n_cols):
                if values[icol]=='-': continue
                val = float(values[icol])
                if val==0.0: continue
                if not filter_channel(channels_for_table[icol]): continue
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                #print cmds[0], k, obsname, procname, values[icol], val
                # add the same parameter (+the factor in the table) as coefficient:
                model.get_coeff(obsname, procname).add_factor('id', parameter = uncertainty)
                n_affected += 1
                n_exp = model.get_histogram_function(obsname, procname).get_nominal_histo()[2][0]
                #print n_exp
                if abs(n_exp - val*k)/max(n_exp, val*k) > 0.03:
                    raise RuntimeError, "gmN uncertainty %s for process %s is inconsistent: the rate expectation should match k*theta but N_exp=%f, k*theta=%f!" % (cmds[0], procname, n_exp, val*k)
            if n_affected > 0:
                k = float(cmds[2])
                hf = HistogramFunction()
                hf.set_nominal_histo((0.0, 1.0, [k]))
                obs_sb = '%s_sideband' % uncertainty
                model.set_histogram_function(obs_sb, 'proc_sb', hf)
                model.set_data_histogram(obs_sb, (0.0, 1.0, [k]))
                model.get_coeff(obs_sb, 'proc_sb').add_factor('id', parameter = uncertainty)
                # the maximum likelihood estimate for the delta parameter is 1.0
                model.distribution.set_distribution(uncertainty, 'gauss', mean = 1.0, width = float("inf"), range = (0.0, float("inf")))
        elif cmds[1] == 'lnN':
            n_affected = 0
            values = cmds[2:]
            for icol in range(n_cols):
                if values[icol]=='-': continue
                if not filter_channel(channels_for_table[icol]): continue
                if '/' in values[icol]:
                    p = values[icol].find('/')
                    lambda_minus = -math.log(float(values[icol][0:p]))
                    lambda_plus = math.log(float(values[icol][p+1:]))
                else:
                    lambda_minus = math.log(float(values[icol]))
                    lambda_plus = lambda_minus
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                n_affected += 1
                #print cmds[0], obsname, procname, lambda_minus, lambda_plus
                model.get_coeff(obsname, procname).add_factor('exp', parameter = uncertainty, lambda_minus = lambda_minus, lambda_plus = lambda_plus)
            if n_affected > 0:
                model.distribution.set_distribution(uncertainty, 'gauss', mean = 0.0, width = 1.0, range = (-7.0, 7.0))
        elif cmds[1] == 'gmM':
            values = cmds[2:]
            values_f = set([float(s) for s in values if float(s)!=0.0])
            if len(values_f)>1: raise RunetimeError, "gmM does not support different uncertainties"
            if len(values_f)==0: continue
            n_affected = 0
            for icol in range(n_cols):
                if not filter_channel(channels_for_table[icol]): continue
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                model.get_coeff(obsname, procname).add_factor('id', parameter = uncertainty)
                n_affected += 1
            if n_affected > 0:
                model.distribution.set_distribution(uncertainty, 'gamma', mean = 1.0, width = float(values[icol]), range = (0.0, float("inf")))
        elif cmds[1] in 'shape':
            factors = cmds[2:]
            n_affected = 0
            for icol in range(n_cols):
                if not filter_channel(channels_for_table[icol]): continue
                if factors[icol] == '-' or float(factors[icol]) == 0.0: continue
                factor = float(factors[icol])
                obsname = transform_name_to_theta(channels_for_table[icol])
                procname = transform_name_to_theta(processes_for_table[icol])
                n_affected += 1
                add_entry(shape_systematics, channels_for_table[icol], processes_for_table[icol], cmds[0], factor)
            if n_affected > 0:
                model.distribution.set_distribution(uncertainty, 'gauss', mean = 0.0, width = 1.0, range = (-4.0, 4.0))
        else: raise RuntimeError, "Line %d: unknown uncertainty type %s" % (lines[0][1], cmds[1])
    # add shape systematics:
    if '*' in shape_observables: shape_observables = set(channel_labels)
    data_done = set()
    for icol in range(n_cols):
        obs = channels_for_table[icol]
        if obs not in shape_observables: continue
        proc = processes_for_table[icol]
        found_matching_shapeline = False
        for l in shape_lines: # l = (process, channel, file, histogram, histogram_with_systematics)
            if l[1]!='*' and l[1]!=obs: continue
            if obs not in data_done and l[0] in ('*', 'data_obs', 'DATA'):
                add_shapes(model, obs, 'DATA', {}, l[2], l[3], '', include_mc_uncertainties)
                data_done.add(obs)
            if l[0]!='*' and l[0]!=proc: continue
            uncs = shape_systematics[obs].get(proc, {})
            #print "adding shapes for channel %s, process %s" % (obs, proc)
            add_shapes(model, obs, proc, uncs, l[2], l[3], l[4], include_mc_uncertainties)
            found_matching_shapeline = True
            break
        if not found_matching_shapeline:
            raise RuntimeError, "did not find a matching 'shapes' specification for channel '%s', process '%s'" % (obs, proc)
    model.set_signal_processes([transform_name_to_theta(proc) for proc in signal_processes])
    if include_mc_uncertainties: model.bb_uncertainties = True
    return model

