#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, os.path, datetime, re, sys, copy, traceback, shutil, hashlib, tempfile

from theta_auto import *

#TODO:
# * support for studies using +-1sigma toy / asimov data as input and re-run the method --> input='toys:scan-nuisance[-asimov]'?
# * support additional_nll_term everywhere:
#   - producers set additional-ll-term
#   - producers add additional terms in override-parameter-distribution
#  => make a base class Producer which contains takes care of these two!
# * switch Function to the one defined in FunctionBase



## \page theta_auto_intro Introduction to theta-auto
# 
#  The theta-auto.py script provides a set of functions to automatically generate theta configuration files, run theta on them, and analyze the output.
#  This page assumes you have some basic working knowledge of theta, concerning the configuration file syntax, and running theta.
#
#  With the theta-auto scripts, some further restrictions apply regarding the statistical model, compared to the capabilities of using theta directly.
#  The prototype model is a model with multiple channels where in each channel, the predicted events yield is the sum of the background and signal processes,
#  where each process is affected by shape uncertainties
#  modeled with cubiclinear_histomorph and rate uncertainties modeled via multiplicative log-normal priors. While some simple extensions to this prototype model
#  are possible (such as using truncated Gaussian priors instead of log-normal ones for the uncertainties, or using a different kind of histogram morphing),
#  this naturally restricts the use case to limit setting, discovery and measurement for a cross-section type parameter. On the other hand, this restriction
#  allows the python scripts to make many assumptions.
#
#  The easiest way to get started is probably to specify the model as a datacard as documented in
#  https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit
#  If you have such a datacard called "datacard.txt", put it in the same directory as a script called "analysis.py" with the content
#  \code
#  model = higgs_datacard.build_model('datacard.txt')
#  result = discovery(model)
#  print result
#  report.write_html('/home/username/public_html/datacard/')
#  \endcode
#  Then, from this same directory, call the theta-auto.py script (without any arguments).
#
#  This example shows the basic work flow of an analysis:
#  <ol>
#   <li>Create the external file(s) required to specify the model. This can be done using the Higgs group datacard format; other formats are discussed later</li>
#   <li>Create the "analysis.py" script for this analysis. This usually consists of three simple steps: (i) read in the model, (ii) call a statistical
#      method such as limit setting, discovery, etc, (iii) write output.</li>
#   <li>Call theta-auto.py. This will create (in a subdirectory called after the analysis script, in this case "analysis") the theta configuration files,
#       call theta, and analyze the result.</li>
#  </ol>
#
# Each function returns the result as a return value. The type of returned object depends on the function, see the function documentation
# for details. Most functions also writes a summary of the result to a global object called "report" which you can instruct to write out
# its content as html in a specified directory; this should be done at the very end of the "analysis.py" script.
#
# Do not execute the "analysis.py" script directly. Instead pass the name of the script as first argument to "theta-auto.py". In this case, also the created
# subdirectory, which contains the automatically produced theta configuration files and other intermediate output, will change to the
# script name without trailing ".py".
#
# The created subdirectory contains intermediate results such as the generated theta configuration files and sqlite files with the result.
# You can (and from time to time should!) delete it, everything will be re-generated. However, it makes sense to keep the cached result because
# theta-auto does not run theta again if there is a result created by theta based on the same configuration file. While this does not
# help when making changes to the model, it can be very helpful if the model and the method are set up and you mainly work on
# the extraction of the result.
#
# The most important statistical methods are:
# <ol>
#  <li>\link theta_auto::ml_fit ml_fit \endlink: make a maximum likelihood fit for all model parameters</li>
#  <li>\link theta_auto::pl_intervals pl_intervals \endlink: calculate profile likelihood intervals for the signal cross section; this includes a maximum likelihood fit</li>
#  <li>\link theta_auto::discovery discovery\endlink: calculate the distribution of a test statistic (usually likelihood ratio) for the background-only hypothesis to get the p-value /
#      number of standard deviations "sigma".</li>
#  <li>\link theta_auto::bayesian_quantiles bayesian_quantiles\endlink: calculate marginal posterior quantiles for the signal
#        strength parameter = upper / lower Bayesian limits on the signal cross section</li>
#  <li>\link theta_auto::posteriors posteriors \endlink: calculate the marginal posteriors for the given parameters (which can be either nuisance parameters or the signal cross section)</li>
#  <li>\link theta_auto::cls_limits cls_limits \endlink: calculate CLs limits, using a certain test statistic (usually a variant of profile likelihood ratio)</li>
#  <li>\link theta_auto::ks_test ks_test \endlink: calculate the KS p-value by generating toys, running a maximum likelihood fit, and calculating the
#    KS test statistic. This is useful to get KS p-values
#    which both (i) are valid after fitting data, (ii) include the effect of systematic uncertainties (both shape and rate).</li>
# </ol>
#
# In addition, there are some more, higher-level methods which typically internally call one or more of the above. For example, the function bayesian_limits
# calculates the limits for different signal hypotheses (such as signals for different mass points) and constructs a limit band plot (= a plot
# with median, +-1sigma and +-2sigma expected and observed limit) as a function of the signal process.
#
#
# \section theta_auto_model The Model in theta-auto
#
# The Model class contains all relevant information of the model, including the observed data, and all the predictions including their dependence on the
# model parameters in the different channels, and which of the processes is to be considered as signal.
# The signal is always added to the prediction of the background yields with a multiplicative signal strength factor called "beta_signal" which corresponds
# to the cross section in units of the signal cross section assumed in the input. All model parameters except "beta_signal" are called "nuisance parameters" throughout the documentation.
#
# Note that the python Model class in theta-auto is slightly different from the C++ theta::Model class you spefify in the %theta
# configuration files: theta::Model does not contain any information about data, nor about which process is considered as signal.
#
# The Model class also includes an instance of Distribution, model.distribution which is the prior distribution
# for all nuisance parameters. Toy data generation is always based on this prior. It is also used to add additional terms to the likelihood function, although
# this can be overridden by the parameter \c nuisance_constraint, which, if present, takes precedence over model.distribution; see \ref theta_auto_params below.
#
# An instance of the Model class is usually not constructed directly. Rather, use one of the functions which generate a model based on some external
# configuration like higgs_datacard.build_model or build_model_from_rootfile.
#
# It is possible to manipulate a Model instance with many functions. Important examples are to add a (log-normal) rate uncertainty for a certain process which
# can be done with
# \code
#  model.add_lognormal_uncertainty('ttbar_rate', math.log(1.12), 'ttbar', '*')
# \endcode
# which adds a 12% log-normal uncertainty controlled by the uncertainty 'ttbar_rate' on the process called 'ttbar', correlated across all channels ('*').
#
# Another important example is the combination of two Models which can be done via
# \code
# model1.combine(model2)
# \endcode
# which will add all channels of model2 to model1. In order for this to work correctly (=with correct correlations induced by shared uncertainties),
# model1 and model2 must have been built using the same convention for the names of nuisance parameters. Also, the signal processes
# should be the same for both models.
#
# For more information about Model manipulation, see the documentation of the Model class.
#
# \section theta_auto_signal Specifying what is 'signal'
#
# The Model and statistical methods in theta-auto are constructed such that it is easy to run them for situations in which there is more than one signal hypothesis, e.g.,
# a Zprime with different masses. In this case, the background model is always the same, and as signal model, one of the signal processes is chosen; the statistical
# procedure is applied for each of the signal processes in turn. For example, one might have the signal processes 'zp1000', 'zp2000', and 'zp3000' for Zprime models
# with different masses.
#
# In order for theta-auto to be able to correctly construct the background-only model, it must know what is signal. Use a call like
# \code
# model.set_signal_processes('zp*')
# \endcode
# to declare which processes exactly are the signal processes. Everything else is considered the background-only model.
#
# Sometimes, it makes sense to split the signal into several components which are conceptually the same process of interest
# (e.g., if the process was generated according to production process or decay or if the uncertainties for different signal components are different). For example,
# one might have split the production and processing of the Zprime signal sample according to production process into 'ggzp1000' and 'qqzp1000', etc., and in the example above
# with three masses one would have six signal processes in total.
# In such a case, one can specify more than one process name, all of which are scaled with the cross section scaling parameter "beta_signal". Such a group of processes to be
# used together at the same time is called "signal process group" and is specified via a list of strings. For different signal hypotheses, such as different masses,
# one has several such signal process groups. Therefore, a complete specification for which processes to run the statistical method on is a list of signal process groups.
# See also the documentation below of the \c signal_processes parameter. To run the statistical model on all three mass hypotheses from the above example, one would
# use
# \code
#  signal_processes = [['qqzp1000', 'ggzp1000'], ['qqzp2000', 'ggzp2000'], ['qqzp3000', 'ggzp3000']]
# \endcode
# For some purposes (such as a key in a dictionary return value or if used as part of a file name), identifying a signal process group via a list of
# strings is not possible. In this case, all process names of a signal process group are concatenated to for a "process group id". The process group
# ids in the above example are
# \code
#  'qqzp1000ggzp1000', 'qqzp2000ggzp2000', 'qqzp3000ggzp3000'
# \endcode
#
# \section theta_auto_params Common Parameters
#
# Almost all statistical methods allow to specify whether to run the method on the actually observed data or on toys. They also
# allow to specify alternative priors/likelihood terms (in theta, no distinction is made between these two). These
# parameters are:
# <ul>
#   <li>\c input: on which set of data to run the method. Can be either "data" or "toys:X" where X is the
#       signal scale factor; for example "toys:0" means no signal=background-only.</li>
#   <li>\c signal_prior: the signal prior to use for the likelihood function definition. Valid values are "fix:X", where X is the value to fix this parameter to, "flat" for
#     a flat prior on [0,inf] and "flat:[X,Y]" for a flat prior on the interval [X, Y]. For more complicated priors, you can directly specify the theta
#     configuration as dictionary.</li>
#   <li>\c nuisance_constraint: the signal prior to use for the likelihood function definition. Valid values are '' (=use model.distributions unmodified),
#     "shape:fix" (=fix shape-changing parameter, leave other parameter unchanged), "shape:free" (=use a flat prior for shape-changing parameters). Similarly, there are
#     "rate:fix", "rate:free", and combinations of "shape:..." and "rate:..." combinations, separated by semicolon, e.g., "shape:fixed;rate:free". For more complicated
#     situations, you can pass here an instance of the Distribution class which will be merged with the model.distribution (see Distribution.merge).</li>
#   <li>\c signal_processes: a list of signal process groups to considered together as signal. A signal process group is a list of strings of
#       the names of the signal processes. They all must have been declared as signal processes (via Model.set_signal_processes) in the model. The default (None)
#       will use all signals individually as a 'trivial' signal process group consisting of only this signal.</li>
# </ul>
#
# \section Some Internals
#
# This section describes some internals of theta-auto which are useful to know even if you do not want to modify theta-auto but only use it.
#
# The generated configuration file name follows the convention: <method>-<signal process group id>-<input>[-<id>]-<hash>.cfg
# where <method> is the python function name of the function generating the configuration file, <input>
# is the 'input' parameter as explained in the previous section, the optional <id> depends on the function, and <hash>
# is (part of a md5-)hash of the configuration file, in order to avoid clashes where otherwise the same values apply.
#
# Each configuration file creates one sqlite output file with the same name in the current directory; the suffix ".cfg"
# is replaced by ".db".


## \file theta-auto.py
#
# Main executable for the theta-auto scripts. See \ref theta_auto_intro for details.

## \namespace theta_auto
# \brief Python scripts of theta-auto
#
    
def clean_workdir():
    shutil.rmtree(config.workdir, ignore_errors = True)
    setup_workdir()
    
def setup_workdir():
    if not os.path.exists(config.workdir): os.mkdir(config.workdir)
    if not os.path.exists(os.path.join(config.workdir, 'plots')): os.mkdir(os.path.join(config.workdir, 'plots'))
            
def main():
    ROOT.gROOT.SetBatch(True)
    scriptname = 'analysis.py'
    tmpdir, profile = False, False
    wd = None
    quiet = False
    for iarg, arg in enumerate(sys.argv[1:]):
        if '.py' in arg: scriptname = arg
        if arg=='--tmpdir': tmpdir = True
        if arg=='--quiet': quiet = True
        if arg=='--profile': profile = True
        if arg=='--workdir': wd = sys.argv[iarg+2]
    if not quiet:
        print "*** WARNING ***"
        print "*** You are using an obsolete version of theta-auto (utils/theta-auto.py)."
        print "*** This version is no longer actively developed and only kept for backward compatibility. It will be removed in a future release."
        print "*** Please use instead version 2 (utils2/theta-auto.py)"
        print "*** (To get rid of this warning, call theta-auto with the '--quiet' option)"
    if tmpdir:
        config.workdir = tempfile.mkdtemp()
    else:
        config.workdir = os.path.join(os.getcwd(), scriptname[:-3])
        config.workdir = os.path.realpath(config.workdir)
        if wd is not None:
            config.workdir = os.path.realpath(wd)
    setup_workdir()
    config.theta_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    config.report = html_report(os.path.join(config.workdir, 'index.html'))
    variables = globals()
    variables['report'] = config.report
    utils.info("executing script %s" % scriptname)
    try:
        if profile:
            import cProfile
            cProfile.runctx("execfile(scriptname, variables)", globals(), locals())
        else:
            execfile(scriptname, variables)
    except Exception as e:
        print "error while trying to execute analysis script %s:" % scriptname
        traceback.print_exc()
        sys.exit(1)
    utils.info("workdir is %s" % config.workdir)
        
if __name__ == '__main__': main()



