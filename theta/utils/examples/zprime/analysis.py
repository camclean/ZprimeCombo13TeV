# for model building:
def get_model():
    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.
    model = build_model_from_rootfile('histograms.root')

    # If the prediction histogram is zero, but data is non-zero, teh negative log-likelihood
    # is infinity which causes problems for some methods. Therefore, we set all histogram
    # bin entries to a small, but positive value:
    model.fill_histogram_zerobins()

    # define what the signal processes are. All other processes are assumed to make up the 
    # 'background-only' model.
    model.set_signal_processes('zp*')

    # Add some lognormal rate uncertainties. The first parameter is the name of the
    # uncertainty (which will also be the name of the nuisance parameter), the second
    # is the 'effect' as a fraction, the third one is the process name. The fourth parameter
    # is optional and denotes the channl. The default '*' means that the uncertainty applies
    # to all channels in the same way.
    # Note that you can use the same name for a systematic here as for a shape
    # systematic. In this case, the same parameter will be used; shape and rate changes 
    # will be 100% correlated.
    model.add_lognormal_uncertainty('vjets_rate', 0.05, 'zjets')
    model.add_lognormal_uncertainty('vjets_rate', 0.05, 'wjets')
    model.add_lognormal_uncertainty('ttbar_rate', 0.15, 'ttbar')
    model.add_lognormal_uncertainty('ttbar_rate', 0.15, 'singletop')
    model.add_lognormal_uncertainty('qcd_rate', 1.0, 'qcd')
    # the qcd model is from data, so do not apply a lumi uncertainty on that:
    for p in model.processes:
        if p == 'qcd': continue
        model.add_lognormal_uncertainty('lumi', 0.045, p)
    # Specifying all uncertainties manually can be error-prone. You can also execute
    # a automatically generated file using python's execfile here
    # which contains these statements, or read in a text file, etc. Remember: this is a
    # python script, so use this power!
    model.add_lognormal_uncertainty('isrfsr', -0.068, 'ttbar', 'mu_htlep')
    model.add_lognormal_uncertainty('isrfsr', -0.064, 'ttbar', 'mu_mtt')
    model.add_lognormal_uncertainty('scale_ttbar', -0.12, 'ttbar', 'mu_htlep')
    model.add_lognormal_uncertainty('scale_ttbar', 0.092, 'ttbar', 'mu_mtt')
    model.add_lognormal_uncertainty('matching_ttbar', 0.061, 'ttbar', 'mu_htlep')
    model.add_lognormal_uncertainty('matching_ttbar', 0.080, 'ttbar', 'mu_mtt')
    model.add_lognormal_uncertainty('scale_vjets', -1.33, 'zjets', 'mu_htlep')
    model.add_lognormal_uncertainty('scale_vjets', -0.537, 'zjets', 'mu_mtt')
    model.add_lognormal_uncertainty('matching_vjets', 0.120, 'zjets', 'mu_htlep')
    model.add_lognormal_uncertainty('matching_vjets', -0.324, 'zjets', 'mu_mtt')
    model.add_lognormal_uncertainty('scale_vjets', -0.660, 'wjets', 'mu_htlep')
    model.add_lognormal_uncertainty('scale_vjets', -0.719, 'wjets', 'mu_mtt')
    model.add_lognormal_uncertainty('matching_vjets', 0.333, 'wjets', 'mu_htlep')
    model.add_lognormal_uncertainty('matching_vjets', 0.085, 'wjets', 'mu_mtt')
    return model

model = get_model()


# first, it is a good idea to generate a summary report to make sure everything has worked
# as expected. The summary will generate quite some information which should it make easy to spot
# errors like typos in the name of uncertainties, missing shape uncertaintie, etc.
model_summary(model)

# 2. apply the methods

# 2.a. Bayesian limits
# Calculate expected and observed Bayesian limits. For faster run time of this example,
# only make a few mass points. (Omitting the 'signal_procsses' parameter completely would
# process all signals defined as signal processes before; see Section "Common Parameters"
# on the theta auto intro doxygen page for details)
plot_exp, plot_obs = bayesian_limits(model, signal_processes = [['zp1000'], ['zp2000'], ['zp3000']])

# plot_exp and plot_obs are instances of plotutil.plotdata. they contain x/y values and
# bands. You can do many things with these objects such as inspect the x/y/ban
# data, pass then to plotutil.plot routine to make pdf plots, ...
# Here, we will just create text files of the plot data. This is useful if you want
# to apply your own plotting routines or present the result in a text Table.
plot_exp.write_txt('bayesian_limits_expected.txt')
plot_obs.write_txt('bayesian_limits_observed.txt')

# 2.b. CLs limits
# calculate cls limit plots. The interface is very similar to bayesian_limits. However, there are a few
# more options such as the definition of the test statistic which is usually a likelihood ratio but can differ in
# which parameters are minimized and which constraints / ranges are applied during minimization.
# Here, we stay with the default which fixes beta_signal=0
# for the background only hypothesis and lets it float freely for the signal+background hypothesis.
# See cls_limits documentation for more options.
plot_exp, plot_obs = cls_limits(model, signal_processes = [['zp1000'], ['zp2000'], ['zp3000']])

# as for the bayesian limits: write the result to a text file
plot_exp.write_txt('cls_limits_expected.txt')
plot_obs.write_txt('cls_limits_observed.txt')

# 3. compatibility tests:
# make a chi2 test by comparing the chi2 value distribution *after* fitting as obtained from background-only toys with
# a fit performed on data. The "signal_process_group" parameter is usually a list
# of the processes to consider together as the signal (which will be all scaled by beta_signal simultaneously).
# Here, we only include background and just use some arbitrary signal process, set the input beta_signal to zero ("toys:0.0")
# and also fix beta_signal to 0.0 during the fit ("fix:0.0").
p = chi2_test(model, signal_process_group = ['zp1000'], input = "toys:0.0", signal_prior = "fix:0.0")
print 'p-value from background-only chi2 test: ', p

# if one wants to include some signal in the comparison, one should make toys accordingly with signal, e.g.,
# with beta_signgl = 0.1 ("toys:0.1"). The signal_prior has a default of "flat" so this will now fit
# also the signal yield for each toy (and data) before computing chi2:
p = chi2_test(model, signal_process_group = ['zp1000'], input = "toys:0.1")
print 'p-value from chi2 test, including 0.1pb zprime at m=1TeV: ', p

# model_summary, bayesian_limits, and cls_limits also write their results to the 'report' object
# which we can ask to write its results as html page to a certain directory. Use an existing, empty
# directory and point your web browser to it.
report.write_html('htmlout')

# After running theta-auto, you probably want to delete the 'analysis' directory which
# contains intermediate results. Keeping it avoids re-running theta unnecessarily for unchanged configurations
# (e.g., because you just want to change the plot). However, this directory can grow very large over time.

