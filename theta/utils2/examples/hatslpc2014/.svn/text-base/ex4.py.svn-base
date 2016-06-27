execfile("common.py")

def get_model():
    model = build_model_from_rootfile('ex4input.root', include_mc_uncertainties = True)

    # define what the signal processes are. All other processes are assumed to make up the  'background-only' model.
    model.set_signal_processes('zp*')
    
    # add a minimum MC stat. uncertainty corresponding to +-1 MC event in each bins (esp. empty bins)
    model.fill_histogram_zerobins(None)
    
    # Specifying all uncertainties. Internally, this adds a factor exp(lambda * p)
    # where p is the parameter specified as first argument and lambda is the constant
    # in the second argument:
    model.add_lognormal_uncertainty('vjets_rate', 0.05, 'zjets')
    model.add_lognormal_uncertainty('vjets_rate', 0.05, 'wjets')
    model.add_lognormal_uncertainty('ttbar_rate', 0.15, 'ttbar')
    model.add_lognormal_uncertainty('ttbar_rate', 0.15, 'singletop')
    model.add_lognormal_uncertainty('qcd_rate', 1.0, 'qcd')
    
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
    
    # the qcd is derived from data, so do not apply a lumi uncertainty on that:
    for p in model.processes:
        if p == 'qcd': continue
        model.add_lognormal_uncertainty('lumi', 0.045, p)
    return model


model = get_model()


# 4.a.:
#model_summary(model)
#report.write_html('htmlout')



# 4.b.
# i.: maximum likelihood fit on data for each signal mass in turn. (See common.py for the implementation of mle_print)
#mle_print(model, 'data', 1, signal_process_groups = {'zp1000': ['zp1000'] })

# ii.: maximum likelihood fit on data without signal (background-only).
#mle_print(model, 'data', 1, signal_process_groups = {'bkgonly': [] })



# 4.c.
# i.
#result = zvalue_approx(model, 'data', 1)
#for signal in result:
#    print signal, result[signal]['Z']


# ii. test the Z-value distribution for toys:
#result = zvalue_approx(model, 'toys:1.0', 1000)
#for signal in result:
#    plot_histogram(result[signal]['Z'], 'histo-toys-Z-%s.pdf' % signal, logy = True, ymin = 0.5)

    
#4.d.

#expected, observed = asymptotic_cls_limits(model)
#print expected, observed
#report.write_html('htmlout')



#options = Options()
#options.set('main', 'n_threads', '10')
#expected, observed = bayesian_limits(model, options = options)
#print expected, observed
#report.write_html('htmlout')



# for the first bonus question in "extras":
#execfile("extras.py")
#plot_mle(model)
