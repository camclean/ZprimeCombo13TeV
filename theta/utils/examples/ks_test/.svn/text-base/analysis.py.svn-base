# Build a model using template morphing from the histograms in the specified root file.
# For details, see build_model_from_rootfile.
model = build_model_from_rootfile('templates.root')

# some bins predict zero while data might be non-zero. This leads to an infinite likelihood value
# which should be avoided. Therefore, make sure that all histogram bin contain at least some (small)
# positive value:
model.fill_histogram_zerobins()

# run the ks_test for each channel (=observable) independently by making toys which include shape
# uncertainties according to the model. fit_yield = True means that first, a maximum likelihood fit for
# the overall yield is done before calculating the KS test statistic. This effectively means to make
# a shape comparison only, not a rate comparison.
chan_to_pval = ks_test_individual_channels(model, fit_yield=True)

# ks_test_individual_channels returns a dictionary with the channel names as key and the p-value
# of the KS-test as value. Here, we just print the result for each channel:
for chan in chan_to_pval:
    print chan, chan_to_pval[chan]

