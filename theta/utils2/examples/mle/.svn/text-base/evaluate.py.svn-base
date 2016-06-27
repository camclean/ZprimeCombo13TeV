model = test_model.simple_counting(s = 2.0, b = 10.0, n_obs = 7, b_uncertainty = 3.0)
signal_process_groups = {'': []}
fit = mle(model, 'data', n = 1, signal_process_groups = signal_process_groups)
parameter_values = {}
for p in model.get_parameters([]):
    parameter_values[p] = fit[''][p][0][0]
histos = evaluate_prediction(model, parameter_values, include_signal = False)
write_histograms_to_rootfile(histos, 'histos-mle.root')


