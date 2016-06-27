execfile("plots.py")

def plot_mle(model):
    # see http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/theta-auto/examples.html#get-the-histograms-at-the-maximum-likelihood-fit
    fit = mle(model, 'data', n = 1, signal_process_groups = {'bkgonly': []}) ['bkgonly']
    parameter_values = {}
    for p in model.get_parameters([]):
        parameter_values[p] = fit[p][0][0]
    parameter_values['beta_signal'] = 1.0 # do not scale signal
    histos = evaluate_prediction(model, parameter_values, include_signal = True)
    write_histograms_to_rootfile(histos, 'histos-mle.root')
    # build the model from the scaled/morphed histograms:
    model2 = build_model_from_rootfile('histos-mle.root')
    model2.set_signal_processes('zp*')
    # note that model2 does not contain data, so use data from model:
    for o in model.get_observables():
        model2.set_data_histogram(o, model.get_data_histogram(o))
    plot_model(model2, fname_prefix = "fitted-")


# a toy-based look-elsewhere effect:
def lee(model):
    # make N background-only toy data:
    N = 1000
    bkgonly = make_data(model, 'toys:0.0', N, signal_process_groups = {'bkgonly': [] })['bkgonly']
    # calculate the (asymptotic) local Z values for all signals, using the generated bkg-only toys as input:
    res0 = zvalue_approx(model, bkgonly, N, eventid_info = True)
    # IMPORTANT: sometimes, fits do not converge for a toy, in which case this output is missing in the res0 list.
    #  To compare really the same background-only toys, we have to use the "eventid" assigned by theta to get this right in that case.
    res = {}
    for signal in res0:
        eventid_to_Z = dict([(e, Z) for e,Z in zip(res0[signal]['__eventid'], res0[signal]['Z'])])
        res[signal] = eventid_to_Z
    # now determine the distribution for the *maximum* Z-value in each toy, maximizing over all signals:
    Z_max_list = []
    for i in range(N):
        Z_max = max([res[signal].get(i, 0.0) for signal in res])
        Z_max_list.append(Z_max)
    plot_histogram(Z_max_list, 'Z_max.pdf')
    # TODO: add your code here to determine how likely it is to see a Z_max >= 2

    
