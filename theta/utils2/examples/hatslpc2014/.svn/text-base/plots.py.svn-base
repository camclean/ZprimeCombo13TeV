def p_to_color(proc):
    ps = {'ttbar': '#cc0000', 'wjets': '#00cc00', 'zjets': '#0000cc', 'singletop': '#cc00cc', 'DATA': '#000000',
          'zp1000': '#ff4444', 'zp2000': '#ff44ff', 'zp3000': '#4444ff', 'qcd': '#cccccc'}
    return ps.get(proc, '#000000')
    
def plot_model(model, fname_prefix = ''): # fname_prefix is used to construct the pdf filename of output plots
    for c in model.get_observables():
        backgrounds = []
        signals = []
        for p in model.get_processes(c):
            pd = plotdata(lw = 1.0, color = p_to_color(p))
            pd.set_histogram(model.get_histogram_function(c, p).get_nominal_histo())
            pd.legend = p
            if p in model.signal_processes:
                signals.append(pd)
            else:
                pd.fill_color = pd.color
                pd.color = '#000000'
                pd.lw = 0.5
                backgrounds.append(pd)
        make_stack(backgrounds)
        data = plotdata()
        data.set_histogram(model.get_data_histogram(c))
        # draw data with sqrt(n) errors and markers:
        data.draw_histo = False
        data.yerrors = [math.sqrt(n) for n in data.y] # NOTE: using sqrt(n) errors is not what the SC statistics committee recommends in this case ...
        data.marker = 'o'
        data.lw = 1.5
        data.markersize = 4.0
        # hack to place legend left for the htlep observable:
        legend_args = {'loc': "upper left" if 'htlep' in c else 'best'}
        plot(backgrounds + signals + [data], c, 'Events', '%s%s.pdf' % (fname_prefix, c), ymin = 0, legend_args = legend_args)


# make plots (as in the slides):
#model = build_model_from_rootfile('ex4input.root')
#model.set_signal_processes('zp*')
#plot_model(model)
