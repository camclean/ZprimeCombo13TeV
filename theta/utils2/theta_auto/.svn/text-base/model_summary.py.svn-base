# -*- coding: utf-8 -*-
import config, os, math

from Report import *
from Model import *
import plotutil
import utils

def model_summary(model, create_plots = True, all_nominal_templates = False, shape_templates = False, lnmode = 'sym', plotargs = {}, dump_histos = False):
    """
    Write a html summary of the statistical model to the ``report`` object (see :ref:`report`).
    
    Parameters:
    
    * ``model`` - the statistical model to summarize
    * ``create_plots`` - if ``True``, create plots. Otehrwise, only create tables.
    * ``all_nominal_templates`` - if ``True``, create a separate plot for each process in each observable
    * ``shape_templates`` - if ``True``, create a plot with nominal/plus/minus templates for each triple (observable, process, uncertainty)
    * ``lnmode`` specifies how rate uncertainties are reported in the rate table. "1sigma" reports the rate change at +-1sigma of the nuisance parameter. This reports
      asymmetric uncertainties for uncertainties created with :meth:`Model.add_lognormal_uncertainty`. The alternative option "sym" reports the coefficient
      used in :meth:`Model.add_lognormal_uncertainty`. 
      
    The main result of this method is in the global report object. However, it also returns some summary tables
    as python object: The return value is a dictionary with two entries
    
    * 'rate_table' the rate table
    * 'sysrate_tables'  is a dictionary with the observable name as key. For each observable, it is a table summarizing the rate change for each process.
    """
    result = {}
    observables = sorted(list(model.observables))
    processes = sorted(list(model.processes))
    #general info
    f = open(os.path.join(config.workdir, "model_summary_general.thtml"), 'w')
    print >>f, '<p>Observables (xmin, xmax, nbins):</p>\n<ul>'
    for o in observables:
        xmin, xmax, nbins = model.observables[o]
        print >> f, '<li>%s (%.5g, %.5g, %d)</li>' % (o, xmin, xmax, nbins)
    print >> f, "</ul>"
    print >>f, '<p>Background processes:</p>\n<ul>'
    for p in processes:
        if p in model.signal_processes: continue
        print >> f, '<li>%s</li>' % p
    print >> f, "</ul>"
    print >> f, "<p>Signal processes:</p><ul>"
    for p in processes:
        if p in model.signal_processes: print >> f, '<li>%s</li>' % p
    print >> f, "</ul>"
    print >> f, '<p>Nuisance parameters (includes only those which apply to the background-only model):</p>\n<ul>'
    parameters = sorted(list(model.get_parameters([])))
    rc_pars, sc_pars = model.get_rate_shape_parameters()
    for par in parameters:
        rate, morph = par in rc_pars, par in sc_pars
        if rate and morph: suffix = '(morph and rate)'
        elif morph: suffix = '(morph only)'
        elif rate: suffix = '(rate only)'
        print >> f, '<li>%s %s</li>' % (par, suffix)
    print >> f, "</ul>"
    f.close()
    
    # rate table
    default_parameters = model.distribution.get_means()
    rate_table = table()
    rate_table.add_column('process', 'process / observable')
    o_bkg_sum = {}
    o_bkg_err2sum = {}
    for o in observables:
        rate_table.add_column(o)
        o_bkg_sum[o] = 0.0
        o_bkg_err2sum[o] = 0.0
    for p in processes:
        if p in model.signal_processes: continue
        rate_table.set_column_multiformat('process', p)
        for o in observables:
           hf = model.get_histogram_function(o,p)
           if hf is None:
               rate_table.set_column_multiformat(o, '---')
               continue
           hf0 = hf.evaluate(default_parameters)
           s = hf0.get_value_sum()
           uncertainties = hf0.get_uncertainties()
           if uncertainties is not None: error = math.sqrt(sum([x**2 for x in uncertainties]))
           else: error = 0.0
           o_bkg_sum[o] += s
           o_bkg_err2sum[o] += error**2
           if error > 0:  rate_table.set_column_multiformat(o, (s, error), html = '%.5g +/- %.5g' % (s, error), tex = '$%.5g \\pm %.5g$' % (s, error))
           else: rate_table.set_column_multiformat(o, s, html = '%.5g' % s, tex = '%.5g' % s)
        rate_table.add_row()
    rate_table.set_column_multiformat('process', 'bkg_tot', html = '<b>total background</b>', tex = r'\textbf{total background}')
    for o in observables:
        error = math.sqrt(o_bkg_err2sum[o])
        s = o_bkg_sum[o]
        if error > 0: rate_table.set_column_multiformat(o, (s, error), html = '%.5g +/- %.5g' % (s, error), tex = '$%.5g \\pm %.5g$' % (s, error))
        else: rate_table.set_column_multiformat(o, s, html = '%.5g' % s, tex = '%.5g' % s)
    rate_table.add_row()
    for p in processes:
        if p not in model.signal_processes: continue
        rate_table.set_column_multiformat('process', p)
        for o in observables:
           hf = model.get_histogram_function(o,p)
           if hf is None:
               rate_table.set_column_multiformat(o, '---')
               continue
           s = sum(hf.get_nominal_histo()[2])
           error = 0
           uncertainties = hf.get_nominal_histo().get_uncertainties()
           if uncertainties is not None: error = math.sqrt(sum([x**2 for x in uncertainties]))
           if error > 0:  rate_table.set_column_multiformat(o, (s, error), html = '%.5g +/- %.5g' % (s, error), tex = '$%.5g \\pm %.5g$' % (s, error))
           else: rate_table.set_column_multiformat(o, s, html = '%.5g' % s, tex = '%.5g' % s)
        rate_table.add_row()
    #always show data row (even if there is no DATA ...):
    rate_table.set_column_multiformat('process', 'data', html = '<b>DATA</b>', tex = r'\textbf{data}')
    for o in observables:
        histo = model.get_data_histogram(o)
        if histo is None: rate_table.set_column_multiformat(o, '---')
        else: rate_table.set_column_multiformat(o, sum(histo[2]), html = '%.5g' % sum(histo[2]), tex = '%.5g' % sum(histo[2]))
        
    rate_table.add_row()
    f = open(os.path.join(config.workdir, "model_summary_rates.thtml"), 'w')
    print >> f, "<p>Rates for all observables and processes as given by the 'nominal' templates. If errors are given, they are MC stat. uncertainties.</p>"
    print >> f, rate_table.html()
    result['rate_table'] = rate_table
    f.close()
    
    result['sysrate_tables'] = {}
    
    # rate impact of systematic uncertainties on the processes. One table per observable:
    f = open(os.path.join(config.workdir, "model_summary_rate_impacts.thtml"), 'w')
    print >> f, """<p>The table below summarises the impact of an nuisance parameter on the rate prediction of a process.</p>
    <p>For a nuisance parameter, (gauss) indicates that this nuisance parameter has a gaussian prior, (gamma) that it has a gamma prior.</p>
    <p>For the individual cells, (r) indicates the 'rate only' part of the uncertainty, (s) indicates the effect on the rate of an uncertainty
    treated via template morphing (i.e., the <em>rate</em> effect of an uncertainty treated as part of the
    template morphing; even if this is zero, the shape effect is still taken into account). Note that both effects are applied seperatly, so
    the total rate change is about the linear sum of these two.<br/>
    The rate change in 'plus' direction of the uncertainty is written as superscript,
    the 'minus' direction as subscript.<br/>All numbers are in percent.</p>"""
    for o in observables:
        print >> f, "<h2>Observable '%s'</h2>" % o
        rate_impact_table = table()
        rate_impact_table.add_column('process', 'process / nuisance parameter')
        parameters = sorted(list(model.get_parameters(list(model.signal_processes))))
        parameters = [p for p in parameters if p!='beta_signal'] 
        for par in parameters:
            d = model.distribution.get_distribution(par)
            rate_impact_table.add_column(par, '%s (%s)' % (par, d['typ']))
        for p in processes:
            rate_impact_table.set_column_multiformat('process', p)
            hf = model.get_histogram_function(o,p)
            if hf is None: continue
            coeff = model.get_coeff(o,p)
            histo_nominal_integral = hf.evaluate(default_parameters).get_value_sum()
            if histo_nominal_integral == 0.0: histo_nominal_integral = float("nan")
            for par in parameters:
                splus, sminus = None, None
                if par in hf.get_parameters():
                    parameters_iplus = dict(default_parameters)
                    parameters_iplus[par] += model.distribution.get_distribution(par)['width']
                    parameters_iminus = dict(default_parameters)
                    parameters_iminus[par] -= model.distribution.get_distribution(par)['width']
                    # set splus and sminus to relative change in rate:
                    splus = hf.evaluate(parameters_iplus).get_value_sum() / histo_nominal_integral - 1.0
                    sminus = hf.evaluate(parameters_iminus).get_value_sum() / histo_nominal_integral - 1.0
                rplus, rminus = 0.0, 0.0
                if par in coeff.factors:
                    if type(coeff.factors[par])==dict and coeff.factors[par]['type'] == 'exp_function':
                        rplus = math.exp(coeff.factors[par]['lambda_plus']) - 1
                        if lnmode == '1sigma': rminus = math.exp(-coeff.factors[par]['lambda_minus']) - 1
                        else: rminus = 1 - math.exp(coeff.factors[par]['lambda_minus'])
                    else:
                       rplus = model.distribution.get_distribution(par)['width']
                       rminus = -model.distribution.get_distribution(par)['width']
                cell = ''
                texcell = ''
                rawcell = {}
                if splus is not None:
                    cell += '<sup>%+.2f</sup><sub>%+.2f</sub> (s) ' % (splus * 100, sminus * 100)
                    texcell = '$^{%+.2f}_{%+.2f}$ (s) ' % (splus * 100, sminus * 100)
                    rawcell['shape_plus'] = splus * 100
                    rawcell['shape_minus'] = sminus * 100
                if (rplus, rminus) != (0.0, 0.0):
                    rawcell['rate_plus'] = (rplus * 100)
                    rawcell['rate_minus'] = (rminus * 100)
                    if rplus==-rminus:
                        if rplus > 0:
                            cell += '&#xb1;%.2f (r)' % (rplus * 100)
                            texcell += "$\\pm %.2f$ (r)" % (rplus * 100)
                        else:
                            cell += '&#x2213;%.2f (r)' % (-rplus * 100)
                            texcell += "$\\mp %.2f$ (r)" % (-rplus * 100)
                    else:
                        cell += '<sup>%+.2f</sup><sub>%+.2f</sub> (r) ' % (rplus * 100, rminus * 100)
                        texcell += '$^{%+.2f}_{%+.2f}$ (r) ' % (rplus * 100, rminus * 100)
                if cell == '': cell, texcell = '---', '---'
                rate_impact_table.set_column_multiformat(par, rawcell, html = cell, tex = texcell)
            rate_impact_table.add_row()
        print >> f, rate_impact_table.html()
        result['sysrate_tables'][o] = rate_impact_table
    f.close()
    
    # nuisance parameter priors
    model_summary_nuisance(model.distribution, os.path.join(config.workdir, "model_summary_nuisance.thtml"))

    # figures:
    if create_plots: model_plots(model, all_nominal_templates = all_nominal_templates, shape_templates = shape_templates, plotargs = plotargs, dump_histos = dump_histos)

    # build the index.html:
    config.report.new_section('General Model Info', file(os.path.join(config.workdir, 'model_summary_general.thtml')).read())
    config.report.new_section('Rate Summary', file(os.path.join(config.workdir, 'model_summary_rates.thtml')).read())
    config.report.new_section('Rate Impact of Systematic Uncertainties', file(os.path.join(config.workdir, 'model_summary_rate_impacts.thtml')).read())
    
    nuisance_priors =  """
    <p>The priors for the nuisance parameters are either Gaussian or gamma distributions. As limit cases, these can have with=0 or width=inf
    which makes them delta or flat distributions, respectively.</p>
    <div class="inner">%(model_summary_nuisance.thtml)s</div>
    </div>"""
    fnames = ['model_summary_nuisance.thtml']
    d = {}
    for fname in fnames:
        f = open(os.path.join(config.workdir, fname), 'r')
        d[fname] = f.read()
        f.close()
    nuisance_priors = nuisance_priors % d
    config.report.new_section('Nuisance Parameter Priors', nuisance_priors)
    if create_plots:
        config.report.new_section('Basic Model Plots', file(os.path.join(config.workdir, 'model_plots.thtml')).read())
        
    return result

def pretty_dict(d):  return '; '.join(['%s = %s' % (k, str(d[k])) for k in d])


def model_summary_nuisance(dist, fname):
    f = open(fname, 'w')
    print >> f, "<h3>Prior Parameters</h3>"
    t = table()
    t.add_column('parameter')
    t.add_column('type', 'distribution type')
    t.add_column('pars', 'distribution parameters')
    for p in dist.get_parameters():
       t.set_column('parameter', p)
       d = dist.get_distribution(p)
       t.set_column('type', d['typ'])
       del d['typ']
       t.set_column('pars', pretty_dict(d))
       t.add_row()
    print >> f, t.html()
    f.close()
    
# creates plots at certain parameter values
def model_plots_at(model, par_values, signal_stacked = False):
    plotdir = os.path.join(config.workdir, 'plots')
    processes = sorted(list(model.processes))
    #TODO: more / better colors
    h = str(hash(str(par_values)))
    background_colors = ['#edd400', '#f57900', '#c17d11', '#73d216', '#3465a4', '#75507b', '#d3d7cf', '#555753']
    signal_colors = ['#ef2929', '#cc0000', '#a40000']
    if not os.path.exists(plotdir): os.mkdir(plotdir)
    text = '<p>Templates evaluated at:<p><ul>'
    for p in par_values:
        text+='<li>%s = %.2g</li>\n' % (p, par_values[p])
    text += '</ul>'
    text += "<p>Everything normalized to expectation, i.e., to the normalization in the template input file, possibly scaled via the python script file.</p>"
    text += "<p>Color Code:</p><ul>"
    i_bkg_col = 0
    i_signal_col = 0
    for p in processes:
        if p in model.signal_processes:
            color = signal_colors[i_signal_col]
            i_signal_col = (i_signal_col + 1) % len(signal_colors)
        else:
            color = background_colors[i_bkg_col]
            i_bkg_col = (i_bkg_col + 1) % len(background_colors)
        text += '<li><span style="background: %s;">&nbsp;&nbsp;&nbsp;</span> %s</li>' % (color, p)
    text += '</ul>'
    templates = get_shifted_templates(model, par_values)
    for o in templates:
        background_pds = []
        signal_pds = []
        i_bkg_col = 0
        i_signal_col = 0
        for p in templates[o]:
            pd = plotutil.plotdata()
            pd.histo_triple(templates[o][p])
            xmin, xmax, data = templates[o][p]
            binwidth = (xmax - xmin) / len(data)
            if p in model.signal_processes:
                pd.color = signal_colors[i_signal_col]
                i_signal_col = (i_signal_col + 1) % len(signal_colors)
                signal_pds.append(pd)
            else:
                pd.fill_color = background_colors[i_bkg_col]
                pd.lw = 1
                pd.color = '#000000'
                i_bkg_col = (i_bkg_col + 1) % len(background_colors)
                background_pds.append(pd)
        data_histo = model.get_data_histogram(o)
        data_pd = None
        if data_histo is not None:
            xmin, xmax, data = data_histo
            data_pd = plotutil.plotdata()
            data_pd.color = '#000000'
            data_pd.histo_triple(data_histo)
            data_pd.yerrors = map(math.sqrt, data_pd.y)
            data_pd.circle = 'o'
        plots = background_pds
        if signal_stacked:
            plots.extend(signal_pds)
            plotutil.make_stack(background_pds)
        else:
            plotutil.make_stack(background_pds)
            plots.extend(signal_pds)
        if data_pd is not None: plots.append(data_pd)
        plotutil.plot(plots, o, '$N / %.4g$' % binwidth, os.path.join(plotdir, '%s_stack%s.png' % (o, h)), xmin=xmin, xmax=xmax)
        text += "<p>Observable '%s':<br /><img src=\"plots/%s_stack%s.png\" /></p>" % (o, o, h)
    config.report.new_section('Model Plots at parameter values', text)
    

# creates plots and model_plots.thtml
# observable units is a dictionary (observable name) --> (caption to use for plots)
def model_plots(model, all_nominal_templates = False, shape_templates = False, plotargs = {}, dump_histos = False):
    plotdir = os.path.join(config.workdir, 'plots')
    observables = sorted(list(model.observables.keys()))
    processes = sorted(list(model.processes))
    #TODO: more / better colors
    background_colors = ['#edd400', '#f57900', '#c17d11', '#73d216', '#3465a4', '#75507b', '#d3d7cf', '#555753']
    signal_colors = ['#ef2929', '#cc0000', '#a40000']
    if not os.path.exists(plotdir): os.mkdir(plotdir)
    f = open(os.path.join(config.workdir, 'model_plots.thtml'), 'w')
    print >> f, "<h2>Stackplots</h2>"
    print >> f, "<p>Everything normalized to expectation, i.e., to the normalization in the template input file, possibly scaled via the python script file.</p>"
    print >> f, "<p>Color Code:</p><ul>"
    i_bkg_col = 0
    i_signal_col = 0
    for p in processes:
        if p in model.signal_processes:
            color = signal_colors[i_signal_col]
            i_signal_col = (i_signal_col + 1) % len(signal_colors)
        else:
            color = background_colors[i_bkg_col]
            i_bkg_col = (i_bkg_col + 1) % len(background_colors)
        print >> f, '<li><span style="background: %s;">&nbsp;&nbsp;&nbsp;</span> %s</li>' % (color, p)
    print >>f, '</ul>'
    default_parameters = model.distribution.get_means()
    for o in observables:
        background_pds = []
        signal_pds = []
        i_bkg_col = 0
        i_signal_col = 0
        for p in processes:
            hf = model.get_histogram_function(o, p)
            if hf is None: continue
            pd = plotutil.plotdata()
            h0 = hf.evaluate(default_parameters)
            pd.set_histogram(h0)
            xmin, xmax, data = h0
            #binwidth = (xmax - xmin) / len(data)
            if p in model.signal_processes:
                pd.color = signal_colors[i_signal_col]
                i_signal_col = (i_signal_col + 1) % len(signal_colors)
                signal_pds.append(pd)
            else:
                pd.fill_color = background_colors[i_bkg_col]
                pd.lw = 1
                pd.color = '#000000'
                i_bkg_col = (i_bkg_col + 1) % len(background_colors)
                background_pds.append(pd)
        data_histo = model.get_data_histogram(o)
        data_pd = None
        if data_histo is not None:
            data_pd = plotutil.plotdata()
            data_pd.color = '#000000'
            data_pd.histo_triple(data_histo)
            data_pd.yerrors = map(math.sqrt, data_pd.y)
            data_pd.circle = 'o'
        plotutil.make_stack(background_pds)
        plots = background_pds + signal_pds
        if data_pd is not None: plots.append(data_pd)
        plotutil.plot(plots, o, '$N$', os.path.join(plotdir, '%s_stack.png' % o), xmin=xmin, xmax=xmax, **plotargs)
        print >> f, "<p>Observable '%s':<br /><img src=\"plots/%s_stack.png\" /></p>" % (o, o)
       
    if all_nominal_templates:
        if dump_histos: df = open('dump.txt', 'w')
        print >> f, "<h2>All 'nominal' Templates</h2>"
        print >> f, "<p>Everything normalized to expectation, i.e., to the normalization in the template input file, possibly scaled via the python script file.</p>"
        for o in observables:
            for p in processes:
                hf = model.get_histogram_function(o,p)
                if hf is None: continue
                xmin, xmax, data = hf.get_nominal_histo()
                binwidth = (xmax - xmin) / len(data)
                pd = plotutil.plotdata()
                pd.x = [xmin + i*binwidth for i in range(len(data))]
                pd.y = data[:]
                pd.color = signal_colors[0]
                xlabel = o
                plotutil.plot([pd], xlabel, '$N / %.4g$' % binwidth, os.path.join(plotdir, '%s_%s.png' % (o, p)), xmin=xmin, xmax=xmax, **plotargs)
                print >> f, '<p>Observable "%s", Process "%s":<br/><img src="plots/%s_%s.png"/></p>' % (o, p, o, p)
            # make also one plot with all signal processes, and normalization versus ordering:
            pd_norm = plotutil.plotdata()
            pd_norm.x = []
            pd_norm.y = []
            pd_norm.as_histo = False
            pd_norm.color = '#000000'
            plots = []
            i_signal_col = 0
            x_to_y = {}
            for p in processes:
                if p not in model.signal_processes: continue
                hf = model.get_histogram_function(o,p)
                if hf is None: continue
                if dump_histos:
                    h = hf.get_nominal_histo()
                    df.write('\n')
                    df.write('%s %s\n' % (o, p))
                    df.write(' '.join(['%.3f' % n for n in h.get_values()]))
                    df.write('\n')
                    if h.get_uncertainties() is not None:
                        df.write(' '.join(['%.3f' % n for n in h.get_uncertainties()]))
                        df.write('\n')
                xmin, xmax, data = hf.get_nominal_histo()
                x_to_y[utils.extract_number(p)] = sum(data)
                binwidth = (xmax - xmin) / len(data)
                pd = plotutil.plotdata()
                pd.x = [xmin + i*binwidth for i in range(len(data))]
                pd.y = data[:]
                pd.color = signal_colors[i_signal_col]
                plots.append(pd)
                i_signal_col = (i_signal_col + 1) % len(signal_colors)
            for x in sorted(x_to_y.keys()):
                pd_norm.x.append(x)
                pd_norm.y.append(x_to_y[x])
            plotutil.plot(plots, o, '$N / %.4g$' % binwidth, os.path.join(plotdir, '%s_signals.png' % o), xmin=xmin, xmax=xmax, **plotargs)
            plotutil.plot([pd_norm], 'signal process', '$N$', os.path.join(plotdir, '%s_norm_vs_signals.png' % o), **plotargs)
            print >> f, '<p>Observable "%s", all signals: <br/><img src="plots/%s_signals.png"/></p>' % (o, o)
            print >> f, '<p>Observable "%s", signal normalization: <br/><img src="plots/%s_norm_vs_signals.png"/></p>' % (o, o)
    # (end if all_nominal_templates)
    if not shape_templates:
        f.close()
        return 
    # shape comparison for morphed templates:
    color_nominal, color_plus, color_minus = '#333333', '#aa3333', '#3333aa'
    print >> f, "<h2>Shape Uncertainty Plots</h2>"
    print >> f, "<p>Color Code:</p><ul>"
    print >> f, "<li><span style=\"background: %s;\">&nbsp;&nbsp;&nbsp;</span> nominal</li><li><span style=\"background: %s;\">&nbsp;&nbsp;&nbsp;</span> plus</li><li><span style=\"background: %s;\">&nbsp;&nbsp;&nbsp;</span> minus</li>" % (color_nominal, color_plus, color_minus)
    print >> f, "</ul>"
    print >> f, "<p>Processes not appearing in the tables do not have any shape uncertainty for this observable.</p>"
    print >> f, "<p>Click on an image to enlarge. If you have javascript, the image will be displayed on this page and you can click through all shape uncertainties of that observable \
                   (instead of clicking, you can also use the left/right key on the keyboard).</p>"
    for o in observables:
        print >> f, '<h3>Observable \'%s\'</h3>' % o
        # save the triples (o,p,u) for which there is a plot:
        opus = []
        for p in model.get_processes(o):
            hf = model.get_histogram_function(o,p)
            for u in hf.syst_histos:
                xmin, xmax, data_nominal = hf.nominal_histo
                xmin, xmax, data_plus = hf.syst_histos[u][0]
                xmin, xmax, data_minus = hf.syst_histos[u][1]
                binwidth = (xmax - xmin) / len(data_nominal)
                pd = plotutil.plotdata(color = color_nominal, legend = 'nominal')
                pd.x = [xmin + i*binwidth for i in range(len(data_nominal))]
                pd.y = data_nominal
                pd_plus = plotutil.plotdata(color = color_plus, legend = 'plus variation')
                pd_plus.x = pd.x
                pd_plus.y = data_plus
                pd_plus.color = color_plus
                pd_minus = plotutil.plotdata(color = color_minus, legend = 'minus variation')
                pd_minus.x = pd.x
                pd_minus.y = data_minus
                name = '%s__%s__%s' % (o,p,u)
                plotutil.plot((pd, pd_plus, pd_minus), o, '$N / %.4g$' % binwidth, os.path.join(plotdir, name + '.png'), xmin=xmin, xmax = xmax)
                opus.append((o,p,u))
        #make a table for this observable:
        t = table()
        t.add_column('process', 'process / uncertainty')
        us = sorted(list(set([u for o,p,u in opus])))
        ps = sorted(list(set([p for o,p,u in opus])))
        for u in us: t.add_column(u)
        for p in ps:
            t.set_column('process', p)
            for u in us:
                if (o,p,u) in opus:
                    t.set_column(u, '<a href="plots/%s__%s__%s.png" rel="lightbox[%s]"><img src="plots/%s__%s__%s.png" width="200"/></a>' % (o,p,u,o,o,p,u))
                else:
                    t.set_column(u, '---')
            t.add_row()
        print >>f, t.html()
        if len(opus)==0: print >> f, '<p>no shape uncertainties for this observable</p>'
    f.close()


