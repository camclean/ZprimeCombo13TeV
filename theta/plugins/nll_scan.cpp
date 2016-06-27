#include "plugins/nll_scan.hpp"
#include "plugins/reduced_nll.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>

using namespace theta;
using namespace std;

void nll_scan::produce(const Data & data, const Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(not start_step_ranges_init){
        const Distribution & d = nll->get_parameter_distribution();
        m_ranges.set_from(d);
        d.mode(m_start);
        m_step.set(asimov_likelihood_widths(model, override_parameter_distribution));
        start_step_ranges_init = true;
    }
    MinimizationResult minres = minimizer->minimize(*nll, m_start, m_step, m_ranges);
    double fval0 = minres.fval;
    ParValues values0(minres.values);
    products_sink->set_product(c_maxl, minres.values.get(pid));
    
    theta::Histogram1D result(n_steps, start, start + n_steps * step);
    if(adaptive_startvalues && re_minimize){
        // instead of going through all x-values, start at the minimum:
        double x0 = minres.values.get(pid);
        unsigned int i0 = max<int>((x0 - start) / step, 0);
        if(i0 >= n_steps) i0 = n_steps-1;
        ParValues tmp_start(m_start);
        ParValues tmp_step(m_step);
        tmp_step.set(pid, 0.0);
        Ranges tmp_ranges(m_ranges);
        
        tmp_start.set(values0);
        for(unsigned int i = i0; i < n_steps; ++i){
            double x = start + i * step;
            tmp_ranges.set(pid, make_pair(x, x));
            tmp_start.set(pid, x);
            minres = minimizer->minimize(*nll, tmp_start, tmp_step, tmp_ranges);
            result.set(i, minres.fval - fval0);
            tmp_start.set(minres.values);
        }
        tmp_start.set(values0);
        for(int i = i0 - 1; i>=0; --i){
            double x = start + i * step;
            tmp_ranges.set(pid, make_pair(x, x));
            tmp_start.set(pid, x);
            minres = minimizer->minimize(*nll, tmp_start, tmp_step, tmp_ranges);
            result.set(i, minres.fval - fval0);
            tmp_start.set(minres.values);
        }
    }
    else{
        ReducedNLL nll_r(*nll, pid, minres.values, re_minimize ? minimizer.get() : 0, m_start, m_step, m_ranges);
        nll_r.set_offset_nll(minres.fval);
        for(unsigned int i=0; i<n_steps; ++i){
            double x = start + i * step;
            result.set(i, nll_r(x));
        }
    }
    products_sink->set_product(c_nll, result);
}


nll_scan::nll_scan(const theta::Configuration & cfg): Producer(cfg), pid(get_parameter(cfg, "parameter")),
   re_minimize(true), adaptive_startvalues(false), start_step_ranges_init(false){
    Setting s = cfg.setting;
    minimizer = PluginManager<Minimizer>::build(theta::Configuration(cfg, s["minimizer"]));
    if(s.exists("re-minimize")){
        re_minimize = s["re-minimize"];
    }
    if(s.exists("adaptive_startvalues")){
        adaptive_startvalues = s["adaptive_startvalues"];
    }
    start = s["parameter-values"]["start"];
    stop = s["parameter-values"]["stop"];
    n_steps = s["parameter-values"]["n-steps"];
    if(n_steps<2){
        throw ConfigurationException("nll_scan: n-steps must be >= 2");
    }
    if(start >= stop){
        throw ConfigurationException("nll_scan: start < stop must hold");
    }
    step = (stop - start) / n_steps;
    c_nll = products_sink->declare_product(*this, "nll", theta::typeHisto);
    c_maxl = products_sink->declare_product(*this, "maxl", theta::typeDouble);
}

REGISTER_PLUGIN(nll_scan)
