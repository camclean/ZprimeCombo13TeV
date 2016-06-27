#include "plugins/deltanll_intervals.hpp"
#include "plugins/reduced_nll.hpp"
#include "interface/secant.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/minimizer.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>
#include <iomanip>

using namespace theta;
using namespace std;

void deltanll_intervals::produce(const theta::Data & data, const theta::Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(not start_step_ranges_init){
        const Distribution & d = nll->get_parameter_distribution();
        ranges.set_from(d);
        d.mode(start);
        try{
            step.set(asimov_likelihood_widths(model, override_parameter_distribution));
        }
        catch(const Exception & ex){
            throw logic_error("could not initiaize step width with asimov likelihood: " + ex.message);
        }
        start_step_ranges_init = true;
    }
    MinimizationResult minres = minimizer->minimize(*nll, start, step, ranges);
    const double value_at_minimum = minres.values.get(pid);
    products_sink->set_product(c_maxl, value_at_minimum);
    bool do_minimization = re_minimize && nll->get_parameters().size() > 1;
    ReducedNLL nll_r(*nll, pid, minres.values, do_minimization ? minimizer.get() : 0, start, step, ranges);
    const pair<double, double> & range = ranges.get(pid);
    for(size_t i=0; i < deltanll_levels.size(); ++i){
        nll_r.set_offset_nll(minres.fval + deltanll_levels[i]);
        //upper value: look for a parameter value with a sign flip:
        double x_low = value_at_minimum;
        double f_x_low = -deltanll_levels[i];
        double initial_step = minres.errors_plus.get(pid);
        if(initial_step <= 0) initial_step = step.get(pid);
        else if(initial_step <= 1e-6 * fabs(x_low)) initial_step = 1e-6 * fabs(x_low);
        if(initial_step < 1e-6) initial_step = 1e-6;
        double step = initial_step;
        const double x_acurracy = step / 100;
        double x_high, f_x_high;
        int k;
        const int k_max = 20;
        for(k=1; k<=k_max; ++k){
            x_high = value_at_minimum + step;
            step *= 2;
            if(x_high > range.second){
                x_high = range.second;
            }
            f_x_high = nll_r(x_high);
            if(f_x_high > 0){
                products_sink->set_product(upper_columns[i], secant(x_low, x_high, x_acurracy, f_x_low, f_x_high, deltanll_levels[i]/1000, nll_r));
                break;
            }
            else if(f_x_high==0.0 || x_high == range.second){
                products_sink->set_product(upper_columns[i], x_high);
                break;
            }
            else{
                x_low = x_high;
                f_x_low = f_x_high;
            }
        }
        if(k > k_max){
            stringstream ss;
            ss << "Could not find upper value for interval after " << k_max << " iterations (min at: " << value_at_minimum << "; initial step: " << initial_step << "; current step: " << step << ")";
            throw Exception(ss.str());
        }
        
        //lower value: same story, just other way round:
        x_high = value_at_minimum;
        f_x_high = -deltanll_levels[i];
        step = initial_step;
        if(step <= 0.0) step = 1.0;
        for(k=1; k <= k_max; ++k){
            x_low = value_at_minimum - step;
            step *= 2;
            if(x_low < range.first){
                x_low = range.first;
            }
            f_x_low = nll_r(x_low);
            if(f_x_low > 0){
                products_sink->set_product(lower_columns[i], secant(x_low, x_high, x_acurracy, f_x_low, f_x_high, deltanll_levels[i] / 1000, nll_r));
                break;
            }
            else if(f_x_low==0.0 || x_low == range.first){
                products_sink->set_product(lower_columns[i], x_low);
                break;
            }
            else{
                x_high = x_low;
                f_x_high = f_x_low;
            }
        }
        if(k > k_max){
            stringstream ss;
            ss << "Could not find lower value for interval after " << k_max << " iterations (min at: " << value_at_minimum << "; initial step: " << initial_step << "; current step: " << step << ")";
            throw Exception(ss.str());
        }
    }
}

deltanll_intervals::deltanll_intervals(const theta::Configuration & cfg): Producer(cfg),
   pid(cfg.pm->get<VarIdManager>()->get_par_id(cfg.setting["parameter"])), re_minimize(true), start_step_ranges_init(false){
    Setting s = cfg.setting;
    minimizer = theta::PluginManager<Minimizer>::build(theta::Configuration(cfg, s["minimizer"]));
    size_t ic = s["clevels"].size();
    if (ic == 0) {
        throw ConfigurationException("deltanll_intervals: empty clevels.");
    }
    for (size_t i = 0; i < ic; i++) {
        clevels.push_back(s["clevels"][i]);
    }
    if(s.exists("re-minimize")){
        re_minimize = s["re-minimize"];
    }
    deltanll_levels.resize(clevels.size());
    for(size_t i=0; i<clevels.size(); ++i){
        if(clevels[i] < 0.0) throw invalid_argument("clevel < 0 not allowed.");
        if(clevels[i] >= 1.0) throw invalid_argument("clevel >= 1.0 not allowed.");
        deltanll_levels[i] = utils::phi_inverse((1+clevels[i])/2);
        deltanll_levels[i] *= deltanll_levels[i]*0.5;
    }
    declare_products();
}

void deltanll_intervals::declare_products(){
    c_maxl = products_sink->declare_product(*this, "maxl", theta::typeDouble);
    for(size_t i=0; i<clevels.size(); ++i){
        stringstream ss;
        ss << "lower" << setw(5) << setfill('0') << static_cast<int>(clevels[i] * 10000 + 0.5);
        lower_columns.push_back(products_sink->declare_product(*this, ss.str(), theta::typeDouble));
        ss.str("");
        ss << "upper" << setw(5) << setfill('0') << static_cast<int>(clevels[i] * 10000 + 0.5);
        upper_columns.push_back(products_sink->declare_product(*this, ss.str(), theta::typeDouble));
    }
    
}

REGISTER_PLUGIN(deltanll_intervals)
