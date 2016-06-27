#include "plugins/nll_der.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/minimizer.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"
#include "interface/variables-utils.hpp"
#include "interface/model.hpp"

using namespace theta;

void nll_der::produce(const theta::Data & data, const theta::Model & model){
    if(not init){
        ParIds model_pars = model.get_parameters();
        if(not (s_plus_b->get_parameters() == model_pars) or not (b_only->get_parameters() == model_pars)){
            throw std::invalid_argument("parameters in s+b / b only distributions do not coincide with model parameters");
        }
        s_plus_b_width.set(asimov_likelihood_widths(model, s_plus_b));
        b_only_width.set(asimov_likelihood_widths(model, b_only));
        if(!epsilon_abs){
            epsilon_abs = s_plus_b_width.get(parameter) * (*epsilon_rel);
        }
        init = true;
    }
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    nll->set_override_distribution(b_only);
    MinimizationResult minres = minimizer->minimize(*nll, b_only_mode, b_only_width, b_only_support);
    // with this minimum evaluate derivative:
    nll->set_override_distribution(s_plus_b);
    ParValues p1(minres.values);
    const double nll0 = (*nll)(minres.values);
    const double p0p = minres.values.get(parameter);
    const std::pair<double, double> & p_range = s_plus_b_support.get(parameter);
    if(p0p + (*epsilon_abs) >= p_range.first and p0p + (*epsilon_abs) <= p_range.second){
        p1.set(parameter, p0p + (*epsilon_abs));
    }
    else if(p0p - (*epsilon_abs) >= p_range.first and p0p - (*epsilon_abs) <= p_range.second){
        p1.set(parameter, p0p - (*epsilon_abs));
    }
    else{
        throw std::invalid_argument("s plus b support neither contains p0 + epsilon or p0 - epsilon");
    }
    const double nll1 = (*nll)(p1);
    const double h = p1.get_unchecked(parameter) - minres.values.get_unchecked(parameter);
    const double der = (nll1 - nll0) / h;
    products_sink->set_product(c_der, der);
    
    double adev_b = 0, adev_sb = 0;
    const ParIds & parameters = s_plus_b->get_parameters();
    for(ParIds::const_iterator pit=parameters.begin(); pit!=parameters.end(); ++pit){
        adev_sb += fabs(minres.values.get(*pit) - s_plus_b_mode.get(*pit)) / s_plus_b_width.get(*pit);
        adev_b += fabs(minres.values.get(*pit) - b_only_mode.get(*pit)) / s_plus_b_width.get(*pit);
    }
    products_sink->set_product(c_adev_sb, adev_sb);
    products_sink->set_product(c_adev_b, adev_b);
}

nll_der::nll_der(const theta::Configuration & cfg): Producer(cfg), parameter(get_parameter(cfg, "parameter")), epsilon_rel(1e-4), init(false){
    Setting s = cfg.setting;
    if(s.exists("override-parameter-distribution")) throw ConfigurationException("override-parameter-distribution not allowed");
    minimizer = theta::PluginManager<Minimizer>::build(theta::Configuration(cfg, s["minimizer"]));
    s_plus_b = theta::PluginManager<Distribution>::build(theta::Configuration(cfg, s["signal-plus-background-distribution"]));
    b_only = theta::PluginManager<Distribution>::build(theta::Configuration(cfg, s["background-only-distribution"]));
    s_plus_b->mode(s_plus_b_mode);
    s_plus_b_support.set_from(*s_plus_b);
    b_only->mode(b_only_mode);
    b_only_support.set_from(*b_only);
    if(s.exists("epsilon_abs")){
        if(s.exists("epsilon_rel")) throw ConfigurationException("Specifying both, 'epsilon_abs' and 'epsilon_rel' is not allowed");
        epsilon_abs = s["epsilon_abs"];
        if(*epsilon_abs <= 0.0){
            throw ConfigurationException("epsilon_abs <= 0.0 not allowed");
        }
    }
    else if(s.exists("epsilon_rel")){
        epsilon_rel = s["epsilon_rel"];
        if(*epsilon_rel <= 0.0){
            throw ConfigurationException("epsilon_rel <= 0.0 not allowed");
        }
    }
    c_der = products_sink->declare_product(*this, "der", theta::typeDouble);
    c_adev_sb = products_sink->declare_product(*this, "adev_sb", theta::typeDouble);
    c_adev_b = products_sink->declare_product(*this, "adev_b", theta::typeDouble);
}

REGISTER_PLUGIN(nll_der)

