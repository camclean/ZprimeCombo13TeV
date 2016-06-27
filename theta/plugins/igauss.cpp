#include "plugins/igauss.hpp"
#include "interface/plugin.hpp"
#include "interface/random.hpp"

#include <boost/foreach.hpp>

using namespace theta;

igauss::igauss(const Configuration & cfg){
    Setting s = cfg.setting;
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    const size_t n = s["parameters"].size();
    if(n != s["mu"].size() or n != s["sigma"].size() or n != s["ranges"].size()){
        throw ConfigurationException("lists 'parameters', 'mu', 'sigma', 'ranges' have not the same length");
    }
    for(size_t i=0; i<n; ++i){
        ParId pid = vm->get_par_id(s["parameters"][i]);
        par_ids.insert(pid);
        double sigma = s["sigma"][i].get_double_or_inf();
        double range_low = s["ranges"][i][0].get_double_or_inf();
        double range_high = s["ranges"][i][1].get_double_or_inf();
        // range must be valid (low <= high); the case low == high must happen iff sigma == 0.0:
        if(range_low > range_high || (range_low==range_high && sigma > 0.0) || (range_low < range_high && sigma <= 0.0)){
            throw ConfigurationException("mu and range setting invalid");
        }
        ranges.set(pid, std::make_pair(range_low, range_high));
        if(s["mu"][i].get_type() == Setting::TypeFloat){
            double mu = s["mu"][i]; // no inf here
            mconst.push_back(parset_muconstant(pid, mu, sigma, range_low, range_high));
        }
        else{
            ParId mu = vm->get_par_id(s["mu"][i]);
            distribution_par_ids.insert(mu);
            mvar.push_back(parset_muvar(pid, mu, sigma, range_low, range_high));
        }
    }
}


void igauss::sample(theta::ParValues & result, theta::Random & rnd) const{
    BOOST_FOREACH(const parset_muconstant & ps, mconst){
        if(ps.sigma == 0.0 or std::isinf(ps.sigma)){
            result.set(ps.parameter, ps.mu);
        }
        else{
            double val;
            do{
                val = rnd.gauss(ps.sigma) + ps.mu;
            }while(val < ps.range_low || val > ps.range_high);
            result.set(ps.parameter, val);
        }
    }
    BOOST_FOREACH(const parset_muvar & ps, mvar){
        const double mu = result.get_unchecked(ps.mu);
        if(ps.sigma == 0.0 or std::isinf(ps.sigma)){
            result.set(ps.parameter, mu);
        }
        else{
            double val;
            do{
                val = rnd.gauss(ps.sigma) + mu;
            }while(val < ps.range_low || val > ps.range_high);
            result.set(ps.parameter, val);
        }
    }
}

void igauss::mode(ParValues & result) const{
    BOOST_FOREACH(const parset_muconstant & ps, mconst){
        result.set(ps.parameter, ps.mu);
    }
    BOOST_FOREACH(const parset_muvar & ps, mvar){
        result.set(ps.parameter, result.get(ps.mu));
    }
}

double igauss::eval_nl(const ParValues & values) const{
    double result = 0;
    BOOST_FOREACH(const parset_muconstant & ps, mconst){
        const double val = values.get_unchecked(ps.parameter);
        if(val < ps.range_low || val > ps.range_high) return std::numeric_limits<double>::infinity();
        if(ps.sigma == 0.0 or std::isinf(ps.sigma)) continue;
        const double d = (val - ps.mu) / ps.sigma;
        result += 0.5 * d * d;
    }
    BOOST_FOREACH(const parset_muvar & ps, mvar){
        const double val = values.get_unchecked(ps.parameter);
        if(val < ps.range_low || val > ps.range_high) return std::numeric_limits<double>::infinity();
        if(ps.sigma == 0.0 or std::isinf(ps.sigma)) continue;
        const double d = (val - values.get_unchecked(ps.mu)) / ps.sigma;
        result += 0.5 * d * d;
    }
    return result;
}

double igauss::eval_nl_with_derivative(const theta::ParValues & values, theta::ParValues & derivative) const{
    double result = 0;
    BOOST_FOREACH(const parset_muconstant & ps, mconst){
        const double val = values.get_unchecked(ps.parameter);
        if(val < ps.range_low || val > ps.range_high){
            result = std::numeric_limits<double>::infinity();
            derivative.set(ps.parameter, 0.0);
        }
        else if(ps.sigma == 0.0 or std::isinf(ps.sigma)){
            derivative.set(ps.parameter, 0.0);
        }
        else{
            const double d = (val - ps.mu) / ps.sigma;
            derivative.set(ps.parameter, d / ps.sigma);
            result += 0.5 * d * d;
        }
    }
    BOOST_FOREACH(const parset_muvar & ps, mvar){
        const double val = values.get_unchecked(ps.parameter);
        if(val < ps.range_low || val > ps.range_high){
            result = std::numeric_limits<double>::infinity();
            derivative.set(ps.parameter, 0.0);
        }
        else if(ps.sigma == 0.0 or std::isinf(ps.sigma)){
            derivative.set(ps.parameter, 0.0);
        }
        else{
            const double d = (val - values.get_unchecked(ps.mu)) / ps.sigma;
            derivative.set(ps.parameter, d / ps.sigma);
            result += 0.5 * d * d;
        }
    }
    return result;
}

const std::pair<double, double> & igauss::support(const ParId & p) const{
    return ranges.get(p);
}

REGISTER_PLUGIN(igauss)
