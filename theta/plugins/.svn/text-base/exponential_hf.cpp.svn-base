#include "plugins/exponential_hf.hpp"

void exponential_hf::add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const{
    double lambda = lambda0 + c * values.get(parameter);
    coeff *= normalize_to / fill_h(lambda);
    h_wu.set(h);
    hres.add_with_coeff(coeff, h_wu);
}

void exponential_hf::add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const{
    double lambda = lambda0 + c * values.get(parameter);
    coeff *= normalize_to / fill_h(lambda);
    hres.add_with_coeff(coeff, h);
}

void exponential_hf::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h.get_nbins();
    xmin = h.get_xmin();
    xmax = h.get_xmax();
}

// return the sum of bincontents of h, which is filled with exp(lambda*x), integrated in each bin using binborders information.
double exponential_hf::fill_h(double lambda) const{
    double res = 0.0;
    
    // exp_low = value at lower bin border; exp_high = value at upper bin border
    double exp_low = exp(lambda * binborders[0]);
    size_t nbins = binborders.size()-1;
    for(size_t i=0; i<nbins; ++i){
        double exp_high = exp(lambda * binborders[i+1]);
        double bincontent = fabs(exp_high - exp_low);
        res += bincontent;
        h.set(i, bincontent);
        // prepare for next iteration:
        exp_low = exp_high;
    }
    return res;
}


exponential_hf::exponential_hf(const theta::Configuration & cfg): parameter(cfg.pm->get<theta::VarIdManager>()->get_par_id(cfg.setting["parameter"])){
    par_ids.insert(parameter);
    lambda0 = cfg.setting["lambda0"];
    c = cfg.setting["c"];
    normalize_to = cfg.setting["normalize_to"];
    theta::Setting s_binborders = cfg.setting["binborders"];
    size_t n = s_binborders.size();
    if(n<2){
        throw theta::ConfigurationException("binborders < 2 not allowed.");
    }
    binborders.reserve(n);
    for(size_t i=0; i<n; i++){
        binborders.push_back(s_binborders[i]);
        // must be strctly ordered:
        if(i > 0){
            if(binborders[i-1] >= binborders[i]){
                throw theta::ConfigurationException("binborder values must increase monotonically!");
            }
        }
    }
    h.reset(n-1, binborders[0], binborders[n-1]);
    h_wu.set(h);
}

REGISTER_PLUGIN(exponential_hf)
