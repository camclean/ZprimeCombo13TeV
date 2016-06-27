#include "plugins/linear-histo-morph.hpp"
#include "interface/utils.hpp"

using namespace std;
using namespace theta;

void linear_histo_morph::fill_h(const ParValues & values) const {
    h.set_all_values(1.0);
    const size_t n_sys = kappa_plus.size();
    //1. interpolate linearly in each bin; also calculate normalization
    double scale_unc = 1;
    for (size_t isys = 0; isys < n_sys; ++isys) {
        const double delta = values.get(parameters[isys]);
        if(delta==0.0) continue;
        const Histogram1D & kappa_sys = delta > 0 ? kappa_plus[isys] : kappa_minus[isys];
        if(kappa_sys.get_nbins() > 0)
           h.add_with_coeff(fabs(delta), kappa_sys);
        double relexp = delta > 0 ? plus_relexp[isys] : minus_relexp[isys];
        // multiply with 1 + |\delta| relexp, but cutoff at 0:
        scale_unc *= max(0.0, 1.0 + fabs(delta) * relexp);
    }
    h *= h0;
    //2. lower bin cutoff
    for(size_t i=0; i < h.size(); ++i){
       if(h.get(i) < 0.0){
         h.set(i, 0.0);
         //throw UnphysicalPredictionException();
       }
    }
    //3.a. rescale to nominal:
    h *= h0exp / h.get_sum();
    //3.b. apply scale uncertainty
    h *= scale_unc;
}

void linear_histo_morph::add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    h_wu.set(h);
    hres.add_with_coeff(coeff, h_wu);
}

void linear_histo_morph::add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    hres.add_with_coeff(coeff, h);
}

void linear_histo_morph::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h0.get_nbins();
    xmin = h0.get_xmin();
    xmax = h0.get_xmax();
}

linear_histo_morph::linear_histo_morph(const Configuration & ctx){
    Setting psetting = ctx.setting["parameters"];
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    size_t n = psetting.size();
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = vm->get_par_id(par_name);
        par_ids.insert(pid);
        parameters.push_back(pid);
        if(ctx.setting.exists(par_name + "-kappa-plus-histogram"))
           kappa_plus.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[par_name + "-kappa-plus-histogram"])).get_values_histogram());
        else
           kappa_plus.push_back(Histogram1D());
        if(ctx.setting.exists(par_name + "-kappa-minus-histogram"))
            kappa_minus.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[par_name + "-kappa-minus-histogram"])).get_values_histogram());
        else
           kappa_minus.push_back(Histogram1D());
        if(ctx.setting.exists(par_name + "-plus-relexp"))
            plus_relexp.push_back(ctx.setting[par_name + "-plus-relexp"]);
        else
            plus_relexp.push_back(0.0);
        if(ctx.setting.exists(par_name + "-minus-relexp"))
            minus_relexp.push_back(ctx.setting[par_name + "-minus-relexp"]);
        else
            minus_relexp.push_back(0.0);
    }    
    const size_t nsys = kappa_plus.size();
    std::set<ParId> pid_set;
    h0 = get_constant_histogram(Configuration(ctx, ctx.setting["nominal-histogram"])).get_values_histogram();
    h = h0;
    for(size_t i=0; i < nsys; i++){
        pid_set.insert(parameters[i]);
        if(kappa_plus[i].get_nbins() > 0)
           h0.check_compatibility(kappa_plus[i]);
        if(kappa_minus[i].get_nbins() > 0)
           h0.check_compatibility(kappa_minus[i]);
    }
    if(pid_set.size()!=nsys){
        throw invalid_argument("linear_histo_morph: duplicate parameter in parameter list.");
    }
    h0exp = ctx.setting["nominal-expectation"];
    h0 *= h0exp / h0.get_sum();
    h = h0;
    h_wu.set(h);
}

REGISTER_PLUGIN(linear_histo_morph)
