#include "plugins/simple_linear_histomorph.hpp"

using namespace std;
using namespace theta;

void simple_linear_histomorph::fill_h(const ParValues & values) const {
    h = h0;
    const size_t n_sys = hplus_diff.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]);
        if(delta==0.0) continue;
        const Histogram1D & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
        if (t_sys.get_nbins() == 0)continue;
        h.add_with_coeff(fabs(delta), t_sys);
    }
    for(size_t i=0; i<h.get_nbins(); ++i){
       h.set(i, max(h.get(i), 0.0));
    }
}

void simple_linear_histomorph::add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    h_wu.set(h);
    hres.add_with_coeff(coeff, h_wu);
}

void simple_linear_histomorph::add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    hres.add_with_coeff(coeff, h);
}


void simple_linear_histomorph::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h0.get_nbins();
    xmin = h0.get_xmin();
    xmax = h0.get_xmax();
}


simple_linear_histomorph::simple_linear_histomorph(const Configuration & cfg){
    Setting psetting = cfg.setting["parameters"];
    //build nominal histogram:
    h0 = get_constant_histogram(Configuration(cfg, cfg.setting["nominal-histogram"])).get_values_histogram();
    size_t n = psetting.size();
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = vm->get_par_id(par_name);
        par_ids.insert(pid);
        vid.push_back(pid);
        string setting_name;
        //plus:
        setting_name = par_name + "-plus-histogram";
        if(cfg.setting.exists(setting_name)){
           hplus_diff.push_back(get_constant_histogram(Configuration(cfg, cfg.setting[setting_name])).get_values_histogram());
           hplus_diff.back().check_compatibility(h0);
           hplus_diff.back().add_with_coeff(-1.0, h0);
        }
        else{
           hplus_diff.push_back(Histogram1D());
        }
        //minus:
        setting_name = par_name + "-minus-histogram";
        if(cfg.setting.exists(setting_name)){
           hminus_diff.push_back(get_constant_histogram(Configuration(cfg, cfg.setting[setting_name])).get_values_histogram());
           hminus_diff.back().check_compatibility(h0);
           hminus_diff.back().add_with_coeff(-1.0, h0);
        }
        else{
           hminus_diff.push_back(Histogram1D());
        }
    }
    h = h0;
    h_wu.set(h);
}

REGISTER_PLUGIN(simple_linear_histomorph)
