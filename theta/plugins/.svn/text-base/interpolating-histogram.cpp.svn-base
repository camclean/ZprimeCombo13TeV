//#include "core.hpp"
#include "plugins/interpolating-histogram.hpp"

using namespace std;
using namespace theta;

void interpolating_histo::fill_h(const ParValues & values) const {
    h.set_all_values(1.0);
    const size_t n_sys = hplus.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]);
        const Histogram1D & t_sys = delta > 0 ? hplus[isys] : hminus[isys];
        if (t_sys.get_nbins() == 0)continue;
        h.multiply_with_ratio_exponented(t_sys, h0, fabs(delta));
    }
    h *= h0;
}

void interpolating_histo::add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    h_wu.set(h);
    hres.add_with_coeff(coeff, h_wu);
}

void interpolating_histo::add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    hres.add_with_coeff(coeff, h);
}

void interpolating_histo::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h0.get_nbins();
    xmin = h0.get_xmin();
    xmax = h0.get_xmax();
}

interpolating_histo::interpolating_histo(const Configuration & ctx){
    Setting psetting = ctx.setting["parameters"];
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    //build nominal histogram:
    h0 = get_constant_histogram(Configuration(ctx, ctx.setting["nominal-histogram"])).get_values_histogram();
    size_t n = psetting.size();
    //note: allow n==0 to allow the user to disable systematics.
    // In case of unintentional type error (parameters="delta1,delta2";), user will get a warning about
    // the unused delta*-{plus,minus}-histogram blocks anyway ...
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = vm->get_par_id(par_name);
        par_ids.insert(pid);
        vid.push_back(pid);
        stringstream setting_name;
        //plus:
        setting_name << par_name << "-plus-histogram";
        hplus.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[setting_name.str()])).get_values_histogram());
        //minus:
        setting_name.str("");
        setting_name << par_name << "-minus-histogram";
        hminus.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[setting_name.str()])).get_values_histogram());
    }
    h = h0;
    
    const size_t nsys = hplus.size();
    std::set<ParId> pid_set;
    for(size_t i=0; i<nsys; i++){
        pid_set.insert(vid[i]);
        h0.check_compatibility(hplus[i]);
        h0.check_compatibility(hminus[i]);
    }
    //to make calculation of derivatives easier, we do not allow the same parameter twice.
    if(pid_set.size()!=nsys){
        throw invalid_argument("interpolating_histo::interpolating_histo: having one parameter parametrizing two interpolations is not supported.");
    }
    h = h0;
    h_wu.set(h0);
}

REGISTER_PLUGIN(interpolating_histo)
