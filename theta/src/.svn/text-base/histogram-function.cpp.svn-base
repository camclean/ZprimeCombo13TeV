#include "interface/histogram-function.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::HistogramFunction);

using namespace theta;

void HistogramFunction::eval_and_add_derivatives(Histogram1D & result, std::map<ParId, Histogram1D> & derivatives, double coeff, const ParValues & values) const{
    throw Exception("not implemented for " + demangle(typeid(this).name()));
}

void HistogramFunction::eval_and_add_derivatives(Histogram1DWithUncertainties & result, std::map<ParId, Histogram1DWithUncertainties> & derivatives, double coeff, const ParValues & values) const{
    throw Exception("not implemented for " + demangle(typeid(this).name()));
}

void ConstantHistogramFunction::add_with_coeff_to(Histogram1D & hres, double coeff, const ParValues & values) const{
    hres.add_with_coeff(coeff, h);
}

void ConstantHistogramFunction::add_with_coeff_to(Histogram1DWithUncertainties & hres, double coeff, const ParValues & values) const{
    hres.add_with_coeff(coeff, h_wu);
}

void ConstantHistogramFunction::eval_and_add_derivatives(Histogram1D & result, std::map<ParId, Histogram1D> & derivatives, double coeff, const ParValues & values) const{
    result = h;
    // derivative is 0.0, so don't add anything ...
}


void ConstantHistogramFunction::eval_and_add_derivatives(Histogram1DWithUncertainties & result, std::map<ParId, Histogram1DWithUncertainties> & derivatives, double coeff, const ParValues & values) const{
    result = h_wu;
    // derivative is 0.0, so don't add anything ...
}

        
void ConstantHistogramFunction::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h.get_nbins();
    xmin = h.get_xmin();
    xmax = h.get_xmax();
}

void ConstantHistogramFunction::set_histo(const Histogram1DWithUncertainties & h_){
    h = h_.get_values_histogram();
    h_wu = h_;
}

ConstantHistogramFunction::ConstantHistogramFunction(){}

Histogram1DWithUncertainties theta::get_constant_histogram(const Configuration & cfg){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(cfg);
    if(hf->get_parameters().size()!=0){
        throw std::invalid_argument("Histogram defined in path '" + cfg.setting.get_path() + "' is not constant");
    }
    size_t nbins;
    double xmin, xmax;
    hf->get_histogram_dimensions(nbins, xmin, xmax);
    Histogram1DWithUncertainties res(nbins, xmin, xmax);
    hf->add_with_coeff_to(res, 1.0, ParValues());
    return res;
}
