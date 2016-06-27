#ifndef EXPONENTIAL_HF_HPP
#define EXPONENTIAL_HF_HPP

#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"


/** \brief A HistogramFunction of an exponential function
 *
 * The configuration is done via:
 * \code
 *   histogram = {
 *      type = "exponential_hf";
 *      lambda0 = -0.01;
 *      c = 0.01;
 *      parameter = "delta1";
 *      normalize_to = 5437.3;
 *      binborders = (0.0, 1.0, 2.0, 3.0);
 *   };
 * \endcode
 * 
 * \c lambda0, \c c, and \c parameter control the shape of the exponential function, which is  exp((lambda0 + c * parameter) * x), so the
 *   coefficient for x in the exponential is lambda0 + c * parameter.
 *  
 * Note that this functional form is chosen instead of the simpler but equivalent exp(parameter*x) to allow the parameter to have a prior
 * distribution with mean 0 and width 1, as for other uncertainties, as this allows a more flexible set of correlations introduced for the
 * exponential function.
 * 
 * \c normalize_to controls what the sum of bin contents should be
 * 
 * \c binborders is a list/tuple of nbins+1 numbers that specify the bin border positions along the x axis.
 */
class exponential_hf : public theta::HistogramFunction {
public:
    exponential_hf(const theta::Configuration & ctx);
    
    virtual void add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const;
    virtual void add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const;
    virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const;
    
private:
    double fill_h(double lambda) const;
    
    theta::ParId parameter;
    double lambda0, c, normalize_to;
    std::vector<double> binborders;
    
    mutable theta::Histogram1D h;
    mutable theta::Histogram1DWithUncertainties h_wu;
};

#endif
