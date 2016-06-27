#ifndef SIMPLE_LINEAR_HISTOMORPH_HPP
#define SIMPLE_LINEAR_HISTOMORPH_HPP

#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"


/** \brief A HistogramFunction which interpolates a "zero" Histogram and several "distorted" Histograms with a simple linear ansatz.
 * 
 * Consider using cubiclinear_histomorph instead which is more flexible.
 *
 * The configuration is very similar to \c interpolating_histo :
 * \code
 *   histogram = {
 *      type = "simple_linear_histomorph";
 *      parameters = ("delta1", "delta2");
 *      nominal-histogram = { / * fixed-histogram-specification * / };
 *      delta1-plus-histogram = { / * fixed-histogram-specification * / };
 *      delta1-minus-histogram = { / * fixed-histogram-specification * / };
 *      delta2-plus-histogram = { / * fixed-histogram-specification * / };
 *      delta2-minus-histogram = { / * fixed-histogram-specification * / };
 *   };
 * \endcode
 * Here, <tt>fixed-histogram-specification</tt> is a Histogram Setting block that returns a Histogram
 * which does not depend on any parameters.
 *
 * No random fluctuations are done: the method \c getRandomFluctuation is not overriden from \link theta::HistogramFunction HistogramFunction \endlink
 * and will therefore return the same values as the usual, non-random, evaluation operator, \c operator().
 *
 * The bin content of the returned histogram is calculated independently for each bin: the bin value of bin k is calculated by adding
 * the content of the nominal histogram with fabs(p[i]) * (hplus[i][k] - h0[k]) or with fabs(p[i]) * (hminus[i][k] - h0[k]). In this formula,
 * the parameters are denoted p[i], the corresponding plus and minus histograms are hplus[i] and hminus[i];
 * where hplus[i][k] is bin k of hplus[i]. Which formula is used depends on the sign of p[i].
 *
 * If for a parameter, no plus or minus histogram is given, the nominal value will be used, i.e., giving no histogram has the same effect
 * as specifying the nominal one.
 *
 * The linear interpolation can in principle result in negative bin entries, i.e. in negative predictions. This does not only make no sense from
 * the physics point of view, but is also problematic for the calculation of the value of the likelihood function. Therefore, the bin content
 * is cut off at zero. Even that can be problematic, but it is not addressed here.
 */
class simple_linear_histomorph : public theta::HistogramFunction {
public:
    
    /** \brief Constructor used by the plugin system to build an instance from settings in a configuration file
     */
    simple_linear_histomorph(const theta::Configuration & ctx);
    
    virtual void add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const;
    virtual void add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const;
    virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const;

private:
    void fill_h(const theta::ParValues & values) const;
    
    theta::Histogram1D h0;
    std::vector<theta::Histogram1D> hplus_diff; // hplus_diff[i] + h0 yields hplus
    std::vector<theta::Histogram1D> hminus_diff;
    //the interpolation parameters used to interpolate between hplus and hminus.
    std::vector<theta::ParId> vid;
    //the Histogram returned by operator(). Defined as mutable to allow operator() to be const.
    mutable theta::Histogram1D h;
    mutable theta::Histogram1DWithUncertainties h_wu;
};

#endif
