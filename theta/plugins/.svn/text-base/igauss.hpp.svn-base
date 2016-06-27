#ifndef PLUGINS_IGAUSS_HPP
#define PLUGINS_IGAUSS_HPP

#include "interface/distribution.hpp"

/** \brief A product of independent Gaussians, including Gaussians with infinity and zero width
 * 
 * This is equivalent to using a product_distribution with gauss1d, delta_distribution and flat_distribution.
 * It is implemented as a single Distribution plugin for performance reasons only.
 * 
 * Configured via
 * \code
 * {
 *   type = "igauss";
 *   parameters = ("p0", "p1");
 *   mu = (0.0, "mean0");
 *   sigma = (0.0, "inf");
 *   ranges = (("-inf", "inf"), (-2.3, 2.));
 * };
 * \endcode
 * 
 * \c parameters is a list of parameter names; this defines the order in which other lists are interpreted
 * 
 * \c mu, \c sigma and \c ranges are the mean, standard deviation and range for the normal distributions. As entries
 * of \c mu, one can specify either a constant or a parameter name. For \c sigma, values 0.0 and infinity
 * are allowed, corresponding to a delta distribution and flat distribution in this parameter, respectively.
 * 
 * For a width of 0.0, \c mu has to be a constant and the range has to contain only the mean. Conversely, specifying a
 * range with only a single point is only allowed for a width of 0.0.
 * 
 * For a distribution with infinite width, the samples value will always be the value given in \c mu. This corresponds
 * to the behavior if specifying "fix-sample-value" in flat_distribution.
 */
class igauss: public theta::Distribution{
public:
    igauss(const theta::Configuration & cfg);
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual const std::pair<double, double> & support(const theta::ParId & p) const;
    virtual double eval_nl_with_derivative(const theta::ParValues & values, theta::ParValues & derivative) const;
private:
    // build two separate structures: one if mu is a parameter, one if mu is a constant
    struct parset_muconstant{
        theta::ParId parameter;
        double mu, sigma, range_low, range_high;
        parset_muconstant(const theta::ParId & parameter_, double mu_, double sigma_,
                          double range_low_, double range_high_): parameter(parameter_), mu(mu_), sigma(sigma_),
                          range_low(range_low_), range_high(range_high_){}
    };
    
    struct parset_muvar {
        theta::ParId parameter, mu;
        double sigma, range_low, range_high;
        parset_muvar(const theta::ParId & parameter_, const theta::ParId & mu_, double sigma_,
                     double range_low_, double range_high_): parameter(parameter_), mu(mu_), sigma(sigma_),
                     range_low(range_low_), range_high(range_high_){}
    };
    std::vector<parset_muconstant> mconst;
    std::vector<parset_muvar> mvar;
    
    // save the ranges redundantly for support():
    theta::Ranges ranges;
};

#endif
