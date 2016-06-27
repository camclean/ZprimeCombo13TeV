#ifndef PLUGINS_GAMMA_DISTRIBUTION_HPP
#define PLUGINS_GAMMA_DISTRIBUTION_HPP

#include "interface/distribution.hpp"
#include "interface/exception.hpp"
#include <vector>
#include <map>


/** \brief gamma distribution
 *
 * The gamma distribution can be used to model an uncertainty without the need of explicit
 * truncation at zero as required e.g. for gauss distribution. (Another common possibility is to use
 * a lognormal distribution.)
 *
 * The gamma distribution has pdf
 * \f[
 * x^{k-1} \frac{\exp(-x / \theta)}{\Gamma(k)\theta^k}
 * \f]
 * with mean \f$ \mu = k \theta \f$ and variance \f$ \sigma^2 = k \theta^2 \f$.
 * 
 * The gamma distribution is configured via the parameters
 * 'mu' and 'sigma'. To model a relative uncertanty for some process just set the 'mu' to 1.0
 * and the 'sigma' to the relative uncertainty and use the parameter as a factor
 * in the corresponding coefficient-function.
 * 
 *
 * Configuration is done with a setting group like
 * \code
 * dist = {
 *   type = "gamma_distribution";
 *   parameter = "factor1";
 *   mu = 1.0;
 *   sigma = 0.2;
 *   range = (0.0, 5.0); // optional, default is (0.0, "inf")
 * };
 * \endcode
 * 
 * Note that all parameters are doubles and must be &gt; 0.
 */
class gamma_distribution: public theta::Distribution{
public:
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual const std::pair<double, double> & support(const theta::ParId & p) const;
    virtual double width(const theta::ParId & p) const;
    gamma_distribution(const theta::Configuration & cfg);
private:
    std::pair<double, double> supp;
    double k, theta, eoverkp1;
};

#endif

