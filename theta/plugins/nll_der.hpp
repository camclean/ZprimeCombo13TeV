#ifndef PLUGIN_NLL_DER_HPP
#define PLUGIN_NLL_DER_HPP

#include "interface/decls.hpp"
#include "interface/producer.hpp"
#include "interface/variables.hpp"
#include "interface/distribution.hpp"

#include <boost/optional.hpp>

/** \brief A producer to create negative log-likelihood derivatives
 *
 * This can be used as an alternative test statistic definition
 *
 *
 * Configuration is done via a setting group like
 * \code
 * hypotest = {
 *   type = "nll_der";
 *   name = "nllder";
 *   minimizer = "@myminuit";
 *   background-only-distribution = {...};
 *   signal-plus-background-distribution = {...};
 *   parameter = "beta_signal";
 *   epsilon_rel = 1e-4; // optional; default is 1e-4
 *   epsilon_abs = 1e-4; // optional; default is None
 * };
 * 
 * myminuit = {...}; // minimizer definition
 * \endcode
 *
 * \c type is always "nll_der" to select this producer.
 *
 * \c name is a unique producer name of your choice; it is used to construct column names in the output database. It may only contain alphanumeric
 *    characters (no spaces, special characters, etc.).
 *
 * \c minimizer is the configuration path to a \link theta::Minimizer minimizer \endlink definition to be used
 *    for minimization of the negative-log-likelihood.
 *
 * \c background-only-distribution defines the Distribution to use to find the minimum of the negative log-likelihood and \c signal-plus-background-distribution
 *  is the Distribution corresponding to the distribution used for evaluating the likelihood derivative as well as the difference from the Asimov case.
 * 
 * \c parameter is the parameter the derivative should be calculated in
 * 
 * \c epsilon_rel and \c epsilon_abs control the used step size: epsilon_rel is the relative step size in terms of the asimov likelihood width in this parameter
 *   The default of 1e-4 usually provides robust results. \c epsilon_abs can be used to specify the step size directly. Specifying both epsilons is not allowed.
 * 
 * The usual use case is using "der" as test statistic with \c background-only-distribution being the null hypothesis, fixing \c parameter ;
 * \c signal-plus-background-distribution is the (composite) alternative in which \c parameter must not be fixed to allow calculating the derivative.
 * If both p0 + epsilon or p0 - epsilon (where p0 is the estimate from \c background-only-distribution) must be within the parameter bounds of \c signal-plus-background-distribution ;
 * otherwise, a logic failure is reported and theta aborts.
 * 
 * The created columns are:
 *   * "der" contining the nll derivative w.r.t. parameter
 *   * "adev_sb", "adev_b" ("adev" = "Asimov deviation"): the deviation of the parameters at the background-only fit w.r.t. the s+b model and w.r.t. the b-only model
 * 
 * The deviation is defined as sum_i |p_i - p_0i| / sigma_p0i
 * where i runs over all parameters, p_i is the mle of p, p_0i is the Asimov value for parameter i, and sigma_p0i is the Asimov width for parameter i. The Asimov
 * values and widths are defined by the \c signal-plus-background-distribution for "adev_sb" and by \c background-only-distribution for "adev_b".
 */
class nll_der: public theta::Producer{
public:
    nll_der(const theta::Configuration & cfg);
    virtual void produce(const theta::Data &, const theta::Model&);
private:
    std::auto_ptr<theta::Minimizer> minimizer;
    boost::shared_ptr<theta::Distribution> s_plus_b;
    boost::shared_ptr<theta::Distribution> b_only;
    theta::ParId parameter;
    boost::optional<double> epsilon_rel;
    boost::optional<double> epsilon_abs;
    
    bool init;
    theta::ParValues s_plus_b_mode, b_only_mode;
    theta::ParValues s_plus_b_width, b_only_width;
    theta::Ranges s_plus_b_support, b_only_support;
    
    theta::Column c_der, c_adev_sb, c_adev_b;
};

#endif
