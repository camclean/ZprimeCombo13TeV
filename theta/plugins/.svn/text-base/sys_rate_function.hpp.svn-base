#ifndef PLUGIN_SYS_RATE_FUNCTION
#define PLUGIN_SYS_RATE_FUNCTION

#include "interface/phys.hpp"

/** \brief Function used for rate uncertainty
 *
 * Configured via a setting like
 * \code
 * coefficient-function = {
 *  type = "sys_rate_function";
 *  factors = ("beta_1", "lumi");
 *  sys_rates = (("delta1", -0.02, 0.02), ("delta2", -0.1, 0.08));
 * };
 * \endcode
 *
 * \c factors is a list of parameters which will be multiplied
 *
 * \c sys_rates is a list of triples with parameter name used to describe the rate uncertainty, the
 *  relative rate uncertainty for the negative parameter values and the relative rate uncertainty for the positive parameter values.
 *
 * Naming the parameters given in \c factors f_1, ... f_n and the \c sys_rates s_1 to s_m with the relative rates
 * given r_{i,-} and r_{i,+} with i=1,...,m, this function returns f_1 * f_2 * ... * f_n * (1 + fabs(s_1) * r_{1,sgn(s_1)})
 * *...* (1 + fabs(s_m) * r_{m,sgn(s_m)}).
 *
 * Note that the function never returns values &lt; 0.0 in order to prevent unphysical model predictions.
 *
 * Usually, the parameters used in \c sys_rates have a (gaussian) prior with mean 0 and width 1. Then, the relative
 * rate changes given in the configuration are the +-1sigma deviations for that systematic uncertainty.
 */
class sys_rate_function: public theta::Function{
private:
    std::vector<theta::ParId> f_pids;
    std::vector<theta::ParId> s_pids;
    std::vector<double> r_plus;
    std::vector<double> r_minus;

public:
    /// constructor used by the plugin system
    sys_rate_function(const theta::Configuration & cfg);
    /// overloaded evaluation operator from theta::Function
    virtual double operator()(const theta::ParValues & v) const;
};

#endif
