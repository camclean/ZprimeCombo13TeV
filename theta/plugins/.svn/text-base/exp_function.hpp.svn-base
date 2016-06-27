#ifndef PLUGIN_EXP
#define PLUGIN_EXP

#include "interface/phys.hpp"

/** \brief Function returning exp(lambda * p) where p is a configurable parameter and lambda a literal constant
 *
 * \code
 *  func = {
 *     type = "exp_function";
 *     parameter = "par0";
 *     lambda = 0.01;
 *     // alternatively:
 *     lambda_minus = 0.01;
 *     lambda_plus = 0.02;
 *     // to multiply several exponentials, use
 *     parameters = ("par0", "par1", "par2");
 *     lambdas_minus = (0.02, 0.03, 0.01);
 *     lambdas_plus = (0.02, 0.03, 0.01);
 *  };
 * \endcode
 *
 * This function is typically used to implement log-normal rate uncertainties by using a nuisanace parameter with a Gaussian distribution around 0 with width 1.
 * In this case, \c lambda is the rate uncertainty. In case of assymetric rate uncertainties, you can supply two lambda values via the "lambdas" setting: the first one
 * is used if the parameter value is below zero, the other one if it is above zero. (Note that "lambda = 0.1;" is completely equivalent to "lambda_minus = 0.1; lambda_plus = 0.1;".
 * In particular, there is no additional sign).
 */
class exp_function: public theta::Function{
private:
    std::vector<theta::ParId> v_pids_pm;
    std::vector<double> lambdas_plusminus;
    
    // symmetric:
    std::vector<theta::ParId> v_pids_s;
    std::vector<double> lambdas_s;

public:
    /// constructor for the plugin system
    exp_function(const theta::Configuration & cfg);
    /// overloaded evaluation operator of theta::Function
    virtual double operator()(const theta::ParValues & v) const;
    
    virtual double eval_with_derivative(const theta::ParValues & v, theta::ParValues & der) const;
};


#endif
