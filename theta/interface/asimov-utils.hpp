#ifndef ASIMOV_UTILS_HPP
#define ASIMOV_UTILS_HPP

#include "interface/decls.hpp"
#include "interface/data.hpp"
#include "interface/phys.hpp"

#include <boost/shared_ptr.hpp>


namespace theta{

/** \brief The negative log-likelihood function for Asimov data
 */
class asimov_data_nll: public Function{
public:
    
    /// Construct from the model and the override parameter distribution, which will be part of the likelihood. values optionally specify alternative values for some parameters.
    asimov_data_nll(const theta::Model & model, const boost::shared_ptr<Distribution> & override_parameter_distribution, const ParValues & values = ParValues());
    
    /// Evaluate the negative log-likelihood.
    double operator()(const ParValues & values) const;
    
private:
    Data asimov_data;
    std::auto_ptr<NLLikelihood> nll;
};


/** \brief Calculate parameter widths from the asimov data for the given model
 * 
 * The main purpose of this function is to get good initial step sizes for minimization and MCMC integration.
 *
 * The model's distribution will be used to determine the most probable values / modes of all parameters.
 * The prediction using these parameter values (without Poisson smearing = Asimov data) is used to construct
 * a nllikelihood function. From this likelihood function, approximate 1sigma parameter widths are extracted by scanning
 * the function in all directions until the function value reaches (nll value at minimum)+1/2. No minimisation is performed
 * in this process.
 * 
 * If the "+1/2 point" cannot be reached within the parameter boundaries given by the Distribution support,
 * the width for that parameter is set to the total width of the parameter's support. This can be 0, e.g., for
 * delta functions. It can also be infinity, e.g., for parameters with a flat distribution which do not really affect the model prediction.
 * Callers of this method should make sure to treat such cases properly. The zero case is often valid (=parameter is fixed), whereas
 * the large value / infinity case usually indicates a problem bacause it means the likelihood function (and thus the model prediction)
 * hardly depends on this parameter.
 *
 * If override_parameter_distribution is given, it will be used instead of the model's Distribution for both the Asimov
 * data construction as well as for the nllikelihood definition.

 */
ParValues asimov_likelihood_widths(const Model & model, const boost::shared_ptr<Distribution> & override_parameter_distribution);

/** \brief Calculate an approximation of the inverse hessian from asimov data of the given model
 * 
 * Matrix version of asimov_likelihood_widths.
 *
 * The function uses numerical second derivatives of the asimov likelihood to calculate the Hessian, and the invserse
 * Hessian is returned. This is a reasonable estimate for the covariance matrix for many purposes.
 *
 * This method performs a "bootstrap" in the sense that first, characteristic lengths for all axes are found
 * via asimov_likelihood_widths. The step size used in the numerical derivative is then this width times epsilon.
 * A good default value for epsilon seems to be around 1e-4.
 */
Matrix asimov_likelihood_matrix(const Model & model, const boost::shared_ptr<Distribution> & override_parameter_distribution, double epsilon);

}//namespace theta

#endif
