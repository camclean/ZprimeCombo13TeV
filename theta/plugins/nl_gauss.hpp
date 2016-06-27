#ifndef PLUGIN_NL_GAUSS
#define PLUGIN_NL_GAUSS

#include "interface/phys.hpp"
#include "interface/matrix.hpp"

/** \brief Negative Logarithm of a multivariate Gauss
 *
 * Configure via a setting group like
 * \code
 *   {
 *     type = "nl_gauss";
 *     rows = ("delta1", "delta2", "@func1");
 *     mu = [1.0, 2.0, 3.0];
 *     covariance = ([2.0, 0.0, 0.2], 
 *                   [0.0, 1.0, 0.0],
 *                   [0.2, 0.0, 2.0]);
 *   };
 * \endcode
 * 
 * \c rows is a list which either contains parameter names or Function specifications. The order given here defines
 *    how the order in the other parameters, \c mu and \c covariance, are interpreted
 * 
 * \c mu is the mean value of the Gaussian
 * 
 * \c covariance is the covariance matrix. It must be symmetric and positive definite.
 * 
 * This function is useful as an additional term to the likelihood function, e.g., as a prior for model parameters or if modeling
 * approximate, external likelihoods.
 */
class nl_gauss: public theta::Function{
private:
    boost::ptr_vector<theta::Function> rows;
    std::vector<double> mu;
    theta::Matrix inv_cov;
    
    mutable std::vector<double> row_values;
public:
    /// constructor used by the plugin system
    nl_gauss(const theta::Configuration & cfg);
    /// overloaded evaluation of the function
    virtual double operator()(const theta::ParValues & v) const;
};

#endif
