#ifndef PLUGINS_NEWTON_HPP
#define PLUGINS_NEWTON_HPP

#include "interface/minimizer.hpp"
#include "plugins/newton-utils.hpp"
#include <vector>
#include <map>

namespace newton_internal{
    class Linesearch;
}


/** \brief Quasi-Newton Minimizer
 * 
 * Minimizer implementing the quasi-Newton method; using a derivative-free interval line search and the symmetric-rank 1
 * update for the covariance matrix.
 * 
 * Configured via a setting group like
 * \code
 * minimizer = {
 *    type = "newton_minimizer";
 * 
 *    // covergence criterion:
 *    par_eps = 1e-6; // optional; default is 1e-4
 *    maxit = 100000; // optional; default is 10,000
 *    improve_cov = true; // optional; default is false
 *    force_cov_positive = true; // optional, default is false
 *    step_cov = 0.01; // optional; default is 0.1
 * };
 * \endcode
 * 
 * Any quasi-Newton method performs minimization in several steps:
 * 1. At the current point x0, using the gradient g and current estimate of the Hessian, propose a search vector v (direction and magnitude). If the
 *    function is actually quadratic with the current Hessian estimate and the gradient is accurate, x0 + v is the actual minimum.
 * 2. using the direction from 1., attempt minimizing the function along this direction [note that this search does not have to be very accurate.];
 *    this leads to an improved estimate of the minimum
 * 3. At the new minimum estimate, calculate the gradient and use it to update the covariance estimate; start over at 1. at the new point.
 * 
 * The ingredients to specify the algorithm are:
 *  * The linesearch routine (which is a 1D minimization). Currently, the only supported type is "brent", which uses Brent's algorithm to find the minimum
 *    along the line with an accuracy better than eps * (line search step).
 *  * The algorithm used to update the estimate of the inverse Hesse matrix. Here, SR1 is used; see http://en.wikipedia.org/wiki/SR1_formula
 *  * The stopping criterion. Here iteration is stopped either after maxit iterations (in this case, a failure is reported), or if the step from the old to the new point
 *    is "small" in the sense that it is smaller than par_eps / step in each parameter, where step is the initial step size for this parameter.
 * 
 * Given that, most parameters have a straight-forward meaning:
 * 
 * \c par_eps controls the stopping criterion: once the step size in all parameters is smaller than this value (compared to this parameter's step size), the iteration stops.
 * 
 * \c maxit is the maximum number of iterations of the algorithm. One iteration executes all steps 1. to 3. as explained above.
 * 
 * \c improve_cov : if true, an additional step is run after the the minimization to calculate the numerical Hesse matrix from the
 *  numerical derivatives. For each parameter, the gradient is evaluated at the minimum plus the \c step_cov times the step size for this parameter.
 *  These points are are used to build the covariance.
 * 
 * \c force_cov_positive is only relevant in case of \c improve_cov set to \c true. If set to \c true and the Hesse matrix is not positive definite, alpha * identity matrix
 * is added to the Hesse matrix to make the Hesse matrix positive definite.  If \c force_cov_positive is \c false and the Hesse matrix is not positive definite,
 * the minimization fails with a MinimizationException.
 * 
 * \c step_cov is only relevant if \c improve_cov is true. In this case, it controls the step size factor to apply to the original step size for
 * evaluating the gradients for the covariance: Each parameter is varied from the minimum by \c step_cov times the initial step size. The
 * gradient at this point is used to calculate the second derivative of the negative log-likelihood function, the Hesse matrix.
 * 
 * Note that the algorithm typically makes n function evaluations at each point to calculate the gradient, where n is the number of parameters of f.
 * For the line search algorithm, another m ~ O(10) function evaluations are required. maxit controls the number of complete iterations of the
 * algorithm, leading to a maximum number of function evaluations of about maxit * (n + m). In case \c improve_cov is true, additional
 * n*(n+1) function evaluations are required.
 */
class newton_minimizer: public theta::Minimizer{
public:
    struct options{
        int maxit;
        double par_eps;
        int debug; // 0 = print nothing; 1 = details of linesearch; 2 = current x and g; 4 = current h (implies 2)   or "or" of those. Any value > 0 will print some minimal info.
        bool improve_cov, force_cov_positive;
        double step_cov;
        bool second_pass;
        options(): maxit(10000), par_eps(1e-4), debug(0), improve_cov(false), force_cov_positive(false), step_cov(0.1), second_pass(false){}
    };

    newton_minimizer(const options & opts_);
    newton_minimizer(const theta::Configuration & cfg);
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                                               const theta::ParValues & step, const theta::Ranges & ranges);
    virtual theta::MinimizationResult minimize2(const theta::Function & f, const theta::FunctionInfo & info, const theta::ParValues & fixed_parameters);
    virtual boost::shared_ptr<theta::FunctionInfo> create_nll_function_info(const theta::Model & m, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
                                                                            const theta::ParValues & fixed_parameters = theta::ParValues());
private:
    std::auto_ptr<newton_internal::Linesearch> ls;
    options opts;
    bool use_nll_der;
};

#endif
