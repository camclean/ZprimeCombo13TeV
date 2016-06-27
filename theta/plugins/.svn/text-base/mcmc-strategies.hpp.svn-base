#ifndef PLUGINS_MCMC_STRATEGIES_HPP
#define PLUGINS_MCMC_STRATEGIES_HPP

#include "interface/decls.hpp"
#include "interface/mcmc.hpp"

/** \brief MCMC strategy using the full covariance matrix from the numerical derivative of asimov data
 * 
 * Configured via a setting group like
 * \code
 *     mcmc_strategy = {
 *         // common settings:
 *         type = "asimov_der_cov";
 *         name = "mcs";
 *         iterations = 10000;
 *         burn-in = 100; // optional
 *         factor = 1.5; // optional, default is 1.0
 *         // specific to this type:
 *         epsilon = 1e-3; // optional, default is 1e-4
 *     };
 * \endcode
 * 
 * \c epsilon is the step size to use for the numerical derivative, in units of the asimov width in the parameter direction. The asimov width
 *    is defined as the step needed to increase the nll by 1 unit from the minimum.
 * 
 * \c factor is a factor for the step size. The default corresponding to factor = 1.0, is to use 2.38 / sqrt(n) times the covariance.
 */
class asimov_der_cov: public theta::MCMCStrategy{
public:
    explicit asimov_der_cov(const theta::Configuration & cfg);
    virtual void init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution);
    virtual void run_mcmc(const theta::Function & nllikelihood, theta::MCMCResult &res) const;
private:
    theta::Matrix sqrt_cov;
    double epsilon;
};

/** \brief MCMC strategy using a gaussian jumping kernel with diagonal covariance matrix derived from the asimov likelihood
 * 
 * It uses a Gaussian jump kernel with a width 2.38 / sqrt(n) * cov, where cov is the diagonal covariance
 * matrix in which each dimension, the width is given by the asimov likelihood width.
 *
 * Configured via a setting group like
 * \code
 *     mcmc_strategy = {
 *         // common settings:
 *         type = "asimov_widths";
 *         name = "mcs";
 *         iterations = 10000;
 *         burn-in = 100; // optional
 *         factor = 1.5; // optional, default is 1.0
 *     };
 * \endcode
 * 
 * See MCMCStrategy for an explanation of those options.
 * 
 */
class asimov_widths: public theta::MCMCStrategy{
public:
    explicit asimov_widths(const theta::Configuration & cfg);
    virtual void init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution);
    virtual void run_mcmc(const theta::Function & nllikelihood, theta::MCMCResult &res) const;
private:
    std::vector<double> widths;
};

/** \brief MCMC strategy going along one random parameter axis at each step.
 *
 * Configured within a mcmc producer via
 * \code
 *   mcmc_strategy = {
 *       type = "asimov_widths_1d";
 *       name = "mcs";
 *       iterations = 10000;
 *       burn-in = 100; // optional
 *       factor = 2.0; // optional, default is 1.0
 *  };
 * \endcode
 * 
 * The proposed mcmc jump is a Gaussian along a random axis, with a width \c factor times the asimov width.
 * 
 * The default of \c factor = 1.0 uses a jump kernel width corresponding to three times the asimov likelihood width.
 */
class asimov_widths_1d: public theta::MCMCStrategy{
public:
    explicit asimov_widths_1d(const theta::Configuration & producer_cfg);
    virtual void init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution);
    virtual void run_mcmc(const theta::Function & nllikelihood, theta::MCMCResult &res) const;
    
private:
    std::vector<double> widths;
};

#endif
