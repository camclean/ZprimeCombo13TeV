#ifndef MCMC_HPP
#define MCMC_HPP

#include "interface/decls.hpp"
#include "interface/matrix.hpp"
#include "interface/random-utils.hpp"

#include <vector>

#include <boost/shared_ptr.hpp>

namespace theta{

/** \brief Abstract base class which receives the result from a Markov Chain
 * 
 * An object of this class is passed as parameter to the Markov-Chain routines, which call the
 * fill method for each distinct point in the Markov Chain. Derived classes can implement
 * this method to perform whatever calculation they want with the information.
 */
class MCMCResult {
public:
    
    /// Returns the number of parameters
    virtual size_t getnpar() const = 0;
    
    /** \brief fill a new chain point with the given parameter values, nll value and weight
     *
     * This method is called by the metropolisHastings routine.
     * \param x contains getnpar() parameter values of the current point in parameter space
     * \param nll is the negative logarithm of the likelihood function at this point
     * \param weight is the multiplicity of this point in the Markov Chain, i.e., the number of rejected proposals to jump away from it, plus one.
     */
    virtual void fill(const double * x, double nll, size_t n) = 0;
    
    /// Virtual destructor for polymorphic access.
    virtual ~MCMCResult();
};


/** \brief MCMCResult class calculating mean and covariance
 * 
 * This class can be either used directly or used as a base class for derived classes
 * which want to provide additional functionality. Derived classes should re-implement
 * the fill2 method and leave the fill method unchanged, as the latter is used to perform
 * the mean and covariance calculation.
 */
class ResultMeanCov: public MCMCResult {
protected:
    /// number of parameters of the likelihood function
    size_t npar;
    
    /// number of total points in the chain (including rejected proposal points)
    size_t count;
    
    /// number of different points in the chain, i.e., not counting rejected proposal points
    size_t count_different_points;
    
    /// sliding mean of the parameter values in the chain
    std::vector<double> means;
    
    /// sliding covariance times count
    theta::Matrix count_covariance;
    
    /// fill method to be implemented by derived classes which want to save more than mean and covariance
    virtual void fill2(const double * p, double nll, size_t weight){}
    
public:
    /** \brief Construct result with \c npar parameters
     */
    explicit ResultMeanCov(size_t npar);
    
    virtual ~ResultMeanCov();
    virtual void fill(const double * x, double nll, size_t weight);
    virtual size_t getnpar() const;
    
    /// Returns the number of points in the chain, also counting rejected proposal points
    size_t getCount() const;
    
    /// Returns the number of different point in the chain, i.e., not including rejected proposals
    size_t getCountDifferent() const;
    
    /// Returns the mean of the parameter values in the chain
    std::vector<double> getMeans() const;
    
    /// Returns the covariance matrix of the parameter values in the chain
    theta::Matrix getCov() const;
};    

/// Small class to factor out common options which usually do not change for the Markov-Chain Monte-Carlo.
struct MCMCOptions{
    /// start values for the parameters
    std::vector<double> startvalues;
    
    /// number of 
    unsigned int iterations;
    
    /// number of burn-in steps
    unsigned int burn_in;
    
    /// additional factor for the step size; the meaning depends on the actual mcmc routine used.
    double factor;
    
    /// whether or not to ignore infinite values for the negative log-likelihood function at the start values.
    bool ignore_inf_nll_start;
    
    /// set some reasonable defaults for all members but startvalues, iterations, burn_in, which should *always* be set manually by whoever uses this class
    MCMCOptions(): factor(1.0), ignore_inf_nll_start(false){}
};


/** \brief Abstract class providing methods for initializing and running MCMC
 *
 * This encapsulates the start point, step size, jump kernel determination, number of burn-in and and
 * actual MCMC iterations, for which different strategies exist and can be used through this common interface.
 * 
 * It is usually used by producer plugins using a MCMC calculation internally. As such, it is usually configured
 * within a mcmc producer with a setting group name "mcmc_strategy", for example
 * \code
 *  producer = {
 *    type = "mcmc_quantiles";
 *    name = "quant";
 *    //... other settings for mcmc_quantiles ...
 *    mcmc_strategy = {
 *        type = "asimov_cov_mcmc";
 *        name = "strat1";
 *        iterations = 10000;
 *        burn-in = 100; // optional
 *        ignore_inf_nll_start = true; // optional, default is false
 *        factor = 2.0; // optional, default is 1.0
 *        // type-dependent parameters ...
 *    };
 *  };
 * \endcode
 * 
 * \c type is -- as usual -- the C++ class name of the concrete class to use
 * 
 * \c name is used to save the random seed in the rndinfo table.
 * 
 * \c iterations is the number of MCMC iterations (excluding burn-in)
 * 
 * \c burn-in is the number of burn-in iterations, i.e., iterations performed at the beginning of the chain which are not used to calculate the result.
 *   It defaults to iterations / 10.
 * 
 * \c ignore_inf_nll_start If set to false (the default), it is considered an error if the value of the negative log-likelihood is infinity at the
 *    start values, and a fatal exception will be thrown. If set to true, the Markov Chain will ignore this condition and pass all Markov-Chain
 *   elements to the result object. In this case, all jumps will be accepted until a point with non-infinite negative log-likelihood function value is found.
 * 
 * \c factor is an additional factor for the step of the Markov Chain. The default choice is factor = 1.0.
 */
class MCMCStrategy: public theta::RandomConsumer{
public:
    /// Typedef to make the plugin system work
    typedef MCMCStrategy base_type;
    
    /** \brief Initialize the parameters for the Markov Chain
     * 
     * This is called at least once before a call to run_mcmc, which performs the actual Markov-Chain integration.
     * It does whatever setup is necessary; usually this means determining the start values for the parameters
     * and the step sizes / jump kernel.
     */
    virtual void init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution) = 0;
    
    /** \brief Run the Markov-Chain Monte-Carlo algorithm
     *
     * Run the actual MCMC integration on the given negativ log-likelihood function, and pass the Markov Chain
     * elements to the given result object \c res.
     * 
     */
    virtual void run_mcmc(const theta::Function & nllikelihood, MCMCResult &res) const = 0;
    
    const MCMCOptions & get_options() const{
        return options;
    }
    
    /// Declare destructor virtual for correct polymorphic destruction
    virtual ~MCMCStrategy();    
protected:
    
    /// constructor to be used by derived classes. Reads "iterations", "burn-in", and "ignore_inf_nll_start" values; also uses "name" for the RandomConsumer
    explicit MCMCStrategy(const theta::Configuration & cfg);
    
    /// derived classes can use this options object for which all members but startvalues (!) is filled.
    MCMCOptions options;
};


/** \brief Strategy factory provided for backward compatibility.
 * 
 * The new way of constructing a MCMCSTrategy would be via its own setting group via the plugin system. However,
 * to more easily keep the old default behaviour of using the asimov_mcmc_cov strategy, this method is provided. It
 * will use the setting group in cfg.setting[setting_key], if present, to construct a MCMCSTrategy via the plugin system.
 * If the specified setting_key does not exist, the "old default", asimov_mcmc_cov, will be constructed, using the
 * producer_cfg as constructor argument. This should work as long as there is at least a setting "iterations".
 * 
 * New MCMC plugins should always use the plugin manager directly to construct the MCMCStrategy.
 */
std::auto_ptr<MCMCStrategy> construct_mcmc_strategy(const Configuration & producer_cfg, const std::string setting_key = "mcmc_strategy");


/** \brief MCMC using the empirical covariance from other MCMCs run on the likelihood from asimov data
 * 
 * Configured via:
 * \code
 * producer = {
 *     ...
 *     mcmc_strategy = {
 *         type = "asimov_mcmc_cov";
 *         name = "mcs";
 *         iterations = 10000;
 *         burn-in = 100;// optional
 *     };
 * };
 * \endcode
 * 
 * Note: for backward compatibility, where the MCMC parameters have been configured directly as settings of the
 * producer rather than forarding them as plugin, producers using MCMC can construct this type directly (rather than through the plugin system).
 * In this case, they should just apss through the configuration object, as this type only requires the settings
 * "name" and "iterations" -- which have been used so far anyway -- to reproduce the old behaviour.
 * 
 * This backward compatibility is the reason that this plugin is in the theta core rather than the plugin code.
 */
class asimov_mcmc_cov: public MCMCStrategy {
public:
    explicit asimov_mcmc_cov(const Configuration & producer_cfg);
    virtual void init(const Model & model, const boost::shared_ptr<Distribution> & override_parameter_distribution);
    virtual void run_mcmc(const Function & nllikelihood, MCMCResult &res) const;
private:
    Matrix sqrt_cov;
};




/** \brief Run the metropolis-Hastings Markov-Chain Monte-Carlo algorithm with a multivariate Gaussian jump kernel
 *
 * @param nll is the negative log-likelihood function to integrate
 * @param options specifies the startvalues, number of burn-in and iterations. The factor
 *   given thereis is to be understood as a factor w.r.t. the default choice of 2.38 / sqrt(n) * cov as jump kernel width
 *   where n is the number of non-fixed parameters.
 * @param res is the result class which will be used to write the Markov Chain to.
 * @param rand is the random number generator to use to compute the next candidate point of the chain.
 * @param sqrt_cov is (apart from a overall factor, see above) the matrix used as jump kernel. It should be set to the Cholesky decomposition of the
 *   covariance matrix of the likelihood (or an approximation thereof), hence the name "square root of covariance".
 */
void metropolis_hastings_multigauss(const theta::Function & nll, const MCMCOptions & options, MCMCResult & res, theta::Random & rand, const theta::Matrix & sqrt_cov);


/** \brief estimate the square root (cholesky decomposition) of the covariance matrix of the likelihood function
 *
 * The method will start a Markov chain at the given startvalues with the \c iterations iterations.
 * The found covariance matrix is used in a next pass where the Markov chain jumping rules
 * are adjusted to the covariance found in the previous step. This procedure is repeated until the
 * estimate is stable.
 *
 * \param rnd is the random number generator to use
 * \param model is the Model to use to get the asimov data from
 * \param[out] startvalues will contain the suggested startvalues. The contents when calling this function will be ignored.
 */
theta::Matrix get_sqrt_cov2(theta::Random & rnd, const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution);

/** \brief Calculate the cholesky decomposition, but allow zero eigenvalues.
 *
 *
 * Writes the cholesky decomposition into \c result. This function treats the case of fixing a parameter through
 * setting its covariance diagonal entry to zero correctly in only decomposing the non-zero part of the matrix
 * and keeping the zero entries where they were.
 *
 * If \c expect_reduced is non-negative, it will be checked whether it matches the determined number
 * of reduced (=non-fixed) dimensions. If it does not, an Exception will be thrown. If \c expect_reduced
 * is negative, no check will be done.
 */
void get_cholesky(const theta::Matrix & cov, theta::Matrix & result, int expect_reduced = -1);


} // namespace theta

#endif
