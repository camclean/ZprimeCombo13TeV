#ifndef PLUGINS_MCMC_POSTERIOR_RATIO_HPP
#define PLUGINS_MCMC_POSTERIOR_RATIO_HPP

#include "interface/decls.hpp"
#include "interface/producer.hpp"
#include "interface/random-utils.hpp"
#include "interface/mcmc.hpp"

#include <string>

/** \brief A producer to create test statistics based on the ratio of the posterior in case of signal search.
 *
 * This producer assumes that you search for a signal and can use special parameter values in your model to define
 * the "background only" and "signal plus background" hypotheses.
 *
 * This producer is similar to the \link deltanll_hypotest deltanll_hypotest \endlink method but integrates the likelihood
 * over free parameters instead of minimizing it.
 *
 * Configuration is done via a setting group like
 * \code
 * hypotest = {
 *   type = "mcmc_posterior_ratio";
 *   name = "ratio";
 *   background-only-distribution = "@bkg-only-dist";
 *   signal-plus-background-distribution = "@default-dist";
 *   iterations = 10000;
 *   burn-in = 100; //optional. default is iterations / 10
 * };
 * \endcode
 *
 * \c type is always "mcmc_posterior_ratio" to select this producer.
 *
 * \c name is a unique producer name of your choice; it is used to construct column names in the output database. It may only contain alphanumeric
 *    characters (no spaces, special characters, etc.).
 *
 * \c background-only-distribution and \c signal-plus-background-distribution deinfe the Distribution instances to use for the
 *   two MCMC integrations.
 *
 * \c iterations is the number of MCMC iterations. Run time is proportional to this values; increasing this value will result in higher precision.
 *     The optimal value to use here is very much problem-dependent. To find out a reasonable value, it is advised to start with a value around 5000 and
 *     run multiple independent  mcmc_posterior_ratio producers which all use the same \c iterations setting. You can then compare the results of
 *     the different producers. As a rule of thumb, the inter-chain variance should be at least two orders of magnitude smaller
 *     compared to the inner-chain variance (i.e., of the posterior_s and posterior_sb within a run. Usually these two have similar variance).
 *
 * \c burn_in is the number of MCMC iterations to do at the beginning and throw away. There is some controversy about whether it is needed at all
 *     and how large it should be. If unsure, just take the default value of iteratios / 10 and if someone asks why, vary it a bit (say, between iterations/100 and iterations)
 *     and show that the result does not depend on it anyway ...
 *     
 * Note that the setting "override-parameter-distribution" is not allowed for this producer.
 *
 * Given data and a model, this producer will construct the negative-loglikelihood for the "signal-plus-background" parameters
 * fixed as specified in the configuration file and for the "background-only" parameters fixed and integrate over all non-fixed parameters.
 *
 * The negative logarithm of the found average values of the likelihood (which takes the role of a posterior here) are
 * saved in the \c nl_posterior_sb and \c nl_posterior_b columns of the result table.
 */
class mcmc_posterior_ratio: public theta::Producer, public theta::RandomConsumer{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    mcmc_posterior_ratio(const theta::Configuration & ctx);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    
    //whether sqrt_cov* and startvalues* have been initialized:
    bool init;
    
    boost::shared_ptr<theta::Distribution> s_plus_b;
    boost::shared_ptr<theta::Distribution> b_only;
    
    unsigned int iterations;
    unsigned int burn_in;
    
    //the matrices and startvalues to use for the Markov chains in the two cases:
    theta::Matrix sqrt_cov_sb, sqrt_cov_b;
    theta::MCMCOptions options_sb, options_b;
    
    theta::Column c_nl_posterior_sb, c_nl_posterior_b;
};

#endif
