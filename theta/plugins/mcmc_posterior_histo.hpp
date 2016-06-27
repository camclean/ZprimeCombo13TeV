#ifndef PLUGIN_MCMC_POSTERIOR_HISTO_HPP
#define PLUGIN_MCMC_POSTERIOR_HISTO_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/random-utils.hpp"
#include "interface/matrix.hpp"

#include <string>

/** \brief Construct a Histogram of the marginalized Posterior in one parameter based on Markov-Chain Monte-Carlo
 *
 * Typically, this producer is called with a large number of iterations for only few pseudo
 * experiments.
 *
 * Configuration is done via a setting group like
 * \code
 * post = {
 *   type = "mcmc_posterior_histo";
 *   name = "post";
 *   parameters = ("s");  //assuming "s" was defined as parameter earlier
 *   iterations = 100000;
 *   burn-in = 100; //optional. default is iterations / 10 
 *   smooth = true; //optional, default is false
 *
 *   histo_s = {
 *      range = [0.0, 100.0];
 *      nbins = 100;
 *   };
 * };
 *
 * \endcode
 *
 * \c type is always "mcmc_posterior_histo" to select this producer.
 *
 * \c name is a unique producer name of your choice; it is used to construct column names in the output database. It may only contain alphanumeric
 *    characters (no spaces, special characters, etc.).
 *
 * \c parameters is a list of parameter names you want to calculate the posterior Histograms for
 *
 * \c iterations is the number of MCMC iterations. See additional comments about runtime and suggested robustness tests
 *     in the documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink.
 *
 * \c burn_in is the number of MCMC iterations to do at the beginning and throw away. See additional comments in the
 *     documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink
 *
 * \c smooth controls whether to make a smooth posterior histogram. In this case, each point in the chain contributes
 *    with a whole histogram in the respective parameter instead of with a single point. This has important consequences
 *    on the runtime, see comment below.
 *
 * For each parameter given in the \c parameters setting, you must specify a setting group
 * of name "histo_<parameter name>" which has two settings:
 * <ul>
 * <li>\c range is an array of two floating point numbers defining the range of the Histogram</li>
 * <li>\c n_bins is an integer defining the number of bins of the Histogram.</li>
 * </ul>
 *
 * For each data given, one chain with the given number of iterations is constructed
 * and the value of the parameters is filled in the histogram. The histogram is
 * written to a column "posterior_<parameter name>".
 *
 * The runtime per dataset is driven by the number of likelihood evaluations. There is an additional
 * initialization overhead for the first toy.
 * 
 * Without smoothing (\c smooth = false, the default), the number of likelihood evaluations is about 
 * \c iterations + \c burn-in. It does not depends on the number of histograms
 * or the number of bins per histogram.
 * 
 * With smoothing, the number of likelihood evaluations increases with respect to the non-smooth version
 * by \c iterations * \c n_nbins * acceptance rate, where "acceptance rate" is the probability of a proposal point to be
 * accepted and is typically between 0.2 and 0.4; \c n_bins is the total number of bins of all histograms requested.
 *
 * Important: if switching smoothing on, the provided range has to be large enough to include the whole posterior.
 */
class mcmc_posterior_histo: public theta::Producer {
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    mcmc_posterior_histo(const theta::Configuration & ctx);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    
    //whether sqrt_cov* and startvalues* have been initialized:
    bool init;
    
    std::vector<theta::ParId> parameters;
    std::vector<std::string> parameter_names;
    std::vector<size_t> ipars; //parameters of the requested index, as in NLLikelihood::operator()(const double*) index convention
    
    //result columns: one per requested parameter:
    std::vector<theta::Column> columns;
    std::vector<double> lower, upper;
    std::vector<size_t> nbins;
    
    std::auto_ptr<theta::MCMCStrategy> mcmc_strategy;
    
    bool smooth;
};

#endif
