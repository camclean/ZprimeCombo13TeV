#ifndef PLUGINS_MCMC_QUANTILES_HPP
#define PLUGINS_MCMC_QUANTILES_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/producer.hpp"
#include "interface/mcmc.hpp"

#include <string>
#include <memory>

/** \brief Construct quantiles of the marginal posterior in one parameter based on Markov-Chain Monte-Carlo
 *
 * The result can be used to give upper limits or to construct symmetric credible intervals.
 *
 * Configuration is done via a setting group like
 * \code
 * hypotest = {
 *   type = "mcmc_quantiles";
 *   name = "quant";
 *   parameter = "s";  //assuming "s" was defined as parameter earlier
 *   quantiles = [0.025, 0.16, 0.5, 0.84, 0.975];
 *   iterations = 10000;
 *   burn-in = 100; //optional. default is iterations / 10
 *   re-init = 1; //optional. Default is 0
 * };
 *
 * \endcode
 *
 * \c type is always "mcmc_posterior_ratio" to select this producer.
 *
 * \c name is a unique producer name of your choice; it is used to construct column names in the output database. It may only contain alphanumeric
 *    characters (no spaces, special characters, etc.).
 *
 * \c parameter is the name of the parameter you want to find the quantiles for
 *
 * \c quantiles is a list or array of floating point values specifying the quantiles you want. Be aware that
 *   specifying very similar quantiles (with a difference less than 0.00001) is not allowed as this would yield same column
 *   name in the result table; see below.
 *
 * \c iterations is the number of MCMC iterations. See additional comments about runtime and suggested robustness tests
 *     in the documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink.
 *
 * \c burn-in is the number of MCMC iterations to do at the beginning and throw away. See additional comments in the
 *     documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink
 *
 * \c re-init is an optional integer which controls re-initialisation of the jumping kernel width. The deault of 0 never re-initialises. For a value N > 0, 
 * re-initialisation is done every N toys.
 *
 * For each data given, one chain will be used to derive all requested quantiles given in the \c quantiles list, so their error
 * from limited chain length is correlated by construction. If you do not want that, use two independent producers of type
 * mcmc_quantiles.
 *
 * The result table contains as many columns as \c quantiles given in the configuration file. The column name
 * will be "quant" + 10000 * quantile, written with leading zeros. For example, if the quantile is 0.5,
 * the column name will be "quant05000", if the 99.9% quantile is requested (i.e., 0.999), the name will be "quant09990".
 * 
 * In addition to the quantiles, the products table will also have a column "accrate" with the acceptance rate for
 * the Markov Chain. This can be used as diagnostical tool in case the performance is not as good as expected: the
 * acceptance rate should be around 15-30%; larger acceptance rates might indicate a step size that's too small,
 * while smaller acceptance rates might indicate a step size chi is too large.
 */
class mcmc_quantiles: public theta::Producer {
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    mcmc_quantiles(const theta::Configuration & ctx);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    void declare_products();
    
    //whether sqrt_cov* and startvalues* have been initialized:
    bool init;
    
    //the requested quantiles:
    std::vector<double> quantiles;
    theta::ParId par_id;
    size_t ipar; //parameter of the requested index, as in NLLikelihood::operator()(const double*) index convention

    int re_init, itoy;
    
    //result columns: one per requested quantile:
    std::vector<theta::Column> columns;
    theta::Column c_accrate;
    
    std::auto_ptr<theta::MCMCStrategy> mcmc_strategy;
};

#endif
