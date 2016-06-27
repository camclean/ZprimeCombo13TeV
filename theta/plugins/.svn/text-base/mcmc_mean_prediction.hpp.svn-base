#ifndef PLUGINS_MCMC_MEAN_PREDICTION_HPP
#define PLUGINS_MCMC_MEAN_PREDICTION_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/random-utils.hpp"
#include "interface/matrix.hpp"

#include <string>

/** \brief Construct the templates for the mean prediction, averaged over the postrior, using MCMC
 *
 * The resultin templates can be used to visualize the (per-bin) mean and width of the MCMC
 * integration and thus the compatibility of the prediction with data considering systematic
 * uncertainties.
 *
 * Configuration is done via a setting group like
 * \code
 * mean_pred = {
 *   type = "mcmc_mean_prediction";
 *   name = "mp";
 *   observables = ("o1", "o2", "o3");
 *   iterations = 10000;
 *   burn-in = 100; //optional. default is iterations / 10
 * };
 * \endcode
 *
 * \c type is always "mcmc_mean_prediction" to select this producer.
 *
 * \c name is a name chosen by the user used to construct unique column names in the result table (this name and two underscores are
 *   prepended to the column names explained below).
 *
 * \c observables is a list of observables for which to calculate the mean and width histograms
 *
 * \c iterations is the number of MCMC iterations. See additional comments about runtime and suggested robustness tests
 *     in the documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink.
 *
 * \c burn_in is the number of MCMC iterations to do at the beginning and throw away. See additional comments in the
 *     documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink
 *
 * For given data, one Markov Chain is constructed. For each point in the Markov Chain,
 * the prediction of the Poisson mean in each bin of each observable is calculated. This mean
 * enters the calculation of the mean and width of the prediction. In addition, the prediction at the highest
 * likelihood value is saved.
 *
 * The results are written to the 'products' table. In this table, 
 * columns per observable of type \link theta::Table::typeHisto \endlink are created.
 * The column names for an observable "obs1" will be "obs1_mean" and "obs1_width" where the former contains the mean
 * and the latter the rms of the predictions in each bin, as seen along the chain. The third column is called
 * "obs1_best" and contains the prediction with the highest likelihood value the chain has seen.
 */
class mcmc_mean_prediction: public theta::Producer {
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    mcmc_mean_prediction(const theta::Configuration & ctx);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    //result columns: one "mean" and one "width" column per observable.
    theta::ObsIds observables;
    std::vector<theta::Column> c_mean;
    std::vector<theta::Column> c_width;
    std::vector<theta::Column> c_best;
    
    std::auto_ptr<theta::MCMCStrategy> mcmc_strategy;
    bool init;
};

#endif
