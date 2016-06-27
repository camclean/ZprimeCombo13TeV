#ifndef PLUGINS_MCMC_QUANTILES_HPP
#define PLUGINS_MCMC_QUANTILES_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/random-utils.hpp"
#include "interface/matrix.hpp"

#include <string>

/** \brief Write out the Markov Chain in text files
 *
 * The result can be used to construct many informations from the chain directly, such as
 * multidimensional posteriors, etc.
 *
 * It will create one text file per run.
 *
 * Configuration is done via a setting group like
 * \code
 * chain = {
 *   type = "mcmc_chain";
 *   iterations = 10000;
 *   burn-in = 100; //optional. default is iterations / 10
 *   re-init = 1; //optional. Default is 0
 *
 *   // option 1: write into text files
 *   outfile_prefix = "chain";
 *
 *   // option 2: write into a database
 *   output_database = {
 *      type = "rootfile_database";
 *      filename = "out.root";
 *   };
 * };
 *
 * \endcode
 *
 * \c type is always "mcmc_chain" to select this producer.
 *
 * \c iterations is the number of MCMC iterations. See additional comments about runtime and suggested robustness tests
 *     in the documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink.
 *
 * \c burn-in is the number of MCMC iterations to do at the beginning and throw away. See additional comments in the
 *     documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink
 *
 * \c re-init is an optional integer which controls re-initialisation of the jumping kernel width. The deault of 0 never re-initialises. For a value N > 0, 
 *   re-initialisation is done every N toys.
 *
 * There are two different options how to write out the result. You have to use exactly one of those options by providing either the \c outfile_prefix or
 * the \c output_database option.
 *
 * The first option is to write out one text file per toy. This text file contains
 * (n+1) columns if n is the number of parameter the likelihood depends on; the extra column is used to save the weight. This option is selected by specifying
 * a value of \c output_prefix; \c outfile_prefix is the filename prefix for the produces .txt files. The files will be constructed
 * by the \c outfile_prefix given here, and an incrementing counter, separated by an underscore ("_"). The filename ends with ".txt".
 *
 * The second option is selected by providing the \c output_database option. In this case, any
 * database plugin from theta can be used such as \c rootfile_database or \c sqlite_database. In this case, one table is created
 * per dataset with the name "chain_&lt;n&gt;" where &lt;n&gt; is a counter starting at 1. Each table is structured like the text file and has (n+1) columns
 * in case of n parameters, the first being the "weight" of the chain element, the others the values of the parameters for this chain element.
 *
 * In general, it is recommended to use the \c output_database as this will usually result in fewer files, smaller total filesize, and faster i/o. For testing --
 * or depending on the which software is used to process the output -- you can also use the text file output.
 *
 * <b>Important:</b> This plugin does not write anything to the output database. As it uses its own way of writing results, it somewhat breaks some concepts that hold
 * for most other producers, i.e. (i) the result is no longer only in the one database file configured in \c main.output_database, but also in the
 * text files/the additional database file; (ii) thread safety when using \c main.type="run_mt" is not guaranteed for this producer, as the creation of output will lead
 * to a race condition.
 */
class mcmc_chain: public theta::Producer {
public:
    mcmc_chain(const theta::Configuration & ctx);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    
    static std::string construct_name();
    
    //whether sqrt_cov* and startvalues* have been initialized:
    bool init;
    

    int re_init, itoy;
    boost::shared_ptr<theta::VarIdManager> vm;
    std::vector<std::string> parameter_names;
    
    // used for text file output:
    std::string outfile_prefix;
    
    // used for database output:
    boost::shared_ptr<theta::Database> database;
    
    //MCMC parameters:
    std::auto_ptr<theta::MCMCStrategy> mcmc_strategy;
};

#endif
