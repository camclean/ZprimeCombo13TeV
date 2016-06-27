#ifndef PLUGINS_ASYMPTOTIC_CLS
#define PLUGINS_ASYMPTOTIC_CLS

#include "interface/decls.hpp"
#include "interface/main.hpp"
#include "interface/variables.hpp"


/** \brief Calculate Asymptotic CLs limits
 * 
 * Calculate expected and observed CLs limits according to the formulas given in
 * http://arxiv.org/abs/1007.1727 with the test statistic for upper limits q_mu-tilde.
 * It uses the likelihood ratio of the asimov dataset for the expected test statistic distribution (called
 * "sigma_A" in the paper).
 *
 * \code
 * main = {
 *   type = "asymptotic_cls";
 *   parameter = "beta_signal";
 *   cl = 0.9; // optional, default is 0.95
 *   model = "@model";
 *   minimizer = { ... };
 *   limit_reltol = 1e-5; // optional, default is 1e-3
 * 
 *   //options for observed limits:
 *   data = {...}; // optional
 *   n = 10; // optional, default is 1
 * 
 *   // options for expected limits:
 *   quantiles_expected = (0.025, 0.16, 0.5, 0.84, 0.975); // optional; this is the default
 *   parameter_value_expected = 1.0; // optional, default is 0.0
 * };
 * \endcode
 * 
 * \c parameter is the parameter name to calculate the asymptotic CLs limits for.
 * 
 * \c cl is the confidence level of the limit to calculate. It has to be larger than 0.6 (the code assumes it to be larger than 0.5 at many places; 0.6
 *    includes some safety margin.).
 * 
 * \c model is the statistical model
 * 
 * \c minimizer is a Minimizer specification for the minimizer to use for the profile likelihood ratio test statistic calculation
 * 
 * \c limit_reltol is a tolerance parameter used in the search of the CLs limit calculation: The CLs limit is searched for by
 *   a 1D root finding algorithm (Brent's algorithm). The search is stopped if the interval containing the true CLs limit is smaller
 *   than \c limit_reltol * \c width, where \c width is the expected uncertainty of the parameter of interest.
 * 
 * \c data is a DataSource specification for the observed data. If missing, no "observed" limit will be calculated.
 * 
 * \c n is ther number of "observed" limits to calculate. This is usually \c 1, but in case the DataSource \c data contains e.g.
 *   many toy datasets, you can specify a larger number here.
 * 
 * \c quantiles_expected are the quantiles of the expected limit distribution to calculate. The default setting given above corresponds to the typical
 *   configuration for the median, central 1sigma, and central 2sigma intervals. All values must be &gt; 0 and &lt; 1. Can be set to an empty list
 *   to suppress calculation of expected limits.
 * 
 * \c parameter_value_expected is the value of the \c parameter to be used for the expected limit calculation. The default of 0.0
 *   means to compute the expected limit without any signal.
 * 
 * The database will contain two tables: a rndinfo table, and a table containing the actual limits,
 * "limits". This table contains three columns: "q" is the quantile as configured in \c quantiles_expected or 0.0 for the observed limit;
 * "limit" is the asymptotic limit; "index" is the 0-based index for \c quantiles_expected. The observed limit is the last, so its
 * index will be 5 for the default \c quantiles_expected.
 * 
 * The reported progress is the number of calculated asymptotic limits. So if calculating both expected and one observed limit, the reported total is 6.
 * 
 * If a MinimizationException occurs during evaluating the expected limit, theta aborts, while for the observed limit(s), such a failure is allowed
 * and the corresponding "index" is skipped. This is useful if calculating "observed" limits for a large ensemble where a failure rate of one permille
 * or so is usually not a problem.
 */
class asymptotic_cls: public theta::Main{
public:
    explicit asymptotic_cls(const theta::Configuration & cfg);
    virtual void run();
    
private:
    theta::ParId parameter;
    double cl;
    std::auto_ptr<theta::Model> model;
    std::auto_ptr<theta::Minimizer> minimizer;
    
    std::auto_ptr<theta::DataSource> data;
    int n;
    
    std::vector<double> quantiles_expected;
    double parameter_value_expected;
    double limit_reltol;
    
    boost::shared_ptr<theta::Database> output_database;
    std::auto_ptr<theta::Table> limits_table;
};


#endif
