#ifndef PLUGIN_DELTANLL_INTERVALS_HPP
#define PLUGIN_DELTANLL_INTERVALS_HPP

#include "interface/plugin.hpp"

#include "interface/variables.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/distribution.hpp"

#include <string>

/** \brief Confidence levels based on likelihood ratio
 *
 * Producer to construct confidence intervals based on asymptotic properties of the likelihood-ratio,
 * or, equvalently, the difference in the negative log-likelihood.
 *
 * Configuration is done with a settings block like:
 * \code
 * {
 *  type = "deltanll_intervals";
 *  name = "intervals";
 *  parameter = "p0";
 *  minimizer = "@myminuit";
 *  clevels = [0.68, 0.95];
 *  re-minimize = false; //Optional. Default is true
 * }
 *
 * myminuit = {...} //see the minimizer documentation.
 * \endcode
 *
 * \c type has always to be "deltanll_intervals" in order to use this producer
 *  
 * \c name is a unique producer name of your choice; it is used to construct column names in the output database. It may only contain alphanumeric
 *    characters (no spaces, special characters, etc.).
 *
 * \c parameter is the name of the parameter for which the interval shall be calculated.
 *
 * \c minimizer is the configuration path to a \link theta::Minimizer Minimizer\endlink definition. It does
 *   not need to provide error estimates although error estimates can help to speed up calculation in some cases.
 *
 * \c re-minimize specifies whether or not to search for a minimum of the negative log-likelihood when scanning
 *    through the parameter of interest or to use the parameter values found at the global minimum. See below for details.
 *
 * \c clevels is an array (or list) of doubles of confidence levels for which the intervals
 *   shall be calculated. Note that an interval for the "0" confidence level (i.e., the
 *   interval containing the minimum) is always determined.
 *
 * This producer uses the likelihood to derive confidence intervals based on the
 * based on the asymptotic property of the likelihood ratio test statistics to be distributed
 * according to a chi^2-distribution.
 *
 * Given a likelihood function \f$ L(\vec p) \f$ which depends on parameters \f$ \vec p = (p_0, ..., p_n) \f$, the
 * profile likelihood function \f$ L_{\mathrm{p}}(p_0)\f$ for parameter \f$ p_0 \f$ is constructed as follows:
 * <ol>
 * <li>maximize \f$ L(\vec p) \f$ to find the parameter values at the maximum \f$ \hat{\vec p} \f$. </li>
 * <li>scan through \f$ p_0 \f$ to find \f$ L_{\mathrm{p}}(p_0)\f$ by either (i) maximizing \f$ L(\vec p) \f$
 *     fixing only \f$ p_0 \f$, or (ii) using the
 *     values found in step 1 for all parameters \f$ \hat{p}_i \f$ but replace \f$ p_0 \f$. Which method is used depends
 *     on the \c re-minimize setting.
 *     The found value is divided by the maximum likelihood value \f$ L(\hat{\vec p}) \f$ such that the maximum of
 *     \f$ L_{\mathrm{p}}(p_0) \f$ is always 1.</li>
 * </ol>
 *
 * The result table always contains the parameter value at the maximum of the likelihood as "maxl". The "maxl" value
 * is the same as you would get for a confidence level of 0.
 *
 * If interval construction is requested through non-empty \c clevels list, the lower and upper bounds are saved which
 * are determined by the point where \f$ L_{\mathrm{p}}(p_0) \f$ takes the values corresponding to these confidence levels determined
 * via the asymptotic chi^2 distribution property of likelihood ratios.
 * For each confidence level c, the result table contains columns "lower" + 10000*c and "upper" + 10000*c, where
 * the numbers are rounded and written with leading zeros. For example, if the \c clevels setting is [0.68, 0.95], the
 * column names will be "lower06800", "upper06800", "lower09500", "upper09500".
 */
class deltanll_intervals: public theta::Producer{
public:

    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    deltanll_intervals(const theta::Configuration & cfg);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    void declare_products();
    
    theta::ParId pid;
    std::vector<double> clevels;
    bool re_minimize;
    //at construction, save here the deltanll values corresponding to
    //clevels:    
    std::vector<double> deltanll_levels;
    std::auto_ptr<theta::Minimizer> minimizer;
    
    bool start_step_ranges_init;
    theta::ParValues start, step;
    theta::Ranges ranges;

    //table columns:
    std::vector<theta::Column> lower_columns;
    std::vector<theta::Column> upper_columns;
    theta::Column c_maxl;
};

#endif

