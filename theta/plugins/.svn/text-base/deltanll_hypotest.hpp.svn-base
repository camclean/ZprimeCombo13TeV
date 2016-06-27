#ifndef PLUGIN_DELTANLL_HYPOTEST_HPP
#define PLUGIN_DELTANLL_HYPOTEST_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/variables.hpp"
#include "interface/distribution.hpp"

#include <boost/optional.hpp>
#include <string>

/** \brief A producer to create test statistics based on likelihood ratio in case of signal search.
 *
 * This producer assumes that you search for a signal and can use special parameter values in your model to define
 * the "background only" and "signal plus background" hypotheses.
 *
 * Configuration is done via a setting group like
 * \code
 * hypotest = {
 *   type = "deltanll_hypotest";
 *   name = "deltanll";
 *   minimizer = "@myminuit";
 *   background-only-distribution = "@bkg-only-dist";
 *   signal-plus-background-distribution = "@default-dist";
 *   //optional:
 *   restrict_poi = "beta_signal";
 *   default_poi_value = 1.0;
 *   write_pchi2 = true; //default is false
 * };
 * 
 * myminuit = {...}; // minimizer definition
 * bkg-only-dist = {...}; //distribution definition
 * default-dist = {...}; //distribution definition
 * \endcode
 *
 * \c type is always "deltanll_hypotest" to select this producer.
 *
 * \c name is a unique producer name of your choice; it is used to construct column names in the output database. It may only contain alphanumeric
 *    characters (no spaces, special characters, etc.).
 *
 * \c minimizer is the configuration path to a \link theta::Minimizer minimizer \endlink definition to be used
 *    for minimization of the negative-log-likelihood.
 *
 * \c background-only-distribution and \c signal-plus-background-distribution define the Distributions which should be used
 *   in the two model variants. <b>Important:</b> the constraints will alter the likelihood function: instead
 *   of the parameter distribution specified in the likelihood function, these distributions will be used. Therefore,
 *   <em>all</em> modl parameters have to be specified for these two cases.
 *
 * \c restrict_poi is the parameter of interest whose range is restricted during minimization of the
 *  signal + background model. The allowed maximum value for thegiven parameter will be set to the current
 *  poi value, as set via setParameterValues.
 *
 * \c default_poi_value only applies if \c restrict_poi is set and is optional in this case. It is the default value to use if no
 * other value is set externally. If \c default_poi_value is not given and no value is set externallty, an exception is thrown.
 * 
 * If \c write_pchi2 is set to true, the chi2 value after the s+b fit is written to the result table.
 *   
 * Note that the setting "override-parameter-distribution" is not allowed for this producer.
 *
 * The result table will contain the columns "nll_sb" and "nll_b", and "nll_diff"
 * which contain the found value of the negative log-likelihood
 * for the "signal-plus-background" and "background-only" hypotheses, and the difference
 * of these two, nll_b - nll_sb, respectively. In case \c restrict_poi is set, an additional column "poi" contains the used value.
 *
 * For a typical application, the "signal-plus-background-distribution" setting is the same as in the model,
 * whereas the "background-only-distribution" setting group includes a delta_distributions which fixes
 * the signal to zero. Only this case is considered in the following.
 * Note that \c nll_sb <= \c nll_b <b>always</b> holds in this case as the minimization
 * is done using a larger set of parameters in the first case and the minimum cannot become larger. Any failure of this
 * inequality can only come from roundoff errors in the minimization process or from the
 * minimizer not finding the correct minimum.
 *
 * If the number of observed events is large, \code sqrt(2 * (nll_b - nll_sb)) \endcode will be a good estimate of
 * the significance (in sigma) with which the "background only" null-hypothesis "background-only" can be rejected.
 * Even if the asymptotic property is not fulfilled, this quantity can still be used as test statistic for the
 * hypothesis test which has the "background-only" case as null hypothesis.
 */
class deltanll_hypotest: public theta::ParameterDependentProducer{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    deltanll_hypotest(const theta::Configuration & cfg);
    virtual void produce(const theta::Data &, const theta::Model&);

    void set_parameter_values(const theta::ParValues & values);
private:
    
    boost::shared_ptr<theta::Distribution> s_plus_b;
    boost::shared_ptr<theta::Distribution> b_only;
    
    bool init;

    boost::optional<theta::ParId> restrict_poi;
    double poi_value;
    double default_poi_value;
    
    theta::ParValues s_plus_b_mode, b_only_mode;
    theta::ParValues s_plus_b_width, b_only_width;
    theta::Ranges s_plus_b_support, b_only_support;
    bool b_only_subset_s_plus_b;
    
    std::auto_ptr<theta::Minimizer> minimizer;
    
    bool write_pchi2;
    theta::Column c_nll_b, c_nll_sb, c_nll_diff, c_poi, c_pchi2;
};

#endif
