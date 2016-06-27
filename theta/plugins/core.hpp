#ifndef PLUGIN_CORE_HPP
#define PLUGIN_CORE_HPP

#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"
#include "interface/phys.hpp"
#include "interface/distribution.hpp"
#include "interface/matrix.hpp"
#include "interface/database.hpp"
#include "interface/random-utils.hpp"

#include <boost/optional.hpp>

/** \brief A polynomial distribution where coefficients do not depend on any parameters
 *
 * Configuration is done with a setting group like
 * \code
 * {
 *  type = "fixed_poly";
 *  observable = "mass";
 *  normalize_to = 1.0;
 *  coefficients = [1.0, 2.0, 5.0];
 *  relative_bb_uncertainty = 0.1; // optional, default is 0.0
 * };
 * \endcode
 *
 * \c type must always be "fixed_poly" to construct an instance of this type
 *
 * \c observable is the name of a defined observable. This is used to construct a Histogram with the correct range and binning
 *
 * \c normalize_to is the sum of bin contents the returned histogram should have
 *
 * \c coefficients is an array (or list) of floating point values which define the polynomial, starting at x^0. The example above defines
 *  a polynomial 1 + 2*x + 5*x^2
 *
 *  \c relative_bb_uncertainty is the relative bin-by-bin uncertainty. The default is 0.0, i.e., no uncertainties.
 */
class fixed_poly: public theta::ConstantHistogramFunction{
public:
    fixed_poly(const theta::Configuration & cfg);
};

/** \brief A normal distribution where mean and width do not depend on any parameters
 *
 * Configuration is done with a setting group like
 * \code
 * {
 *  type = "fixed_gauss";
 *  observable = "mass";
 *  normalize_to = 1.0;
 *  mean = 1.0;
 *  width = 0.2;
 *  relative_bb_uncertainty = 0.1; // optional, default is 0.0
 * };
 * \endcode
 *
 * \c type must always be "fixed_gauss" to construct an instance of this type
 *
 * \c observable is the name of a defined observable. This is used to construct a Histogram with the correct range and binning
 *
 * \c normalize_to is the sum of bin contents the returned histogram should have
 *
 * \c mean and \c width are the mean value and standard deviation for the distribution to construct.
 *
 *  \c relative_bb_uncertainty is the relative bin-by-bin uncertainty. The default is 0.0, i.e., no uncertainties.
 */
class fixed_gauss: public theta::ConstantHistogramFunction{
public:
   fixed_gauss(const theta::Configuration & cfg);
};

/** \brief A lognormal distribution in one dimension.
 *
 * It is configured with a setting group like
 * \code
 * {
 *  type = "log_normal";
 *  parameter = "p0";
 *  mu = 2.0;
 *  sigma = 0.5;
 * };
 * \endcode
 *
 * \c parameter specifies the parameter the normal distribution depends on
 *
 * \c mu and \c sigma are floating point constants used to define the distribution.
 *
 * In the parametrization used, the density is proportional to:
 *   \f$ \exp( - (\ln(x) - \mu)^2 / (2 * \sigma^2) ) \f$
 * for x > 0 and 0 otherwise. x is the parameter the density depends on, mu and sigma are the configuration parameters.
 *
 */
class log_normal: public theta::Distribution{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    log_normal(const theta::Configuration & cfg);
    
    ///@{
    /** \brief Implementation of the pure methods of theta::Distribution
     *
     * See documentation of theta::Distribution.
     */
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual const std::pair<double, double> & support(const theta::ParId&) const;
    ///@}
private:
    
    double mu, sigma;
    std::pair<double, double> support_;
};

/** \brief The delta distribution in one or more dimensions
 *
 * This distribution always returns the same values when sampled, which is equal to the mode.
 * The width is zero, the support only contains the fixed value.
 * 
 * Configuration is done via a setting group like
 * \code
 * {
 *   type = "delta_distribution";
 *   s = 1.0; 
 *   b = 2.0; //assuming "s" and "b" have been declared as parameters
 * };
 * \endcode
 * 
 * \c type must always be "delta_distribution" to create an instance of this class.
 * 
 * Further, for every parameter, the fixed value is given.  
 */
class delta_distribution: public theta::Distribution{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    delta_distribution(const theta::Configuration & cfg);
    
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual const std::pair<double, double> & support(const theta::ParId&) const;
    virtual double eval_nl_with_derivative(const theta::ParValues & values, theta::ParValues & derivative) const;
    
private:
    theta::ParValues values;
    std::map<theta::ParId, std::pair<double, double> > supports;
};

/** \brief A flat Distribution
 * 
 * This distribution sets ranges for parameters and returns all parameters with the same
 * probability on the given range.
 * 
 * Configuration is done via a setting group like
 * \code
 * distribution = {
 *   type = "flat_distribution";
 *   s = {
 *      range = (0.0, "inf");
 *      fix-sample-value = 7.0; // must be given here, as the range is infinite
 *   };
 *   
 *   b = {
 *     range = (0.0, 5.0);
 *     fix-sample-value = 2.0; //optional for finite range
 *   };
 * };
 * \endcode
 * 
 * \c type is always "flat_distribution" to select this type.
 * 
 * For every parameter this Distribution will be defined for, a setting group is given with
 * <ul>
 *   <li>a \c range setting which specifies, as a list of two doubles, the range of the variable.
 *     The special strings "inf" and "-inf" are allowed here.</li>
 *   <li>A \c fix-sample-value setting which will be used by the sample routine. In case of finite
 *     intervals, the default is to sample from a flat distribution on this interval. For infinite intervals, \c fix-sample-value must be given and is the
 *     value returned when sampling from the distribution.</li>
 * </ul>
 * Note that \c fix-sample-value is also used as the mode (=most probable value) for this parameter and as such
 * is used as starting point for minimization, MCMC, etc. Therefore, it is advisable to set \c fix-sample-value even for
 * finite intervals to make sure the this value makes sense. In case of finite intervals, the mode is the midpoint of the interval.
 */
class flat_distribution: public theta::Distribution{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    flat_distribution(const theta::Configuration & cfg);
    ///@{
    /** \brief Implementation of the pure methods of theta::Distribution
     *
     * See class documentation and documentation of theta::Distribution for details.
     * Note that width, sample and mode might throw an exception for infinite intervals, as discussed
     * in the class documentation.
     */
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual const std::pair<double, double> & support(const theta::ParId&) const;
    virtual double eval_nl_with_derivative(const theta::ParValues & values, theta::ParValues & derivative) const;
    ///@}
    
    
private:
    theta::ParValues fix_sample_values;
    theta::ParValues modes;
    std::map<theta::ParId, std::pair<double, double> > ranges;
};

/** \brief A (possibly truncated) normal distribution in one or more dimensions
 *
 * A rectangular-truncated normal distribution.
 *
 * A one-dimensional case is configured with a setting group like 
 * \code
 * { 
 *  type = "gauss";
 *  parameter = "p0";
 *  range = (0.0, 5.0);
 *  mean = 2.0;
 *  width = 0.5;
 * };
 * \endcode
 *
 * \c parameter specifies the parameter the normal distribution depends on
 *
 * \c range is a list of doubles (or the special strings "inf", "-inf") specifying the truncation range.
 *
 * \c mean is a floating point value specifying the mean value of the distribution, \c width is its standard deviation.
 *
 * A multi-dimensional normal distribution can be specified with a setting like 
 * \code
 * { 
 *  type = "gauss";
 *  parameters = ("p0", "p1");
 *  ranges = ((0.0, 5.0), (0.0, "inf"));
 *  mean = [2.0, 3.0];
 *  covariance = ([1.0, 0.2], [0.2, 1.0]);
 * } 
 * \endcode
 *
 * \c parameters is a list of parameters this normal distribution depends on
 *
 * \c ranges is a list of ranges, in the same order as the parameters, which control the truncation
 *
 * \c mean specifies, in the same order as parameters, the mean values to use for the gaussian. 
 * 
 * \c covariance is the (symmetric) covariance matrix for the normal distribution. Note that 
 *     you give it as list of arrays (as it is symmetric anyway, it is left open what the "rows" and "columns" are). 
 */
class gauss: public theta::Distribution{
   public:
        /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
        gauss(const theta::Configuration & cfg);

        virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
        virtual void mode(theta::ParValues & result) const;
        virtual double eval_nl(const theta::ParValues & values) const;
        virtual const std::pair<double, double> & support(const theta::ParId&) const;
        
    private:
        
        std::vector<theta::ParId> v_par_ids;
        std::vector<double> mu;
        theta::Matrix sqrt_cov; //required for sampling
        theta::Matrix inverse_cov;//required for density evaluation
        std::vector<std::pair<double, double> > ranges;

        // temporary variables for sampling:
        mutable std::vector<double> x, x_trafo;
};

/** \brief A one-dimensional gauss-distribution
 *
 * It is configured via a setting group like
 * \code
 * gauss = {
 *  type = "gauss1d";
 *  parameter = "p0";
 *  range = (0.0, 5.0);
 *  mean = 2.0;
 *  width = 0.5;
 * };
 * \endcode
 *
 * \c parameter specifies the parameter the normal distribution depends on
 *
 * \c range is a list of doubles (or the special strings "inf", "-inf") specifying the truncation range.
 *
 * \c mean is a floating point value specifying the mean value of the distribution, \c width is its standard deviation.
 *
 * As special case, the \c mean can also be a string. In this case, it is interpreted as the name of a parameter which
 * will be used as mean.
 */
class gauss1d: public theta::Distribution{
   public:
        /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
        gauss1d(const theta::Configuration & cfg);

        virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
        virtual void mode(theta::ParValues & result) const;
        virtual double eval_nl(const theta::ParValues & values) const;
        virtual const std::pair<double, double> & support(const theta::ParId&) const;
        
    private:
        double mu;
        double sigma;
        boost::optional<theta::ParId> mu_pid;
        theta::ParId pid;
        std::pair<double, double> range;
};

/** \brief A Distribution product of other distributions
 *
 * Product of zero or more Distributions. Configured via a setting group
 * like
 * \code
 *  {
 *     type = "product_distribution";
 *     distributions = ("@d1", "@d2");
 *  };
 *
 *  d1 = {...}; //some Distribution
 *  d2 = {...}; //some other Distribution
 * \endcode
 *
 * \c distributions contains a list of distribution specifications.
 *
 * It is not allowed that a parameter is provided by more than one distribution.
 */
class product_distribution: public theta::Distribution{
public:

    /// Constructor from a Configuration for the plugin system
    product_distribution(const theta::Configuration & cfg);

    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual double eval_nl_with_derivative(const theta::ParValues & values, theta::ParValues & derivative) const;
    virtual const std::pair<double, double> & support(const theta::ParId & p) const;
    

private:
    void add_distributions(const theta::Configuration & cfg, const theta::Setting & s, int depth);
    
    boost::ptr_vector<theta::Distribution> distributions;
    std::map<theta::ParId, size_t> parid_to_index;
};

/** \brief A data source using a model
 *
 * Configured via a setting group like
 * \code
 * source = {
 *   type = "model_source";
 *   model = "@some-model-path";
 *   dice_poisson = false; // optional; default is true
 *   dice_template_uncertainties = false; // optional; default is true
 *   dice_rvobs = false; // optional; default is true
 *   override-parameter-distribution = "@some-dist"; // optional
 *   parameters-for-nll = { p1 = 0.0; p2 = 1.0; p3 = "diced_value"; }; //optional; assuming p1, p2, p3 are parameters
 *   rnd_gen = { seed = 123; }; // optional
 * };
 * \endcode
 *
 * \c type must be "model_source" to select this DataSource type
 *
 * \c model specifies the model to use for data creation
 *
 * \c dice_poisson controls whether the pseudo data this plugin provides dices a poisson around the model prediction in each bin
 *   or just returns the model prediction directly; see below for details.
 *
 * \c dice_template_uncertainties controls whether there will be a random smearing of the templates within their uncertainties prior to sampling pseudo data;
 *    see below for details.
 *
 * \c dice_rvobs controls whether or not to dice the real-valued observables.
 *
 * \c override-parameter-distribution is an optional setting overriding the distribution from which parameter values are drawn, see below. It
 *     has to provide exactly the same parameters as the model.
 *
 * \c parameters-for-nll defines the parameter values used for saving the value of the negative log-likelihood. They are
 *     either given directly as floating point, or the special string "diced_value". Each parameter of the model
 *     must be specified here. If not given, the negative log-likelihood will not be saved.
 *
 * \c rnd_gen defined the random number generator to use. The default is to construct a \link theta::Random Random \endlink class
 *    with a empty setting group which will use defaults documented there.
 *
 * For each call of fill
 * <ol>
 *   <li>parameter values are sampled from the parameter distribution of the model. If the setting \c override-parameter-distribution
 *       is given, it will be used for this step instead.</li>
 *   <li>the parameter values from the previous step are used to calculate the Model prediction. If \c dice_template_uncertainties is true,
 *      this will be done using the method HistogramFunction::getRandomFluctuation instead of the evaluation operator of \link theta::HistogramFunction HistogramFunction\endlink.
 *      What this means exactly depends on the HistogramFunction instances used; see their documentation of \c getRandomFlucutation method for details.</li>
 *   <li>if \c dice_poisson is \c true, a Poisson random number is drawn for each bin around the mean of the prediciton from the previous step. Otherwise,
 *      the bin contents from the model prediction are used as pseudo data directly (somtimes called "Asimov data").</li>
 *   <li>if \c parameters-for-nll is given, the likelihood function is built for the Data sampled
 *      in the previous step and evaluated at the given parameter values.</li>
 * </ol>
 *
 * \c model_source will create a column in the event table for each model parameter and save the values used to sample the pseudo
 * data. If \c parameters-for-nll is given, a column "nll" will be created where the value of the negative log-likelihood is saved.
 *
 * The \c parameters-for-nll setting can be used to contruct Feldman-Cousins intervals which uses the likelihood ratio as
 * ordering principle. One of the likelihood values uses the true value of the parameter of interest and is calculated here in \c model_source.
 * The other likelihood value can be calculated used the \link mle mle \endlink producer which saves the nll value at the
 * minimum of the negative log-likelihood.
 */
class model_source: public theta::DataSource, public theta::RandomConsumer {
public:

    /// Construct from a Configuration; required by the plugin system
    model_source(const theta::Configuration & cfg);

    /** \brief Fills the provided Data instance with data from the model
     */
    virtual void fill(theta::Data & dat);
    
private:
    theta::ParValues parameters_for_nll;
    
    bool save_nll;
    theta::Column c_nll;
    
    theta::ParIds par_ids;
    std::vector<theta::Column> parameter_columns;
    
    std::auto_ptr<theta::Model> model;
    std::auto_ptr<theta::Distribution> override_parameter_distribution;
    
    bool dice_poisson;
    bool dice_template_uncertainties;
    bool dice_rvobs;

};

#endif

