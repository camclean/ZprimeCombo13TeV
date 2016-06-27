#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"

#include <map>

namespace theta{

    /** \brief A probability distribution of real-values random variables in one or more dimensions
     *
     * Implementations of this class provide methods for generating random numbers according to the distributions
     * they are representing. Further, it provides the negative logarithm of the probability density, including
     * derivatives at a given point.
     *
     * The intended use of this class is twofold:
     * <ol>
     *   <li>As priors in the pseudo data generation. Here, only the sample method is important.</li>
     *   <li>As additional terms in the likelihood function. Here, the other methods
     *       (mainly mpv, width, evalNL, support) are used.</li>
     * </ol>
     *
     * Each Distribution is defined for a set of real-valued random variables (implemented as ParIds); these
     * are the parameters as returned by Distribution::getParameters().
     *
     * The Distribution itself can depend on parameters, e.g., a Gaussian distribution can depend on the mean;
     * these are the parameters returned by Distribution::getDistributionParameters(). If using a Distribution
     * which depends on such parameters, the mathematical distribution is only completely defined once the
     * values for these parameters are specified and these values must be set in the arguments
     * to sample, mode, evalNL, etc.
     */
    class Distribution{
    public:
        /// Define us as the base_type for derived classes; required for the plugin system
        typedef Distribution base_type;
        
        /** Sample values from this distribution using \c rnd as random number generator
         *  and respecting limits set for the parameters in \c vm.
         *
         * The ParValues instance \c result will be used to set the random values.
         * This will set values for all variables this distribution
         * is defined on (i.e. all returned by getVariables()). Other parameter
         * values will not be touched.
         *
         * \param[out] result Fill the sampled values here.
         * \param rnd Proxy to the random number generator to use for sampling.
         */
        virtual void sample(ParValues & result, Random & rnd) const = 0;

        /** \brief Provides the mode (most probable values)
         *
         * All parameters returned by getParameters() are set to their most
         * probable value.
         *
         * Derived classes must ensure that if calling evalNL with all parameter values
         * set to their mode as returned by this function has a non-zero probability
         * (i.e., a non-infinite evalNL result).
         *
         * This function is mainly used to select valid start values for
         * algorithms like minimizations or markov chains.
         */
        virtual void mode(ParValues & result) const = 0;

        /** \brief The negative logarithm of the probability density.
         * 
         * The density is not guaranteed to be normalized, i.e., the negative
         * logarithm can be shifted by an arbitrary (but constant) value.
         *
         * \c values must contain (at least) the variable ids in getParameters().
         * Otherwise, a NotFoundException is thrown.
         * 
         * \param values The point in parameter space the density is to be evaluated.
         * \return The density at \c values.
         */
        virtual double eval_nl(const ParValues & values) const = 0;
        
        virtual double eval_nl_with_derivative(const ParValues & values, ParValues & derivative) const;

        /** \brief Get the support of a parameter
         *
         * The support is the set of values on which the density is non-vanishing. If the support is not
         * an interval, this method should return the smallest interval containing the support.
         *
         * This is mainly used to set constraints for that parameter in a minimization procedure.
         *
         * If \c p is not in getParameters(), the behaviour is undefined (i.e., derived one-dimensional classes
         *  need not check whether p is the correct ParId).
         */
        virtual const std::pair<double, double> & support(const ParId & p) const = 0;

        /** \brief Get the random variables of this Distribution
         */
        const ParIds & get_parameters() const{
            return par_ids;
        }
        
        /** \brief Get the parameters of this Distribution
         */
        const ParIds & get_distribution_parameters() const{
            return distribution_par_ids;
        }
        
        /// declare destructor as virtual, as polymorphic access will happen
        virtual ~Distribution();
        
    protected:
        ParIds par_ids;
        ParIds distribution_par_ids;
    };
    
    class EmptyDistribution: public Distribution{
    public:
        virtual void sample(ParValues & result, Random & rnd) const{}
        virtual void mode(ParValues & result) const{}
        virtual const std::pair<double, double> & support(const ParId & p) const;
        virtual double eval_nl(const ParValues & values) const { return 0.0;}
        virtual double eval_nl_with_derivative(const ParValues & values, ParValues & derivative) const{return 0.0; }
        virtual ~EmptyDistribution();
    };


    /** \brief Class saving parameter ranges
     *
     */
    class Ranges{
       std::map<ParId, std::pair<double, double> > ranges;
    public:
       Ranges(){}
       Ranges(const std::map<ParId, std::pair<double, double> > & ranges);

       explicit Ranges(const theta::Distribution & dist);

       /// set ranges from the support of the distribution
       void set_from(const theta::Distribution & dist);

       /// set the range
       void set(const ParId & pid, const std::pair<double, double> & range);

       /// Get the range for the specified parameter. Throws invalid_argument in case no range is set for this parameter
       const std::pair<double, double> & get(const ParId & pid) const;

       /// Returns whether or not the parameter is fixed, i.e. whether lower and upper bound coincide.
       bool fixed(const ParId & pid) const;

       /// truncate the values to the ranges. Parameters for which no range is set are not changed.
       void trunc(ParValues & values) const;
       
    };
}

#endif

