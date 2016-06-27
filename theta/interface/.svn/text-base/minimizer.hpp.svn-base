#ifndef MINIMIZER_HPP
#define MINIMIZER_HPP

#include "interface/decls.hpp"
#include "interface/distribution.hpp"
#include "interface/variables.hpp"
#include "interface/matrix.hpp"

#include <boost/shared_ptr.hpp>

namespace theta{

    /** \brief The result of a minimization process, returned by \c Minimizer::minimize().
     */
    struct MinimizationResult{
        
        /** \brief The function value at the minimum.
         */
        double fval;

        /** \brief The parameter values at the function minimum.
         *
         * Contains exactly the parameters the function to minimize depends
         * on.
         */
        ParValues values;

        /** \brief The errors at the function minimum, in positive direction.
         *
         * Contains all parameters the function to minimize depends
         * on. How this is calculated depends on the minimizer used. In cases where
         * the minimizer provides symmetrical errors, the entries are equal to \c errors_minus.
         *
         * Empty if not provided by the minimizer.
         */
        ParValues errors_plus;

        /** \brief The errors at the function minimum, in negative direction.
         *
         * Contains all parameters the function to minimize depends
         * on. How this is calculated depends on the minimizer used. In cases where
         * the minimizer provides symmetrical errors, the entries are equal to \c errors_plus.
         *
         * Note that while these are the errors in negative direction, the
         * entries are always positive in case it contains valid errors.
         *
         * Empty if not provided by the minimizer.
         */
        ParValues errors_minus;

        /** \brief Contains the error matrix at the minimum.
         *
         * It is quadratic and has values.size() rows.
         * The index convention is such that 0 corresponds to the first ParId in
         * the sorted (=as iterated) list of parameter ids contained in \c values.
         *
         * It is set to a Matrix of size 0 in case the
         * minimization does not provide errors.
         */
        Matrix covariance;
        
        /// Define explicitely as ParValues::operator= is private
        void operator=(const MinimizationResult& rhs);
    };
    
    
    /** \brief Class providing information useful for minimization for a function to minimize.
     * 
     * It stores at least the start step and range information for the parameters. However, Minimizers
     * can add additional information based on their requirements.
     */
    class FunctionInfo{
    public:
        /// The start values for the minimization.
        const ParValues & get_start() const{
            return start;
        }
        
        /// The parameter ranges. Guaranteed to contain all parameters of f.
        const Ranges & get_ranges() const{
            return ranges;
        }
        
        /// the step values contain approximately the standard deviation of the parameters to estimate
        const ParValues & get_step() const{
            return step;
        }
        
        /// The fixed parameters
        const ParIds & get_fixed_parameters() const{
            return fixed_parids;
        }
        
        /// Declare destructor virtual as we allow deriving; make it pure to avoid creating a FunctionInfo directly.
        virtual ~FunctionInfo() = 0;
        
    private:
        ParValues start;
        ParValues step;
        Ranges ranges;
        ParIds fixed_parids;
        
    protected:
        // The value in fixed_parameters overrides the value in start
        FunctionInfo(const ParValues & start, const ParValues & step, const Ranges & ranges, const ParValues & fixed_parameters);
    };


    /** \brief Abstract interface to different minimizers.
     *
     * The possible settings are documented at derived classes. The interface provides two methods:
     *  - the simpler "minimize" method just receives the Function to minimize, with values for the start values, step size, and parameter ranges.
     *  - the more complicated interface using FunctionInfo, which is (currently) only useful for likelihood functions. A clarrer would:
     *     * call create_nll_function_info first
     *     * call minimize2 with the negative log-likelihood and the FunctionInfo from the previous step
     * The second method the minimizer to save arbitrary data derived from the model about the nll function which can be used to optimize
     * the minimization. For example, a minimizer could calculate the covariance matrix of the Asimov data and store it in the FunctionInfo class.
     * As the covariance matrix usually does not have a strong dependence on the data, this provides a good estimate for the covariance for subsequent
     * minimizations of a likelihood function derived from the same model.
     * 
     * For backward compatibility, the minimize2 is implemented generically for all minimizers here which just call minimize.
     */
    class Minimizer{
    private:
        // to detect interface violations, make our own version of FunctionInfo here; create_nll_info actually
        // returns a DefFunctionInfo and minimize2 checks for that.
        class DefFunctionInfo: public FunctionInfo{
        public:
            DefFunctionInfo(const ParValues & start_, const ParValues & step_, const Ranges & ranges_, const ParValues & fixed_parameters): FunctionInfo(start_, step_, ranges_, fixed_parameters){ }
            ~DefFunctionInfo();
        };
        
    public:
        
        /// Define this class as the base_type for derived classes; required for the plugin system
        typedef Minimizer base_type;

        /// declare destructor virtual as we expect polymorphic access to derived classes
        virtual ~Minimizer();

        /** \brief Find the minimum of a function.
         *
         * The minimum of f is found, starting at the specified start values, step sizes and
         * ranges. The step size should be an estimate of the parameter uncertainty.
         *
         * If the minimization fails,a MinimizationException is thrown. The reasons for such a failure 
         * depend on the particular minimization algorithm and should be documented in derived classes.
         * Typical failues include not being able to reach the desired accuracy in the maximum number
         * of iterations.
         * 
         * If for a parameter either step is 0.0 or the range contains only one value, this parameter is considered
         * as constant and it is not varied during minimization.
         */
        virtual MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                const theta::ParValues & step, const Ranges & ranges) = 0;
        
        /** \brief Alternative minimize method with more information
         *
         * Similar to minimize, but can make use of more information about the minimization problem
         * via the FunctionInfo instance.
         *
         * \c fixed_parameters specify all parameters that should be fixed in the minimization. This has the same effect
         * as using a step size of 0, and setting start and according (trivial) range in "minimize".
         * The parameters specified in \c fixed_parameters must be a subset of info.get_fixed_parameters().
         * 
         * \c info must be an instance created previously with create_nll_function_info from the same Minimizer.
         * Otherwise, the behavior is undefined.
         * 
         * Derived classes should implement this method if they want to make use of the FunctionInfo mechanism. The default
         * implementation just calls "minimize" using \c info for all parameters to minimize.
         */
        virtual MinimizationResult minimize2(const Function & f, const FunctionInfo & info, const ParValues & fixed_parameters = ParValues());
        
        
        /** \brief Create a FunctionInfo instance suitable for minimizing negative log-likelihood functions from the given Model
         * 
         * \c fixed_parameters are values for fixing parameters for the negative log-likelihood function; this overrides the parameter distribution
         *   from the model / from \c override_parameter_distribution.
         * 
         * \c override_parameter_distribution can be a null pointer in which case the parameter distrbution from the Model is used.
         * 
         * Important: the set of fixed parameters for which FunctionInfo instance is created here and used for minimize2 must be the same (only the
         *  value of those fixed parameters is allowed to change).
         * 
         * The default implementation constructs a FunctionInfo instance with start and ranges from the parameter Distribution;
         * step is derived from asimov_likelihood_widths.
         */
        virtual boost::shared_ptr<FunctionInfo> create_nll_function_info(const Model & m, const boost::shared_ptr<Distribution> & override_parameter_distribution,
                                                                         const ParValues & fixed_parameters = ParValues());
    };
    
}

#endif
