#ifndef HISTOGRAM_FUNCTION_HPP
#define HISTOGRAM_FUNCTION_HPP

#include "interface/decls.hpp"
#include "interface/histogram.hpp"
#include "interface/histogram-with-uncertainties.hpp"
#include "interface/variables.hpp"

namespace theta {
    
    /** \brief A Histogram-valued function which depends on zero or more parameters.
     *
     * This class is used extensively for model building: a physical model is given by specifying
     * the expected observation in one or more observables and this expectation in turn is specified
     * by histograms which depend on the model parameters. As this can be seen as a histogram-valuesd function,
     * the class is called \c HistogramFunction.
     */
    class HistogramFunction{
    public:
        
        /// Define us as the base_type for derived classes; required for the plugin system
        typedef HistogramFunction base_type;
        
        //@{
        /** \brief Calculate h += coeff * hf(values)
         *
         * This somewhat complicated construction is done for reasons of efficiency and avoiding allocations.
         * 
         * h must be initialized with the correct range and binning, according to get_histogram_dimensions.
         */
        virtual void add_with_coeff_to(Histogram1DWithUncertainties & h, double coeff, const ParValues & values) const = 0;
        virtual void add_with_coeff_to(Histogram1D & h, double coeff, const ParValues & values) const = 0;
        //@}
        
        /*
         * Calculate
         * result = hf(values)   and
         * derivatives[p] += coeff * d / dp  hf(values)    for all parameters p in get_parameters()
         * 
         * All Histograms "derivatives[p]" for p in get_parameters() MUST be initialized with the correct range and binning; result does not have to be.
         * 
         * The version with uncertainties adds the derivative of the *squared* uncertainties to the uncertainties of derivaties[p], i.e. it calculates the above and
         * for the uncertainties, it does
         *   derivatives[p].unc2 += coeff^2 * d / dp  hf(values).unc2
         * where ".unc2" refers to the histogram of squared uncertainties.
         * 
         * [Note: for many HistogramFunctions, d / dp hf(values).unc2 is identical 0.]
         */
        virtual void eval_and_add_derivatives(Histogram1D & result, std::map<ParId, Histogram1D> & derivatives, double coeff, const ParValues & values) const;
        virtual void eval_and_add_derivatives(Histogram1DWithUncertainties & result, std::map<ParId, Histogram1DWithUncertainties> & derivatives, double coeff, const ParValues & values) const;

        /** \brief Returns the parameters which this HistogramFunction depends on.
         */
        const ParIds & get_parameters() const{
            return par_ids;
        }

        /** \brief Get the dimensions of the Histogram (nbins, xmin, xmax) filled by the evaluation operators
         *
         * This function is used as part of the consistency checks to make sure that the Histogram dimensions match; to save
         * time, it is not usually not used during usual likelihood evaluation, etc.
         */
        virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const = 0;

        /// Declare the destructor virtual as there will be polymorphic access to derived classes
        virtual ~HistogramFunction(){}
        
    protected:
        /// To be filled by derived classes:
        ParIds par_ids;
    };
    

    /** \brief A simple HistogramFunction which always returns the same Histogram, independent of any parameters.
     */
    class ConstantHistogramFunction: public HistogramFunction{
    public:

        virtual void add_with_coeff_to(Histogram1DWithUncertainties & h, double coeff, const ParValues & values) const;
        virtual void add_with_coeff_to(Histogram1D & h, double coeff, const ParValues & values) const;
        virtual void eval_and_add_derivatives(Histogram1D & result, std::map<ParId, Histogram1D> & derivatives, double coeff, const ParValues & values) const;
        virtual void eval_and_add_derivatives(Histogram1DWithUncertainties & result, std::map<ParId, Histogram1DWithUncertainties> & derivatives, double coeff, const ParValues & values) const;
        virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const;

    protected:
        /** \brief Set the constant Histogram to return
         *
         * This method is meant for derived classes which can use it to set the constant Histogram to
         * be returned by operator()
         */
        void set_histo(const Histogram1DWithUncertainties & h);
        
        /** \brief Default constructor to be used by derived classes
         */
        ConstantHistogramFunction();
        
     private:
        Histogram1DWithUncertainties h_wu;
        Histogram1D h;
    };
 
    /** \brief Build a HistogramFunction according to the given configuration and return the result
     *
     * This assumes that the HistogramFunction specified by cfg does not depend on any parameters.
     * If this is not the case, an invalid_argument exception will be thrown.
     */
    Histogram1DWithUncertainties get_constant_histogram(const Configuration & cfg);
    
    
}

#endif
