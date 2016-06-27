#ifndef HISTOGRAM_WITH_UNCERTAINTIES_HPP
#define HISTOGRAM_WITH_UNCERTAINTIES_HPP

#include "interface/decls.hpp"
#include "interface/histogram.hpp"
#include "interface/exception.hpp"
#include "interface/variables.hpp"

#include <cmath>
#include <vector>

namespace theta{

/** \brief A Histogram class holding binned, 1D data, without overflow and underflow bins, where each bin also has an uncertainty.
 *
 */
//This class implements similar methods as Histogram1D. However, it does not inherit from Histogram1D
//as the overall meaning is different and it has the risk of information truncation in the case of
// h *= hwu
// where h is a Histogram1D and hwu is a Histogram1DWithUncertainties.
//
//The histogram includes optimizations for the case where all uncertainties are 0, so that the overhead
// in this case compared to Histogram1D should be small.
class Histogram1DWithUncertainties {
private:
    
    double xmin, xmax;
    // sq_uncertainties contain the *squared* uncertainty in this bin (easier to calculate this way ...)
    DoubleVector values, sq_uncertainties;
    
    bool nontrivial_unc;
    
    // set nontrivial_unc to true and allocate sq_uncertainties (if necessary).
    void set_nontrivial_unc();
    
    void fail_check_compatibility(const Histogram1DWithUncertainties & h) const;
public:

    /// create a Histogram with the given range and number of bins
    explicit Histogram1DWithUncertainties(size_t bins=0, double xmin=0, double xmax=1);
    
    
    /// Construct from Histogram1D, setting all uncertainties to zero
    explicit Histogram1DWithUncertainties(const Histogram1D & h);
    
    void reset(size_t nbins, double xmin_, double xmax_){
        xmin = xmin_;
        xmax = xmax_;
        values.realloc(nbins);
        values.set_all_values(0.0);
        sq_uncertainties.realloc(0);
        nontrivial_unc = false;
    }
    
    //@{
    /// Get metadata
    
    /// Get the number of bins of this Histogram
    size_t get_nbins() const{
       return values.size();
    }

    /// Get the minimum x value for this Histogram
    double get_xmin() const{
       return xmin;
    }

    /// Get the maximum x value of this Histogram
    double get_xmax() const{
       return xmax;
    }
    //@}
    
    
    //@{
    /// Get and set data.
    double get_value(size_t i) const{
        return values.get(i);
    }
    
    // alias for get_value:
    double get(size_t i) const{
        return values.get(i);
    }
    
    double get_uncertainty(size_t i) const{
        if(nontrivial_unc){
            return std::sqrt(sq_uncertainties.get(i));
        }
        else return 0.0;
    }
    
    double get_uncertainty2(size_t i) const{
        if(nontrivial_unc){
            return sq_uncertainties.get(i);
        }
        else return 0.0;
    }
    
    double * get_data(){
        return values.get_data();
    }
    
    const double * get_data() const{
        return values.get_data();
    }

    const DoubleVector & get_values() const{
        return values;
    }
    
    Histogram1D get_values_histogram() const{
        return Histogram1D(xmin, xmax, values);
    }
    
    Histogram1D get_uncertainty2_histogram() const{
        if(nontrivial_unc){
            return Histogram1D(xmin, xmax, sq_uncertainties);
        }
        else{
            return Histogram1D(get_nbins(), xmin, xmax);
        }
    }

    // setting uncertainty to NAN leaves it unchanged.
    void set(size_t i, double value, double uncertainty = NAN){
        values.set(i, value);
        if(!std::isnan(uncertainty)){
            if(nontrivial_unc || uncertainty!=0){
                set_nontrivial_unc();
                sq_uncertainties.set(i, uncertainty * uncertainty);
            }
        }
    }

    void set_unc2(size_t i, double value, double uncertainty2 = NAN){
        values.set(i, value);
        if(!std::isnan(uncertainty2)){
            if(nontrivial_unc || uncertainty2!=0){
                set_nontrivial_unc();
                sq_uncertainties.set(i, uncertainty2);
            }
        }
    }

    // set values, uncertainties are set to zero.
    void set(const Histogram1D & rhs){
        xmin = rhs.get_xmin();
        xmax = rhs.get_xmax();
        values = rhs;
        sq_uncertainties = DoubleVector();
        nontrivial_unc = false;
    }

    void set(const Histogram1D & values_, const Histogram1D & uncertainties){
    	values_.check_compatibility(uncertainties);
        xmin = values_.get_xmin();
        xmax = values_.get_xmax();
        values = values_;
        sq_uncertainties = uncertainties;
        for(size_t i=0; i<sq_uncertainties.size(); ++i){
            double newval = sq_uncertainties.get(i);
            sq_uncertainties.set(i, newval * newval);
        }
        nontrivial_unc = true;
    }

    //@}

    //@{
    /** \brief Arithmetic Operations
     * 
     * Operations are provided for both cases: with and without uncertainties. In the "no uncertainty" case,
     * uncertainties of 0 are assumed.
     */
    void operator*=(double a){
        values *= a;
        if(nontrivial_unc){
            sq_uncertainties *= a * a;
        }
    }
    
    void operator+=(const Histogram1DWithUncertainties & other){
        values += other.values;
        if(other.nontrivial_unc){
            set_nontrivial_unc();
            sq_uncertainties += other.sq_uncertainties;
        }
    }
    
    void add_with_coeff(double c, const Histogram1DWithUncertainties & other){
        values.add_with_coeff(c, other.values);
        if(other.nontrivial_unc){
            set_nontrivial_unc();
            sq_uncertainties.add_with_coeff(c * c, other.sq_uncertainties);
        }
    }
    
    void operator+=(const Histogram1D & other){
        values += other;
    }
    
    void add_with_coeff(double c, const Histogram1D & other){
        values.add_with_coeff(c, other);
    }
    //@}
    
    /** \brief Add with coefficient, but only for the value, not the uncertainties
     * 
     * Calculate this += c * rhs, but only changing the values of this, and only reading the values of rhs, not
     * the uncertainties.
     */
    void add_with_coeff_values(double c, const Histogram1DWithUncertainties & rhs){
        values.add_with_coeff(c, rhs.values);
    }
    
    /** \brief Add the squared uncertainties
     * 
     * calculate this->unc2 += c * rhs.unc2  where "unc2" is the histogram of squarted uncertainties
     *
     */
    void add_with_coeff_unc2(double c, const Histogram1DWithUncertainties & rhs){
        if(rhs.nontrivial_unc){
            set_nontrivial_unc();
            sq_uncertainties.add_with_coeff(c, rhs.sq_uncertainties);
        }
    }

    /** \brief check compatibility of \c this to the \c other Histogram.
     *
     * A Histogram is considered compatible if it has the exact same number of bins and range. In case if incompatibility, an
     * \c InvalidArgumentException is thrown.
     */
    // see Histogram1D
    void check_compatibility(const Histogram1DWithUncertainties & h) const{
        if (get_nbins() != h.get_nbins() || xmin!=h.xmin || xmax != h.xmax){
           fail_check_compatibility(h);
        }
    }
};

}


#endif
