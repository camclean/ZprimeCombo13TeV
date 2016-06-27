#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "interface/decls.hpp"
#include "interface/utils.hpp"
#include <cstring>
#include <vector>

namespace theta{

/** \brief Container for a vector of doubles
 * 
 * Container for a vector of doubles, with some common operations, to be used
 * by different Histogram classes as internal storage for the data.
 * 
 * As optimization, no memory is allocated for a size of 0.
 */
class DoubleVector{
private:
    double * data;
    size_t n_data, n_alloc; // n_alloc was actually allocated, n_data is the (user-visible) number of elements. Like std::vector.size versus capacity
    
public:
    
    /// Create a new DoubleVector with \c n entries, all set to zero.
    explicit DoubleVector(size_t n=0);
    
    /// Copy constructor
    DoubleVector(const DoubleVector & rhs);
    
    /// Destructor
    ~DoubleVector();
    
    /// Assignment
    void operator=(const DoubleVector & rhs){
        if(&rhs == this) return;
        if(n_data != rhs.n_data){
            realloc(rhs.n_data);
        }
        utils::copy_fast(data, rhs.data, n_data);
    }
    
    /// reallocate (unitialized!) to new size n
    void realloc(size_t n);
    
    /// Set all values to the given value
    void set_all_values(double value){
        std::fill(data, data + n_data, value);
    }
    
    /// Multiply with a scalar
    void operator*=(double a){
        utils::mul_fast(data, a, n_data);
    }

    /// Add another DoubleVector. No size checking is done.
    void operator+=(const DoubleVector & rhs){
        utils::add_fast(data, rhs.data, n_data);
    }

    /// this += c*rhs. No range checking is done.
    void add_with_coeff(double c, const DoubleVector & rhs){
        utils::add_fast_with_coeff(data, rhs.data, c, n_data);
    }
    
    // calculate this += c1*v1 + c2*v2
    void add_with_coeff2(double c1, const DoubleVector & v1, double c2, const DoubleVector & v2){
    	utils::add_fast_with_coeff2(data, v1.data, c1, v2.data, c2, n_data);
    }

    /// The current number of stored values
    size_t size() const{
        return n_data;
    }
    
    /** \brief Returns the entry number i
     *
     * It must hold: i>=0 && i&lt;n; this function does not do any  range checking.
     */
    double get(size_t i) const{
        return data[i];
    }

    /** \brief Set entry i to data.
     *
     * It must hold: i>=0 && i&lt;n; this function does not do any  range checking.
     */
    void set(size_t i, double value){
       data[i] = value;
    }
    
    /// Get the sum of all data values.
    double get_sum() const{
        double result = 0.0;
        for(size_t i=0; i < n_data; ++i){
            result += data[i];
        }
        return result;
    }
    
    ///@{
    /** \brief Return a pointer to the raw double data
     * 
     * In case of size()==0, returns a NULL pointer.
     */
    const double* get_data() const{
        return data;
    }
    double* get_data() {
        return data;
    }
    ///@}
};
    

/** \brief A Histogram class holding binned, 1D data, without overflow and underflow bins.
 *
 * Only equally spaced binning is supported. Binning is therefore fully determined
 * by a range (\c xmin, \c xmax) and the number of bins.
 * 
 * Bins are indexed from 0 up to and including nbins-1; no data for "underflow"
 * or "overflow" is kept.
 */
class Histogram1D: public DoubleVector {
private:
    double xmin, xmax;
    
    // throws the actual exception, with the message.
    void fail_check_compatibility(const Histogram1D & h) const;
    
public:
    
    using DoubleVector::operator*=;

    /** \brief Create an empty Histogram with \c bins bins with range (\c xmin, \c xmax )
     */
    explicit Histogram1D(size_t bins=0, double xmin=0, double xmax=1);
    
    Histogram1D(double xmin_, double xmax_, const DoubleVector & dv): DoubleVector(dv), xmin(xmin_), xmax(xmax_){}
    
    /// Get the number of bins of this Histogram
    size_t get_nbins() const{
       return size();
    }

    /// Get the minimum x value for this Histogram
    double get_xmin() const{
       return xmin;
    }

    /// Get the maximum x value of this Histogram
    double get_xmax() const{
       return xmax;
    }

    /** \brief Add weight to the bin corresponding to xvalue
     * 
     * filling xvalues outside of the range will have no effect.
     */
    void fill(double xvalue, double weight);

    /** \brief Returns the x value for the bin center of bin ibin
     *
     *  The returned value is only valid for ibin>=1 && ibin&lt;=nbins
     */
    double get_bincenter(int ibin) const{
        return xmin + (ibin+0.5) * (xmax - xmin) / get_nbins();
    }

    /** \brief Multiply the Histogram \c other with this, bin-by-bin
     *
     * Throws an \c InvalidArgumentException if \c other is not compatible with this Histogram.
     * \c &other must be different from \c this.
     */
    void operator*=(const Histogram1D & other);

    /** \brief check compatibility of \c this to the \c other Histogram.
     *
     * A Histogram is considered compatible if it has the exact same number of bins and range. In case if incompatibility, an
     * \c InvalidArgumentException is thrown.
     */
    //note: this is split into checking and reporting in another routine fail_check_compatibility
    // to increase inlining probability
    void check_compatibility(const Histogram1D & h) const{
        if (get_nbins() != h.get_nbins() || xmin!=h.xmin || xmax != h.xmax){
           fail_check_compatibility(h);
        }
    }
    
    /// Set all values to 0.0
    // provided for compatibility with DataT
    void reset(){
        set_all_values(0.0);
    }
    
    /// Set to a zero Histogram with the given bin number and range
    void reset(size_t nbins, double xmin_, double xmax_){
        realloc(nbins);
        xmin = xmin_;
        xmax = xmax_;
        set_all_values(0.0);
    }
    
    /** \brief Calculate this = this * (nominator/denominator)^exponent, bin-by-bin.
     *
     * An \c InvalidArgumentException is thrown if either \c nominator or \c denominator are not compatible with this.
     */
    void multiply_with_ratio_exponented(const Histogram1D & nominator, const Histogram1D & denominator, double exponent);
};

}

#endif
