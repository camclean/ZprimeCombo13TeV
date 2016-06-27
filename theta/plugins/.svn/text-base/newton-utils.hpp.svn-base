#ifndef PLUGINS_NEWTON_UTILS_HPP
#define PLUGINS_NEWTON_UTILS_HPP

// This header defines some utility function used in the newton_minimizer.
// They could be kept internal to newton.cpp/hpp, but exposing them here enables
// easier testing and re-use. Still, being somewhat internal to the newton_minimizer, the classes
// are in the newton_internal namespace for now.

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/matrix.hpp"

#include <vector>
#include <algorithm>
#include <cmath>

#include <boost/function.hpp>


namespace newton_internal {
    

/** \brief Abstract base class for functions R^n -> R
 *
 * This class provides   
 *  * a simple vector&lt;double&gt; interface for evaluating function values and derivatives
 *  * range constraints for the parameters
 *  * numerical derivative (which respects the ranges)
 * 
 * Derived classes have to:
 *  * implement the evaluation operator() and the number of dimensions ndim()
 *  * initialize all protected members correctly.
 */
class RangedFunction{
protected:
    std::vector<double> step; // the typical scale the function changes on. Used in the numerical derivative
    std::vector<double> range_min;
    std::vector<double> range_max;
    double epsilon_f;
    
public:
    
    /// Evaluate the function at x; x.size()==ndim() must hold, and x must be within the range (i.e. trunc_to_range(x)==0 must hold).
    virtual double operator()(const std::vector<double> & x) const = 0;
    
    /// The dimensionality
    virtual size_t ndim() const = 0;
    
    /// Evaluate the function and its derivatives. ipar is -1 to calculate all derivatives, or the index to x0 / grad
    virtual double eval_with_derivative(const std::vector<double> & x0, std::vector<double> & grad) const;
    
    /// Truncate x to the range this function is defined on. Returns the number of entries changed in x.
    size_t trunc_to_range(std::vector<double> & x) const;
    
    /// Mark the destructor virtual, as this is an abstract base class.
    virtual ~RangedFunction();
};

/** \brief Wrapper for a theta::Function as RangedFunction
 * 
 * This is used to "convert" a theta::Function instance as a RangedFunction object.
 * Only the non-fixed parameters of the theta::Function are exposed.
 */
class RangedThetaFunction: public RangedFunction{
private:
    const theta::Function & f;
    theta::ParIds fixed_parameters, nonfixed_parameters;
    bool use_f_derivative;
    
    std::vector<theta::ParId> nonfixed_pids;
    
    // pvs for *all* parameters of f
    mutable theta::ParValues pvs, pv_der;
public:
    
    /** \brief Construct from a theta::Function and ranges
     * 
     * The range and steps of the RangedFunction are set according to step and ranges. Any parameter given in \c fixed_values
     * will make this parameter fixed; the value given in \c fixed_values overrides the settings via step and ranges.
     */
    RangedThetaFunction(const theta::Function & f, const theta::ParValues & fixed_values, const theta::ParValues & step, const theta::Ranges & ranges, bool use_f_derivative);
    
    /// Set the value for epsilon_f explicitly.
    void set_epsilon_f(double newval);
    
    /// Get the subset of the theta::Function parameters that are fixed
    const theta::ParIds & get_fixed_parameters() const;
    
    /// Get the subset of the theta::Function parameters that are not fixed
    const theta::ParIds & get_nonfixed_parameters() const;
    
    /// Get the number of (non-fixed) parameters
    virtual size_t ndim() const;
    
    virtual double eval_with_derivative(const std::vector<double> & x0, std::vector<double> & grad) const;
    
    // evaluate the theta::Function, as a function of the non-fixed parameters.
    virtual double operator()(const std::vector<double> & vals) const;
};


/** \brief Determine the (machine/rounding) accuracy for a RangedFunction.
 * 
 * f is evaluated around x0 -- changing only x0[i] -- and returns the smallest non-zero absolute difference in function
 * values observed. This is an estimate of the machine/rounding accuracy for f.
 * Note that while the return value is not very accurate (the points are tested increasing/decresing their distance by factor of 2),
 * this is usually good enough for the purpose of calculating rounding errors for f.
 * 
 * f(x) has to change if changing x[i]; otherwise, an invalid_argument exception is thrown.
 *
 * Note that the actual algorithm slightly differs from the one outlied above:
 * For special f and x0, such as f(x) = x**2 and x0 = 0, the algorithm above would return extremely
 * small epsilon which are not really valid at any other point but x0 = 0. Therefore, not actually f(x) and
 * f(x + d) are compared, but f_scale + |f(x)| and f_scale + |f(x + d)|. f_scale should be the (positive) "typical scale"
 * of (changes in) f in the region one wants to apply the resulting f_accuracy in; the default of 1.0 makes sense
 * e.g. for log-likelihood functions.
 */
double f_accuracy(const RangedFunction & f, const std::vector<double> & x0, size_t i, double f_scale = 1.0);


/** \brief Simple structure holding status of minimum search in an interval
 * 
 * If a minimum is sought with an inteerval search routine, the current status
 * is summarized as the working interval [a,b] in which the minimum is contained
 * and a point c in this interval with f(c) <= f(a) and f(c) <= f(b) is a cancidate
 * for the position of the minimum.
 */
struct min_triplet{
    double a, b, c; // a < c < b
    double fa, fb, fc; // fa >= fc <= fb
};

/// Definition for find_argmin: The one-dimensional function double -> double
typedef boost::function<double (double)> f1d;

/// Definition for find_argmin: The proposal-generating function
typedef boost::function<double (const min_triplet&)> fproposal;

/// Definition for find_argmin: The stopping criterion. Should return true if iteration should be stopped.
typedef boost::function<bool (const min_triplet&)> fstop;

/** \brief Simple proposal function for find_argmin
 *
 * Always proposes the midpoint of the larger subinterval of (a,c,b) in the triplet tr.
 */
double quadratic_interpolation(const min_triplet & tr);

/** \brief Quadratic interpolation function for find_argmin
 *
 * Interpolates the points ith a quadratic function and returns the position of the minimum.
 */
double bisect_larger(const min_triplet & tr);

/** \brief Find the position of the minimum of a 1D function.
 * 
 * The initial position of the minimum and search interval is passed via tr, which is also used
 * to communicate the result of the search.
 * 
 * The algorithm is as follows:
 * 
 * 1. stop(tr) is evaluated. If it returns true, the iteration is stopped.
 * 2. fp(tr) is called to generate a proposal \c p. If \c p is is not a finite number in the interval (a,b), p = bisect_larger(tr) it used
 *    as a fallback. \c p is used to update the current triplet tr.
 * 3. There is one pathological case requiring special attention: if the proposal p or the point c is very close to a,b
 *    then the interval (a,b) in the triplet hardly shrinks in one iteration. If this is the case, an additional
 *    proposal is generated using bisect_larger on the new triplet.
 * 
 * \c fp is usually set to \c quadratic_interpolation.
 * 
 * \c maxit is the number of maximum iterations (where one iteration consists of executing 1. to 3. as described above).
 *    If the minimum is not found within this number of iterations, the calculation aborts with a theta::Exception.
 */
void find_argmin(const f1d & f, min_triplet & tr, const fproposal & fp, const fstop & stop, unsigned int maxit = 10000);


/** \brief Various arithmetic functions for vectors and matrices.
 * 
 * As the newton_minimizer internally uses vector&lt;double&gt; and theta::Matrix, we provide some
 * utility functions for these types here:
 * \code
 *   vector<double> x, y;
 *   Matrix A;
 *   double c;
 *   ...
 *   add_with_coeff(x, y, c); // calculates x += c*y
 *   mul(x, A, y); // x = A*y
 *   mul(x, c); // x*= c
 *   double res = norm(x); // euclidean norm: res = sqrt(sum_i x_i**2)
 *   res = pnorm(x, y); // a pseudo-norm: res = max_i fabs(x[i] / y[i]); the y[i] are non-zero step sizes used to scale down x before taking the infinity-norm
 * \endcode
 * 
 * Note that the methods having a vector as result do not return the vector as return value; rather, they manipulate
 * the first argument. This allows avoiding allocations by the client re-using those result vectors.
 */
//@{
inline void add_with_coeff(std::vector<double> & lhs, const std::vector<double> & rhs, double coeff){
    theta_assert(lhs.size() == rhs.size());
    for(size_t i=0; i<rhs.size(); ++i){
        lhs[i] += coeff * rhs[i];
    }
}

inline void mul(std::vector<double> & l, const theta::Matrix & m, const std::vector<double> & x){
    size_t n_rows = m.get_n_rows();
    size_t n_cols = m.get_n_cols();
    theta_assert(x.size() == n_cols);
    l.resize(n_rows);
    for(size_t i=0; i<n_rows; ++i){
        l[i] = 0.0;
        for(size_t j=0; j<n_cols; ++j){
            l[i] += m(i,j) * x[j];
        }
    }
}

inline void mul(std::vector<double> & l, double c){
    const size_t n = l.size();
    for(size_t i=0; i<n; ++i){
        l[i]*=c;
    }
}

inline double pnorm(const std::vector<double> & x, const std::vector<double> & step){
    const size_t n = x.size();
    theta_assert(step.size()==n);
    double result = 0.0;
    for(size_t i=0; i<n; ++i){
        result = std::max(result, fabs(x[i] / step[i]));
    }
    return result;
}

inline double norm(const std::vector<double> & x){
    const size_t n = x.size();
    double s2 = 0.0;
    for(size_t i=0; i<n; ++i){
        s2 += x[i] * x[i];
    }
    return sqrt(s2);
}
//@}

}

#endif
