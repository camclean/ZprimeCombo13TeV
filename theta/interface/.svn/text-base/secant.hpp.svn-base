#ifndef PLUGINS_SECANT_HPP
#define PLUGINS_SECANT_HPP

#include "interface/exception.hpp"
#include <cmath>

/** \brief The secant method to find the root of a one-dimensional function
 *
 * \param x_low The lower end of the start interval
 * \param x_high The higher end of the start interval
 * \param x_accuracy If the found interval is shorter than this, the iteration will stop
 * \param f_x_low is function(x_low).Used to save one function evalutation.
 * \param f_x_high is function(x_high). Used to save one function evaluation.
 * \param f_accuracy If the absolute function value is lower than this, the iteration will stop
 * \param function The function object to use.
 *
 * The iteration will stop if either one if the x_accuracy and f_accuracy criteria is fulfilled.
 * To use only one criterion, set the other value to 0.0 to prevent it from being fulfilled.
 *
 * Note that the function values at x_low and x_high must have different sign. Otherwise,
 * an std::invalid_argument will be thrown.
 * All x and function values must be finite.
 */
template<typename T>
double secant(double x_low, double x_high, double x_accuracy, double f_x_low, double f_x_high, double f_accuracy, const T & function, int limit_depth = 500){
    theta_assert(std::isfinite(x_low) && std::isfinite(x_high));
    theta_assert(x_low <= x_high);
    theta_assert(std::isfinite(f_x_low) && std::isfinite(f_x_high));
    theta_assert(limit_depth > 0);
    if(f_x_low * f_x_high >= 0) throw std::invalid_argument("secant: function values have the same sign!");
    if(std::fabs(f_x_low) <= f_accuracy) return x_low;
    if(std::fabs(f_x_high) <= f_accuracy) return x_high;

    const double old_interval_length = x_high - x_low;    
    //calculate intersection point for secant method:
    double x_intersect = x_low - (x_high - x_low) / (f_x_high - f_x_low) * f_x_low;
    theta_assert(x_intersect >= x_low);
    theta_assert(x_intersect <= x_high);
    if(old_interval_length < x_accuracy){
        return x_intersect;
    }
    double f_x_intersect = function(x_intersect);
    double f_mult = f_x_low * f_x_intersect;
    //fall back to bisection if the new interval would not be much smaller:
    double new_interval_length = f_mult < 0 ? x_intersect - x_low : x_high - x_intersect;
    if(new_interval_length > 0.5 * old_interval_length){
        x_intersect = 0.5*(x_low + x_high);
        f_x_intersect = function(x_intersect);
        f_mult = f_x_low * f_x_intersect;
    }
    if(f_mult < 0){
        return secant(x_low, x_intersect, x_accuracy, f_x_low, f_x_intersect, f_accuracy, function, limit_depth - 1);
    }
    else if(f_mult > 0.0){
        return secant(x_intersect, x_high, x_accuracy, f_x_intersect, f_x_high, f_accuracy, function, limit_depth - 1);
    }
    //it can actually happen that we have 0.0. In this case, return the x value for
    // the smallest absolute function value:
    else{
        f_x_intersect = fabs(f_x_intersect);
        f_x_low = fabs(f_x_low);
        f_x_high = fabs(f_x_high);
        if(f_x_low < f_x_high && f_x_low < f_x_intersect) return x_low;
        if(f_x_high < f_x_intersect) return x_high;
        return x_intersect;
    }
}

// use brent's algorithm for root finding of the 1D function f.
// f_x_low * f_x_high < 0.0 or f_x_low==0.0 or f_x_high==0.0 must hold.
//
// f_accuracy is the estimated (machine / other) accuracy to which f(x) can be evaluated. Any
// absolute value of f equal to or below that -- or any function value difference below that -- will be considered as zero.
// Setting to 0.0 should work for well-behaved functions ...
template<typename FT>
double brent(const FT & f, double a, double b, double xtol, double fa, double fb, double f_accuracy = 1e-7, unsigned int limit = 10000){
    if(fabs(fa) <= f_accuracy) return a;
    if(fabs(fb) <= f_accuracy) return b;
    theta_assert(fa * fb < 0.0);
    if(fabs(fa) < fabs(fb)){
        std::swap(fa, fb);
        std::swap(a, b);
    }
    double d = 0.0;
    double c = a;
    double fc = fa;
    bool mflag = true;
    for(unsigned int i=0; i < limit; ++i){
        if(fabs(a-b) <= xtol){
            return (a + b) / 2;
        }
        if(fabs(fa) <= f_accuracy) return a;
        if(fabs(fb) <= f_accuracy) return b;
        double s;
        if(fabs(fa - fc) < f_accuracy && fabs(fb - fc) < f_accuracy){
            s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
        }
        else{
            s = b - fb * (b - a) / (fb - fa);
        }
        double tmp2 = (3 * a + b) / 4;
        if (!((s > tmp2 && s < b) || (s < tmp2 && s > b))
             || (mflag && fabs(s - b) >= fabs(b - c) / 2)
             || (!mflag && fabs(s - b) >= fabs(c - d) / 2)
             || (mflag && fabs(b - c) < xtol)
             || (!mflag && (fabs(c - d) < xtol))
           ){
            s = (a + b) / 2;
            mflag = true;
        }
        else{
            mflag = false;
        }
        double fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0) { b = s; fb = fs; }
        else { a = s; fa = fs; }
        if(fabs(fa) < fabs(fb)){
            std::swap(fa, fb);
            std::swap(a, b);
        }
    }
    throw theta::Exception("brent: Too many iterations");
}

#endif

