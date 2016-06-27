#include "plugins/newton-utils.hpp"
#include "interface/phys.hpp"
#include "interface/distribution.hpp"
#include "interface/redirect_stdio.hpp"
#include <limits>
#include <iostream>

using namespace std;
using namespace theta;
using namespace newton_internal;

double RangedFunction::eval_with_derivative(const vector<double> & x0, vector<double> & grad) const{
    const size_t n = ndim();
    theta_assert(x0.size()==n);
    //theta_assert(ipar < static_cast<int>(n));
    vector<double> x(x0);
    grad.resize(n);
    double f0 = operator()(x);
    for(size_t i=0; i<n; ++i){
        //if(ipar >= 0 and static_cast<int>(i)!=ipar) continue;
        double x0i = x0[i];
        const double sigma = step[i];
        theta_assert(sigma > 0);
        // we choose h ~ sqrt(epsilon_f / f'') where f'' is the second derivative.
        // Here, f''  ~  1 / step**2
        double h_ = sqrt(epsilon_f) * sigma;
        double x0prime = x0i + h_;
        if(x0prime > range_max[i] || x0prime < range_min[i]){
            x0prime = x0i - h_;
            if(x0prime > range_max[i] || x0prime < range_min[i]){
                throw invalid_argument("parameter range is to small");
            }
        }
        theta_assert(!std::isnan(x0prime));
        x[i] = x0prime;
        double fprime = operator()(x);
        x[i] = x0[i];
        theta_assert(!std::isnan(fprime));
        volatile double h = x0prime - x0i;
        double g = (fprime - f0) / h;
        // estimate the error: roundoff error + truncation error:
        //double g_err = epsilon_f / h + h / pow(sigma, 2);
        grad[i] = g;
    }
    return f0;
}

size_t RangedFunction::trunc_to_range(vector<double> & x) const{
    const size_t n = ndim();
    x.resize(n);
    size_t result = 0;
    for(size_t i=0; i<n; ++i){
        if(x[i] > range_max[i]){
            x[i] = range_max[i];
            ++result;
        }
        else if(x[i] < range_min[i]){
            x[i] = range_min[i];
            ++result;
        }
    }
    return result;
}

RangedFunction::~RangedFunction(){}

RangedThetaFunction::RangedThetaFunction(const Function & f_, const ParValues & pv_fixed_values, const ParValues & pv_step, const Ranges & ranges, bool use_f_derivative_):
                f(f_), use_f_derivative(use_f_derivative_){
    epsilon_f = numeric_limits<double>::epsilon();
    const ParIds & pids = f.get_parameters();
    step.reserve(pids.size());
    nonfixed_pids.reserve(pids.size());
    range_min.reserve(pids.size());
    range_max.reserve(pids.size());
    size_t i=0;
    for(ParIds::const_iterator pit=pids.begin(); pit!=pids.end(); ++pit, ++i){
        if(pv_fixed_values.contains(*pit)){
            fixed_parameters.insert(*pit);
            pvs.set(*pit, pv_fixed_values.get_unchecked(*pit));
        }
        else{
            const pair<double, double> & range = ranges.get(*pit);
            if(pv_step.get(*pit)==0.0){
                theta_assert(range.first == range.second);
                fixed_parameters.insert(*pit);
                pvs.set(*pit, range.first);
            }
            else{
                //nonfixed_indices.push_back(i);
                nonfixed_parameters.insert(*pit);
                nonfixed_pids.push_back(*pit);
                step.push_back(pv_step.get(*pit));
                theta_assert(step.back() > 0.0);
                range_min.push_back(range.first);
                range_max.push_back(range.second);
            }
        }
    }
}
    
void RangedThetaFunction::set_epsilon_f(double newval){
    theta_assert(newval > 0);
    epsilon_f = newval;
}

const ParIds & RangedThetaFunction::get_fixed_parameters() const{
    return fixed_parameters;
}

const ParIds & RangedThetaFunction::get_nonfixed_parameters() const{
    return nonfixed_parameters;
}
    
size_t RangedThetaFunction::ndim() const{
    return nonfixed_pids.size();
}

// vals are the ndim() non-fixed parameters in conversion order defined by pids_ as passed to the constructor.
double RangedThetaFunction::operator()(const vector<double> & vals) const {
    //theta_assert(vals.size() == nonfixed_indices.size());
    for(size_t i=0; i<nonfixed_pids.size(); ++i){
        pvs.set(nonfixed_pids[i], vals[i]);
        //values[nonfixed_indices[i]] = vals[i];
    }
    return f(pvs);
}

double RangedThetaFunction::eval_with_derivative(const std::vector<double> & x0, std::vector<double> & grad) const{
    if(use_f_derivative){
        for(size_t i=0; i<nonfixed_pids.size(); ++i){
            pvs.set(nonfixed_pids[i], x0[i]);
        }
        double result = f.eval_with_derivative(pvs, pv_der);
        // convert back:
        for(size_t i=0; i<nonfixed_pids.size(); ++i){
            grad[i] = pv_der.get(nonfixed_pids[i]);
        }
        return result;
    }
    else{
        return RangedFunction::eval_with_derivative(x0, grad);
    }
}

namespace newton_internal{
    
// update the triplet tr using the new point d with f(d) = fd.
// The new point must be within the interval: d > tr.a and d < tr.b
void add_point(min_triplet & tr, double d, double fd){
    theta_assert(d > tr.a and d < tr.b);
    if(d < tr.c){
        // a < d < c < b
        if(fd <= tr.fc){
            // new triplet is a < d < c:
            tr.b = tr.c;
            tr.fb = tr.fc;
            tr.c = d;
            tr.fc = fd;
        }
        else{
            // new triplet is d < c < b:
            tr.a = d;
            tr.fa = fd;
        }
    }
    else{
        // a < c < d < b
        if(fd <= tr.fc){
            // new triplet: c < d < b:
            tr.a = tr.c;
            tr.fa = tr.fc;
            tr.c = d;
            tr.fc = fd;
        }
        else{
            // new triplet: a < c < d
            tr.b = d;
            tr.fb = fd;
        }
    }
}
    
    
}


double newton_internal::bisect_larger(const min_triplet & tr){
    if(tr.b - tr.c > tr.c - tr.a){
        return tr.c + 0.5 * (tr.b - tr.c);
    }
    else{
        return tr.a + 0.5 * (tr.c - tr.a);
    }
}


double newton_internal::quadratic_interpolation(const min_triplet & tr){
    // shift the problem such that c = 0, fc = 0
    double a = tr.a - tr.c;
    double b = tr.b - tr.c;
    double fa = tr.fa - tr.fc;
    double fb = tr.fb - tr.fc;
    // minimum of the parabola in shifted coordinates:
    double x = (b*b * fa - a*a*fb) / (b * fa - a*fb) / 2;
    return x + tr.c;// note: returning inf or nan is ok here; this is handled by find_argmin
}


void newton_internal::find_argmin(const f1d & f, min_triplet & tr, const fproposal & fp,  const fstop & stop, unsigned int maxit){
    unsigned int it = 0;
    const double max_factor = 0.8;
    for(; it < maxit; ++it){
        if(stop(tr)){
            return;
        }
        theta_assert(tr.a <= tr.c and tr.c <= tr.b);
        theta_assert(tr.fa >= tr.fc and tr.fb >= tr.fc);
        // generate the proposal with fp:
        double c_trial = fp(tr);
        // if the proposal is a finite number within a,b then accept it:
        if(isfinite(c_trial) and c_trial > tr.a and c_trial < tr.b){
        }
        else{
            // Make a better suggestion by bisection the larger sub-interval:
            c_trial = bisect_larger(tr);
        }
        const double fct = f(c_trial);
        const double isize_old = tr.b - tr.a;
        add_point(tr, c_trial, fct);
        const double isize_new = tr.b - tr.a;
        if(isize_new > max_factor * isize_old){
            double c_trial2 = bisect_larger(tr);
            double fct2 = f(c_trial2);
            add_point(tr, c_trial2, fct2);
            // There is one problematic situation left: if both c and c_trial are close to an interval end,
            // say a < c < c_trial with only very small differences. Then c could be the new a in tr_tmp. Bisecting
            // the larger sub-interval in this if condition bisects c_trial--b and the final interval of this iteration
            // could pick c_trial as the final new a which means that the interval has hardly shrinked during this iteration.
            // However, in this case, c split a,b into two equal sub-intervals; this guarantees that the next
            // iteration is free from this pathology.
        }
    }
    throw Exception("find_argmin2: too many iterations to find minimum");
}


double newton_internal::f_accuracy(const RangedFunction & f, const vector<double> & x0, size_t i, double f_scale){
    theta_assert(f_scale > 0.0);
    theta_assert(i < x0.size());
    theta_assert(f.ndim() == x0.size());
    vector<double> x(x0);
    size_t f_eval = 0;
    theta_assert(f.trunc_to_range(x) == 0);
    const double f0 = fabs(f(x)) + f_scale;
    ++f_eval;
    double h = numeric_limits<double>::epsilon() * fabs(f0);
    x[i] = x0[i] + h;
    if(f.trunc_to_range(x) > 0){
        h = -h;
        x[i] = x0[i] + h;
        if(f.trunc_to_range(x) > 0) throw invalid_argument("f_accuracy: parameter range to small");
    }
    double f1 = fabs(f(x)) + f_scale;
    ++f_eval;
    // two cases: either f0==f1, then we have to increase h, or f0!=f1, then we have to decrease h:
    bool increase_h = f0 == f1;    
    double f_eps = increase_h ? numeric_limits<double>::infinity() : fabs(f1 - f0);
    for(int j=0; j<1025; ++j){ // 2**+-1024 = double_max / double_min
        if(increase_h) h *= 2;
        else h /= 2;
        x[i] = x0[i] + h;
        if(f.trunc_to_range(x) > 0){
            throw invalid_argument("f_accuracy: parameter range to small");
        }
        double f1 = fabs(f(x)) + f_scale;
        ++f_eval;
        if(f1 != f0){
            f_eps = min(f_eps, fabs(f1 - f0));
            // in increasing mode, this is the first non-zero difference, and we are done:
            if(increase_h){
                return f_eps;
            }
        }
        else{
            // if in decresing mode, this is the first zero difference, and we are done:
            if(!increase_h){
                return f_eps;
            }
        }
    }
    throw logic_error("too many iterations to find f_eps (f constant?)");
}


