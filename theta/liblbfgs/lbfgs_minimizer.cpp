#include "lbfgs_minimizer.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/variables.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include <iostream>

#include <boost/ptr_container/ptr_vector.hpp>

using boost::ptr_vector;

using namespace theta;
using namespace std;

namespace{
    
    class par_trafo{
    public:
        virtual double operator()(double x, double & der) const = 0;
        virtual ~par_trafo(){}
    };
    
    
    class restricted_lower: public par_trafo{
    private:
        double xmin;
    public:
        restricted_lower(double xmin_): xmin(xmin_){}
        
        using par_trafo::operator();
        
        virtual double operator()(double x, double & der) const {
            double xdiff = x - xmin - M_PI_2;
            if(xdiff > 0){
                der = 1.0;
                return x;
            }
            der = 1 / (1 + xdiff * xdiff);
            return xmin + M_PI_2 + atan(xdiff);
        }
    };
    
    class restricted_upper: public par_trafo{
    private:
        double xmax;
    public:
        restricted_upper(double xmax_): xmax(xmax_){}
        
        using par_trafo::operator();
        
        virtual double operator()(double x, double & der) const {
            double xdiff = xmax - x - M_PI_2;
            if(xdiff > 0){
                der = 1.0;
                return x;
            }
            der = 1 / (1 + xdiff * xdiff);
            return xmax - M_PI_2 - atan(xdiff);
        }
    };
    
    // for the case where xmin, xmax is larger than M_PI * 1.1
    class restricted_both_large: public par_trafo{
    private:
        double xmin, xmax;
    public:
        restricted_both_large(double xmin_, double xmax_): xmin(xmin_), xmax(xmax_){
            assert(xmax - xmin > M_PI * 1.09);
        }
        
        using par_trafo::operator();
        
        virtual double operator()(double x, double & der) const {
            double xdiff_lower = x - xmin - M_PI_2;
            double xdiff_upper = xmax - x - M_PI_2;
            if(xdiff_lower > 0){
                if(xdiff_upper > 0){
                    return x;
                    der = 1.0;
                }
                else {
                    der = 1 / (1 + xdiff_upper * xdiff_upper);
                    return xmax - M_PI_2 - atan(xdiff_upper);
                }
            }
            der = 1 / (1 + xdiff_lower * xdiff_lower);
            return xmin + M_PI_2 + atan(xdiff_lower);
        }
    };
    
    
    
    class restricted_both_small: public par_trafo{
    private:
        double xmin, xmax, xwidth;
    public:
        restricted_both_small(double xmin_, double xmax_): xmin(xmin_), xmax(xmax_), xwidth(xmax - xmin){
        }
        
        using par_trafo::operator();
        
        virtual double operator()(double x, double & der) const {
            double x_xwidth = x / xwidth;
            der = M_1_PI / (1 + x_xwidth * x_xwidth);
            return xmin + xwidth / M_PI * atan(x_xwidth);
        }
    };
    
    // identical function for unrestricted parameters:
    class id: public par_trafo{
    public:
        virtual double operator()(double x, double & der) const {
            der = 1.0;
            return x;
        }
    };
    
    
    struct lbfgs_info{
        
        // the non-fixed parameters ...
        vector<ParId> pars;
        ptr_vector<par_trafo> trafos;
        vector<pair<double, double> > ranges;

        // the fixed parameters are always already set in values:
        ParValues values, derivatives;
        const Function & func;
        
        // pars, par_ranges refer to the non-fixed parameters.
        lbfgs_info(const Function & f, const vector<ParId> & pars_, const vector<pair<double, double> > & ranges_, ptr_vector<par_trafo> & trafos_, const ParValues & fixed_values): pars(pars_), ranges(ranges_),
           values(fixed_values), func(f){
            theta_assert(pars.size() == trafos_.size());
            trafos.transfer(trafos.end(), trafos_.begin(), trafos_.end(), trafos_);
        }
        
        // transform the "outer parameters" x (unrestricted) to the "inner parameters" (restricted to par_ranges), stored in values.
        // store the derivatives of this transformation in derivatives.
        void transform(const double * x){
            size_t i=0;
            for(vector<ParId>::const_iterator it=pars.begin(); it!=pars.end(); ++it, ++i){
                double der = NAN;
                if(isnan(x[i])) throw MinimizationException("x was NAN during function evaluation");
                double x_inner = trafos[i](x[i], der);
                values.set(*it, x_inner);
                derivatives.set(*it, der);
            }
        }
    };
    
    
    double evaluate(void * instance, const double *x, double *g, const int n, const double step){
        lbfgs_info * info = static_cast<lbfgs_info*>(instance);
        info->transform(x);
        double result = info->func(info->values);

        size_t i=0;
        const double eps = numeric_limits<double>::epsilon();
        const double se = sqrt(eps);
        for(vector<ParId>::const_iterator it=info->pars.begin(); it!=info->pars.end(); ++it, ++i){
            const ParId & pid = *it;
            double x0 = info->values.get_unchecked(pid);
            double sign = x0 > 0.0 ? 1.0 : -1.0;
            double x0prime = x0 * (1 + se);
            if(fabs(x0prime) < eps) x0prime = sign * eps;
            if(x0prime > info->ranges[i].second || x0prime < info->ranges[i].first){
                x0prime = x0 * (1 - se);
                if(fabs(x0prime) < eps) x0prime = -sign * eps;
                if(x0prime > info->ranges[i].second || x0prime < info->ranges[i].first){
                    throw invalid_argument("lbfgs_minimizer: relative parameter range is to small");
                }
            }
            info->values.set(pid, x0prime);
            double fprime = info->func(info->values);
            volatile double h = x0prime - x0;
            g[i] = (fprime - result) / h;
            g[i] *= info->derivatives.get(*it);
            info->values.set(pid, x0);
        }
        theta_assert(i==(size_t)n);
        //TODO: what is step?
        return result;
    }

}




lbfgs_minimizer::lbfgs_minimizer(const theta::Configuration & cfg){
    lbfgs_parameter_init(&params);
    params.epsilon = 0.0;
    params.delta = 1e-2;
    params.past = 3;
    params.ftol = 1e-4;
    //params.m = 10;
    if(cfg.setting.exists("delta")){
        params.delta = cfg.setting["delta"];
    }
    params.max_linesearch = 200; // default is 20
    params.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    //params.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
    params.max_iterations = 10000;
}

theta::MinimizationResult lbfgs_minimizer::minimize(const theta::Function & f, const theta::ParValues & start,
                                   const theta::ParValues & step, const Ranges & ranges){
    ptr_vector<par_trafo> trafos;
    vector<ParId> non_fixed_pars;
    ParValues fixed_values;
    const ParIds & parameters = f.get_parameters();
    vector<pair<double, double> > ranges_nonfixed;
    for(ParIds::const_iterator it = parameters.begin(); it!=parameters.end(); ++it){
        const ParId & pid = *it;
        const std::pair<double, double> & range = ranges.get(pid);
        if(range.first == range.second){
            double val = start.get(pid);
            if(val != range.first){
                throw invalid_argument("lbfgs_minimizer: start value outside of allowed range given");
            }
            fixed_values.set(pid, val);
        }
        else{
            if(step.get(pid) <= 0.0){
                throw invalid_argument("lbfgs_minimizer: non-empty range given, but step width is <=0.0");
            }
            non_fixed_pars.push_back(pid);
            ranges_nonfixed.push_back(range);
            if(isinf(range.first)){
                if(isinf(range.second)){
                    trafos.push_back(new id());
                }
                else{
                    trafos.push_back(new restricted_upper(range.second));
                }
            }
            else{
                if(isinf(range.second)){
                    trafos.push_back(new restricted_lower(range.first));
                }
                else{
                    if(range.second - range.first > 1.1 * M_PI){
                        trafos.push_back(new restricted_both_large(range.first, range.second));
                    }
                    else{
                        trafos.push_back(new restricted_both_small(range.first, range.second));
                    }
                }
            }
        }
    }
    theta_assert(non_fixed_pars.size() == trafos.size());
    theta_assert(non_fixed_pars.size() == ranges_nonfixed.size());
    const size_t n_nonfixed = non_fixed_pars.size();
    
    vector<double> x(n_nonfixed);
    for(size_t i=0; i<n_nonfixed; ++i){
        x[i] = start.get(non_fixed_pars[i]);
    }
    
    theta::MinimizationResult result;
    lbfgs_info info(f, non_fixed_pars, ranges_nonfixed, trafos, fixed_values);
    int status = lbfgs(n_nonfixed, &x[0], &result.fval, evaluate, 0, &info, &params);
    if(status < 0){
        stringstream ss;
        ss << "liblbfgs returned error code " << status;
        throw MinimizationException(ss.str());
    }
    info.transform(&x[0]);
    result.values.set(info.values);
    return result;
}

REGISTER_PLUGIN(lbfgs_minimizer)
