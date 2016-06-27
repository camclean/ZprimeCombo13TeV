#include "plugins/newton.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/phys.hpp"
#include "interface/plugin.hpp"
#include "interface/variables.hpp"
#include "interface/matrix.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/redirect_stdio.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <limits>

using namespace theta;
using namespace std;
using namespace newton_internal;


namespace newton_internal {

class Linesearch{
public:
    // do a line search for a minimum starting at x0, with an estimate for the minimum at x0 + step.
    // The result is filled in x.
    // x does not have to be on the line defined by x0 and step. (this allows an easier
    // handling of constraints). The result point x must be in the range of f.
    //
    // returns the function value at x.
    virtual double do_linesearch(const RangedFunction & f, const std::vector<double> & x0, const std::vector<double> & step, std::vector<double> & x, bool precise = false) const = 0;
    virtual ~Linesearch(){}
};


}


namespace{


// Function info for the newton Minimizer. So far, we only add epsilon_f.
class NewtonFunctionInfo: public FunctionInfo{
private:
    double epsilon_f;
    Matrix inverse_hessian;
    
public:
    NewtonFunctionInfo(const ParValues & start, const ParValues & step, const Ranges & ranges, const ParValues & fixed_parameters, double epsilon_f_):
        FunctionInfo(start, step, ranges, fixed_parameters), epsilon_f(epsilon_f_){
        theta_assert(epsilon_f > 0.0);
    }
    
    NewtonFunctionInfo(const FunctionInfo & info, double epsilon_f_): FunctionInfo(info), epsilon_f(epsilon_f_){
        theta_assert(epsilon_f > 0.0);
    }
    
    double get_epsilon_f() const{
        return epsilon_f;
    }
    
    void set_inverse_hessian(const Matrix & m){
        inverse_hessian = m;
    }
    
    Matrix get_inverse_hessian() const {
        return inverse_hessian;
    }
};


class IntervalLinesearch: public Linesearch {
private:
    double eps;
    bool debug;
    
    struct eval_f{
        const RangedFunction & f;
        const vector<double> & x0, step;
        vector<double> x;
        
        eval_f(const RangedFunction & f_, const vector<double> & x0_, const vector<double> & step_): f(f_), x0(x0_), step(step_){
            x.resize(x0.size());
        }
        
        double operator()(double c){
            x = x0;
            add_with_coeff(x, step, c);
            f.trunc_to_range(x);
            return f(x);
        }
    };
    
    struct ab_small{
        const RangedFunction & f;
        const vector<double> & x0, step;
        double eps;
        vector<double> x, y;
        bool precise;
        
        ab_small(const RangedFunction & f_, const vector<double> & x0_, const vector<double> & step_, double eps_, bool precise_): f(f_), x0(x0_), step(step_), eps(eps_), precise(precise_){
            x.resize(x0.size());
            y.resize(x0.size());
        }
        
        bool operator()(const min_triplet & tr){
            // small enough means that the interval a,b is known to a precision of 10%.
            // This only applies if the interval does not switch signs (or is at 0.0):
            if(!precise){
                if(tr.a * tr.b > 0.0){
                    double denom = min(fabs(tr.a), fabs(tr.b));
                    if((tr.b - tr.a) / denom < 0.1) return true;
                }
            }
            // Otherwise, test the distance (pnorm) in truncated parameter space:
            x = x0;
            add_with_coeff(x, step, tr.a);
            f.trunc_to_range(x);
            y = x0;
            add_with_coeff(y, step, tr.b);
            f.trunc_to_range(y);
            add_with_coeff(x, y, -1.0);
            return pnorm(x, step) < eps;
        }
    };    
    
public:
    explicit IntervalLinesearch(double eps_, bool debug_ = false): eps(eps_), debug(debug_){}
    
    double do_linesearch(const RangedFunction & f, const vector<double> & x0, const vector<double> & step, vector<double> & x_new, bool precise) const{
        const double f0 = f(x0);
        double fbest = f0;
        double cbest = 0.0; // best coefficient c found so far for   x = x0 + c * step
        // 1. find start of coefficients (a,b) in which the minimum is contained.
        // We assume that step is already a good start, so test nearby points:
        // In general, delta0 could be tuned; a quick test suggests that either the choice
        // does not matter much, or this choice is already pretty good.
        const double delta0 = 0.2;
        const size_t n = f.ndim();
        double a = 1 - delta0;
        
        // 1.a. make sure the lower point a is within the function range:
        vector<double> x(x0);
        add_with_coeff(x, step, a);
        // search for a value of a which is not outside the function range (that should not happen very often):
        bool found_a = false;
        for(int i=0; i < int(1/delta0 + 2); ++i){
            if(f.trunc_to_range(x) < n){
                found_a = true;
                break;
            }
            a -= delta0;
            x = x0;
            add_with_coeff(x, step, a);
        }
        theta_assert(found_a);
        
        // 1.b. find a, cbest such that a < cbest   and f(cbest) <= f(a).
        // the 1d function, depending on the coefficient of the step vector:
        eval_f f1d(f, x0, step);
        double fa = f1d(a);
        double delta = delta0;
        // find  a  such that fa > fbest and a < cbest
        while(fa < fbest || a >= cbest){
            if(fa < fbest){
                fbest = fa;
                cbest = a;
            }
            a -= delta;
            fa = f1d(a);
            delta *= 2;
        }
        theta_assert(fa >= fbest);
        theta_assert(a < cbest);
        // note that in rare cases, a can become negative, but that's Ok.
        
        // 1.c. find b > cbest with f(xb) > f(xbest)
        delta = max(cbest - a, delta0);
        double b = cbest + delta;
        double fb = f1d(b);
        while(fb < fbest){
            a = cbest;
            fa = fbest;
            fbest = fb;
            cbest = b;
            b += delta;
            fb = f1d(b);
            delta *= 2;
        }
        // 2. using a < c < b with fa >= fc <= fb, search for the minimum of f.
        ab_small small_enough(f, x0, step, eps, precise);
        min_triplet tr;
        tr.a = a; tr.b=b; tr.c = cbest;
        tr.fa = fa; tr.fb = fb; tr.fc = fbest;
        find_argmin(boost::ref(f1d), tr, quadratic_interpolation, boost::ref(small_enough), 10000);
        x_new = x0;
        add_with_coeff(x_new, step, tr.c);
        f.trunc_to_range(x_new);
        return fbest;
    }
};



} // anon. namespace

   
std::ostream & operator<<(std::ostream & out, const std::vector<double> & x){
    const size_t n = x.size();
    for(size_t i=0; i<n; ++i){
        out << x[i] << " ";
    }
    return out;
}


std::ostream & operator<<(std::ostream & out, const Matrix & m){
    const size_t nr = m.get_n_rows();
    const size_t nc = m.get_n_cols();
    for(size_t i=0; i<nr; ++i){
        for(size_t j=0; j<nc; ++j){
            out << m(i,j) << " ";
        }
        out << endl;
    }
    return out;
}

// PROBLEM: operator<< is defined globally


newton_minimizer::newton_minimizer(const newton_minimizer:: options & opts_): opts(opts_), use_nll_der(false){
    ls.reset(new IntervalLinesearch(opts.par_eps, opts.debug & 1));
}

newton_minimizer::newton_minimizer(const Configuration & cfg): use_nll_der(false){
    Setting s = cfg.setting;
    if(s.exists("maxit")){
        opts.maxit = s["maxit"];
    }
    if(s.exists("step_cov")){
        opts.step_cov = s["step_cov"];
    }
    if(s.exists("improve_cov")){
        opts.improve_cov = s["improve_cov"];
    }
    if(s.exists("force_cov_positive")){
        opts.force_cov_positive = s["force_cov_positive"];
    }
    if(s.exists("par_eps")){
        opts.par_eps = s["par_eps"];
    }
    if(s.exists("debug")){
        opts.debug = s["debug"];
    }
    if(s.exists("use_nll_der")){
        use_nll_der = s["use_nll_der"];
    }
    if(s.exists("second_pass")){
        opts.second_pass = s["second_pass"];
    }
    ls.reset(new IntervalLinesearch(opts.par_eps, opts.debug & 1));
}


MinimizationResult newton_minimizer::minimize2(const theta::Function & f_, const theta::FunctionInfo & info_, const theta::ParValues & fixed_parameters){
    const NewtonFunctionInfo & info = dynamic_cast<const NewtonFunctionInfo &>(info_);
    RangedThetaFunction f(f_, fixed_parameters, info.get_step(), info.get_ranges(), use_nll_der);
    f.set_epsilon_f(info.get_epsilon_f());
    const ParIds & all_pids = f_.get_parameters();
    const ParIds & fixed_pids = f.get_fixed_parameters();
    theta_assert(fixed_pids == info.get_fixed_parameters());
    ParIds nonfixed_pids;
    for(ParIds::const_iterator it=all_pids.begin(); it!=all_pids.end(); ++it){
        if(fixed_pids.contains(*it)) continue;
        nonfixed_pids.insert(*it);
    }
    const size_t n = f.ndim();
    theta_assert(n == nonfixed_pids.size());
    // fill the inverse hessian either from info or using step**2 at the diagonals:
    ParValues pv_step = info.get_step();
    Matrix inverse_hessian = info.get_inverse_hessian();
    bool ih_info = inverse_hessian.get_n_rows() == n;
    if(!ih_info){
        inverse_hessian.reset(n,n);
    }
    vector<double> step0(n);
    {
        size_t i=0;
        for(ParIds::const_iterator it=nonfixed_pids.begin(); i<n; ++i, ++it){
            double st = pv_step.get(*it);
            if(!ih_info){
                inverse_hessian(i,i) = pow(st, 2);
            }
            step0[i] = st;
        }
    }
    vector<double> x(n), grad(n), z(n);
    vector<double> direction(n), next_x(n), next_grad(n), dx(n), step(step0);
    ParValues pv_start = info.get_start();
    pv_start.set(fixed_parameters);
    pv_start.fill(&x[0], nonfixed_pids);
    double fx = f.eval_with_derivative(x, grad);
    int it = 0;
    for(; it < opts.maxit; ++it){
        if(opts.debug){
            out << "Iteration " << it << " starting.";
            if(opts.debug >=2){
                out << endl << "x = " << x << " --> f(x) = " << fx << endl;
                out << "g = " << grad;
            }
            if(opts.debug >= 4){
                out << endl << "h = " << endl;
                char s[200];
                for(size_t k=0; k<n; ++k){
                    for(size_t l=0; l<n; ++l){
                        snprintf(s, 200, "%12.4g", inverse_hessian(k, l));
                        out << s;
                    }
                    out << endl;
                }
            }
            out << endl;
        }
        
        // the estimated minimum is   -inverse_hessian * grad
        mul(direction, inverse_hessian, grad);
        mul(direction, -1.0);
        if(opts.debug >= 2){
            out << "Iteration " << it << ": search direction = " << direction << endl;
        }
        if(pnorm(direction, step) < opts.par_eps){
            if(opts.debug){
                out << "Iteration " << it << ": estimated minimum close enough." << endl;
            }
            break;
        }
        if(opts.debug) out << "Iteration " << it << ": doing linesearch now" << endl;
        double next_fx = ls->do_linesearch(f, x, direction, next_x);
        if(opts.debug >= 2){
            out << "Iteration " << it << ": linesearch proposes:" << endl << "x =" << next_x << endl;
        }
        // calculate dx = next_x - x
        dx = next_x;
        add_with_coeff(dx, x, -1.0);
        if(opts.debug){
            out << "Iteration " << it << ": step size of linesearch was pnorm = " << pnorm(dx, step) << "; norm = " << norm(dx) << endl;
        }
        // we are done if the actual step size was small. However, the line search is only accurate to 10% of the search direction, so
        // make a more precise linesearch again:
        if(pnorm(dx, step) < opts.par_eps){
            if(opts.debug){
                 out << "Iteration " << it << ": actual dx by linesearch was smaller than par_eps, doing a more precise linesearch." << endl;
            }
            ls->do_linesearch(f, x, direction, next_x, true);
            dx = next_x;
            add_with_coeff(dx, x, -1.0);
            if(pnorm(dx, step) < opts.par_eps){
                if(opts.debug){
                    out << "Iteration " << it << ": precise linesearch also below par_eps, stopping." << endl;
                }
                break;
            }
            else{
                if(opts.debug){
                    out << "Iteration " << it << ": precise linesearch was above par_eps, continuing." << endl;
                }
            }
        }
        // do the update to inverse_hessian:
        next_fx = f.eval_with_derivative(next_x, next_grad);
        // calculate z = (Delta x - inverse_hessian * (next_grad - grad))
        for(size_t i=0; i<n; ++i){
            z[i] = dx[i];
            for(size_t j=0; j<n; ++j){
                z[i] -=  inverse_hessian(i, j) * (next_grad[j] - grad[j]);
            }
        }
        // calculate the denominator, z * (next_grad - grad):
        double denom = 0.0;
        double norm_graddiff = 0.0;
        for(size_t i=0; i<n; ++i){
            double graddiff = next_grad[i] - grad[i];
            norm_graddiff += graddiff * graddiff;
            denom += z[i] * graddiff;
        }
        double norm_z = norm(z);
        norm_graddiff = sqrt(norm_graddiff);
        if(fabs(denom) > 1e-7 * norm_graddiff * norm_z){
            // calculate the update to inverse_hessian: add z * z^T / denom:
            for(size_t i=0; i<n; ++i){
                for(size_t j=0; j< n; ++j){
                    inverse_hessian(i,j) += z[i] * z[j] / denom;
                }
            }
        }
        // update step:
        for(size_t i=0; i<n; ++i){
            if(inverse_hessian(i,i) > pow(step0[i], 2)){
                step[i] = sqrt(inverse_hessian(i,i));
            }
        }
        // prepare for next iteration:
        swap(x, next_x);
        swap(grad, next_grad);
        fx = next_fx;
    }
    if(it==opts.maxit){
        if(opts.debug) out << "maximum number of iterations reached; throwing exception" << endl;
        throw MinimizationException("maximum number of iterations reached");
    }
    // cov is the best estimate for the covariance matrix: either inverse_hessian directly, or
    // the explicit inverse of the numerical hessian in case of improve_cov:
    Matrix cov = inverse_hessian;
    if(opts.improve_cov){
        if(opts.debug) out << "improve_cov is true, calculating Hesse matrix now" << endl;
        Matrix hessian(n,n);
        // calculate better estimate for covariance here by going to opts.step_cov times the step size in the
        // direction of the parameter i and evaluate the gradient there; repeat for all i.
        theta_assert(f.trunc_to_range(x)==0);
        for(size_t i=0; i<n; ++i){
            next_x = x;
            double stepi = opts.step_cov * step[i];
            next_x[i] = x[i] + stepi;
            if(f.trunc_to_range(next_x) > 0){
                stepi *= -1.0;
                next_x[i] = x[i] + stepi;
                if(f.trunc_to_range(next_x) > 0){
                    throw MinimizationException("for calculating covariance: range too small (< step_cov * sigma)");
                }
            }
            if(opts.debug) out << " step" << i << " = " << stepi << endl;
            f.eval_with_derivative(next_x, next_grad);
            for(size_t j=0; j<n; ++j){
                hessian(i,j) = (next_grad[j] - grad[j]) / stepi;
            }
        }
        // symmetrize the off-diagonals: any difference here comes from rounding/truncation errors:
        for(size_t i=0; i<n; ++i){
            for(size_t j=i+1; j<n; ++j){
                hessian(i,j) = hessian(j,i) = 0.5 * (hessian(i,j) + hessian(j,i));
            }
        }
        // make covariance matrix positive definite by adding alpha * identify matrix with small alpha:
        double min_ = numeric_limits<double>::infinity(), max_ = -numeric_limits<double>::infinity();
        for(size_t i=0; i<n; ++i){
            min_ = min(min_, hessian(i,i));
            max_ = max(max_, hessian(i,i));
        }
        if(opts.force_cov_positive){
            // min is the smallest eigenvalue of inverse_hessian. If it is negative (or positive but very small), add a small
            // value to make it > sqrt(eps) * max:
            double desired_min = sqrt(numeric_limits<double>::epsilon()) * fabs(max_);
            if(min_ < desired_min){
                if(opts.debug){
                    out << "hesse not positive definite; hesse = " << hessian << endl;
                    out << "adding " << (desired_min - min_) << " to diagonal of Hesse to make it positive definite" << endl;
                }
                for(size_t i=0; i<n; ++i){
                    hessian(i,i) += desired_min - min_;
                }    
            }
        }
        try{
            hessian.invert_cholesky();
            cov = hessian;
        }
        catch(range_error & re){
            stringstream s;
            s << "error inverting matrix to calculate covariance matrix: " << re.what();
            throw MinimizationException(s.str());
        }
        
    }
    MinimizationResult res;
    res.fval = fx;
    res.values.set(pv_start); // to get the fixed parameters right
    res.values.set(ParValues(&x[0], nonfixed_pids));
    // write the covariance we have so far; make sure to use 0 for fixed parameters:
    res.covariance = Matrix(all_pids.size(), all_pids.size());
    size_t i_nf = 0, j_nf = 0; // count non-fixed parameters
    size_t i=0; // index to all parameters
    for(ParIds::const_iterator pit=all_pids.begin(); pit!=all_pids.end(); ++pit, ++i){
        size_t j=0;
        j_nf = 0;
        for(ParIds::const_iterator pit2=all_pids.begin(); pit2!=all_pids.end(); ++pit2, ++j){
            if(nonfixed_pids.contains(*pit) && nonfixed_pids.contains(*pit2)){
                res.covariance(i,j) = cov(i_nf, j_nf);
            }
            if(i==j){
                double error = sqrt(fabs(cov(i_nf, i_nf)));
                if(cov(i_nf, i_nf) < 0.0){
                    error *= -1.0; // let the user see if there's something weird
                }
                res.errors_plus.set(*pit, error);
                res.errors_minus.set(*pit, error);
            }
            if(nonfixed_pids.contains(*pit2)) ++j_nf;
        }
        theta_assert(j_nf == nonfixed_pids.size());
        theta_assert(j == all_pids.size());
        if(nonfixed_pids.contains(*pit)) ++i_nf;
    }
    theta_assert(i == all_pids.size());
    theta_assert(i_nf == nonfixed_pids.size());
    return res;
}


MinimizationResult newton_minimizer::minimize(const theta::Function & f_, const theta::ParValues & start,
                                              const theta::ParValues & step, const Ranges & ranges){
    ParValues fixed_parameters;
    RangedThetaFunction f(f_, fixed_parameters, step, ranges, use_nll_der);
    // get epsilon_f:
    vector<double> x0(f.ndim());
    start.fill(&x0[0], f.get_nonfixed_parameters());
    double epsilon_f = f_accuracy(f, x0, 0);
    if(opts.debug) out << "epsilon_f = " << epsilon_f << endl;
    NewtonFunctionInfo info(start, step, ranges, fixed_parameters, epsilon_f);
    try{
        MinimizationResult minres = minimize2(f_, info, fixed_parameters);
        if(opts.second_pass){
            NewtonFunctionInfo info2(minres.values, step, ranges, fixed_parameters, epsilon_f);
            minres = minimize2(f_, info2, fixed_parameters);
        }
        return minres;
    }
    catch(const logic_error & ex){
        // log_errors can occur in case the function to minimize returns inf or nan. Treat this as a MinimizationException:
        throw MinimizationException(string("logic errors during minimization: ") + ex.what());
    }
}

boost::shared_ptr<theta::FunctionInfo> newton_minimizer::create_nll_function_info(const theta::Model & m,
                                                                                  const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
                                                                                  const theta::ParValues & fixed_parameters){
    boost::shared_ptr<FunctionInfo> info = Minimizer::create_nll_function_info(m, override_parameter_distribution, fixed_parameters);
    // get epsilon_f for the asimov likelihood:
    asimov_data_nll nll(m, override_parameter_distribution, fixed_parameters);
    RangedThetaFunction f(nll, fixed_parameters, info->get_step(), info->get_ranges(), use_nll_der);
    const size_t n = f.ndim();
    vector<double> x0(n);
    ParValues start = info->get_start();
    start.fill(&x0[0], f.get_nonfixed_parameters());
    double epsilon_f = f_accuracy(f, x0, 0);
    boost::shared_ptr<NewtonFunctionInfo> result(new NewtonFunctionInfo(*info, epsilon_f));
    // build approximate inverse hesse from step values:
    Matrix inverse_hessian(n,n);
    size_t i=0;
    ParIds nonfixed_pids = f.get_nonfixed_parameters();
    ParValues step = info->get_step();
    for(ParIds::const_iterator pit=nonfixed_pids.begin(); pit!=nonfixed_pids.end(); ++pit, ++i){
        inverse_hessian(i,i) = pow(step.get(*pit), 2);
    }
    result->set_inverse_hessian(inverse_hessian);
    return result;
}


REGISTER_PLUGIN(newton_minimizer)

