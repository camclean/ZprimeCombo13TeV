#include <boost/test/unit_test.hpp>

#include "interface/matrix.hpp"
#include "interface/phys.hpp"
#include "interface/variables.hpp"
#include "interface/random.hpp"

#include "plugins/newton.hpp"

#include <vector>
#include <iostream>


using namespace theta;
using namespace std;
using namespace newton_internal;

/*
namespace{
    void dump(ostream & out, const ParValues & vals, const VarIdManager & vm){
        ParIds pids = vm.get_all_parameters();
        for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
            try{
                double value = vals.get(*it);
                std::string name = vm.get_name(*it);
                out << " " << name << " = " << value;
            }
            catch(invalid_argument &){
                //ignore
            }
        }
    }
}*/

namespace{

class quadratic_function: public theta::Function{
    Matrix hesse;
    ParValues x0;
    mutable size_t n_eval;
public:
    quadratic_function(const ParIds & pids, const ParValues & x0_, const Matrix & hesse_): hesse(hesse_), x0(x0_), n_eval(0){
        theta_assert(hesse.get_n_cols() == hesse.get_n_rows());
        theta_assert(hesse.get_n_cols() == pids.size());
        theta_assert(x0.contains_all(pids));
        par_ids = pids;
    }

    size_t get_n_eval() const{
        return n_eval;
    }
    void reset_n_eval(){
        n_eval = 0;
    }

    using Function::operator();

    virtual double operator()(const ParValues & x) const{
        ++n_eval;
        const size_t n = par_ids.size();
        std::vector<double> xx0(n); // x-x0
        size_t i=0;
        for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
                xx0[i] = x.get(*it) - x0.get(*it);
                theta_assert(std::isfinite(xx0[i]));
        }
        double result =0.;
        for(i=0; i<n; ++i){
                for(size_t j=0; j<n; ++j){
                        result += 0.5 * hesse(i,j) * xx0[j] * xx0[i];
                }
        }
        theta_assert(!std::isnan(result));
        return result;
    }

};


// x = x**2
class QuadraticRangedFunction: public RangedFunction{
size_t n;
public:
    explicit QuadraticRangedFunction(size_t n_, double min_, double max_): n(n_){
        epsilon_f = numeric_limits<double>::epsilon();
        step.resize(n);
        fill(step.begin(), step.end(), 1.0);
        range_min.resize(n);
        fill(range_min.begin(), range_min.end(), min_);
        range_max.resize(n);
        fill(range_max.begin(), range_max.end(), max_);
    }
    
    size_t ndim() const{
        return n;
    }
    
    double operator()(const vector<double> & x) const {
        theta_assert(x.size() == n);
        double result = 0.0;
        for(size_t i=0; i<n; ++i){
            theta_assert(x[i] >= range_min[i] and x[i] <= range_max[i]);
            result += x[i] * x[i];
        }
        return result;
    }
};


const double inf = std::numeric_limits<double>::infinity();

vector<ParId> create_pids(VarIdManager & vm, const std::string & basename, size_t n){
    vector<ParId> result;
    result.reserve(n);
    for(size_t i=0; i<n; ++i){
        stringstream ss;
        ss << basename << i;
        result.push_back(vm.create_par_id(ss.str()));
    }
    return result;
}

}


BOOST_AUTO_TEST_SUITE(newton_utils)

BOOST_AUTO_TEST_CASE(test_add_with_coeff){
    const size_t n = 9;
    vector<double> x(n), y(n), expected(n);
    double c = M_PI;
    for(size_t i=0; i<n; ++i){
        x[i] = i + 0.1;
        y[i] = -i + 10.1;
        expected[i] = x[i] + c*y[i];
    }
    add_with_coeff(x, y, c);
    for(size_t i=0; i<n; ++i){
        BOOST_CHECK_CLOSE(x[i], expected[i], 1e-10);
    }
    x.resize(n/2);
    BOOST_CHECK_THROW(add_with_coeff(x, y, c), logic_error);
}

BOOST_AUTO_TEST_CASE(test_mul){
    const size_t n = 9;
    Matrix m(n, n);
    vector<double> expected(n);
    vector<double> x(n);
    for(size_t i=0; i<n; ++i){
        x[i] = i*i-M_PI;
    }
    for(size_t i=0; i<n; ++i){
        expected[i] = 0.0;
        for(size_t j=0; j<n; ++j){
            m(i, j) = i+j+0.2;
            expected[i] += m(i,j) * x[j];
        }
    }
    
    vector<double> y;
    BOOST_REQUIRE(y.size()==0);
    mul(y, m, x);
    BOOST_REQUIRE(y.size()==n);
    for(size_t i=0; i<n; ++i){
        BOOST_CHECK_CLOSE(y[i], expected[i], 1e-10);
    }
    x.resize(n-1);
    BOOST_CHECK_THROW(mul(y, m, x), logic_error);
    
    // mul with vector and coefficient:
    double c = M_PI;
    for(size_t i=0; i<n; ++i){
        expected[i] *= c;
    }
    mul(y, c);
    BOOST_REQUIRE(y.size()==n);
    BOOST_REQUIRE(expected.size()==n);
    for(size_t i=0; i<n; ++i){
        BOOST_CHECK_CLOSE(y[i], expected[i], 1e-10);
    }
}

BOOST_AUTO_TEST_CASE(test_norm){
    const size_t n = 9;
    vector<double> v1(n), step(n);
    for(size_t i=0; i<n; ++i){
        v1[i] = i + 0.5 * n;
        step[i] = i + 1.0;
    }
    double expected_norm = 0.0;
    double expected_pnorm = 0.0;
    for(size_t i=0; i<n; ++i){
        expected_norm += v1[i] * v1[i];
        expected_pnorm = max(expected_pnorm, fabs(v1[i] / step[i]));
    }
    expected_norm = sqrt(expected_norm);
    
    BOOST_CHECK_CLOSE(expected_norm, norm(v1), 1e-10);
    BOOST_CHECK_CLOSE(expected_pnorm, pnorm(v1, step), 1e-10);
    step.resize(n-1);
    BOOST_CHECK_THROW(pnorm(v1, step), logic_error);
}


BOOST_AUTO_TEST_CASE(test_rangedf){
   const size_t n = 2;
   QuadraticRangedFunction f(n, -5.0, 6.0);
   vector<double> x(n);
   for(size_t i=0; i<n; ++i){
       x[i] = 0.0;
   }
   BOOST_CHECK(f(x)==0.0);
   x[0] = 1.0;
   BOOST_CHECK(f(x)==1.0);
   
   // check truncation:
   BOOST_CHECK(f.trunc_to_range(x) == 0);
   x[0] = -10.0;
   BOOST_CHECK(f.trunc_to_range(x) == 1);
   BOOST_CHECK(x[0] == -5.0);
   x[0] = 11.0;
   BOOST_CHECK(f.trunc_to_range(x) == 1);
   BOOST_CHECK(x[0] == 6.0);
   
   vector<double> x_before(x);
   BOOST_CHECK(f.trunc_to_range(x) == 0);
   for(size_t i=0; i<n; ++i){
       BOOST_CHECK_EQUAL(x_before[i], x[i]);
   }
}

BOOST_AUTO_TEST_CASE(test_accuracy){
   const size_t n = 2;
   double min_ = -5.0, max_ = 6.0;
   QuadraticRangedFunction f(n, min_, max_);
   vector<double> x(n, 0.0);
   double acc = f_accuracy(f, x, 0);
   const double eps = numeric_limits<double>::epsilon();
   // acc should be around eps. We don't care about a factor 2, so tolerate up to 3 here:
   BOOST_CHECK_GT(acc, eps / 3);
   BOOST_CHECK_LT(acc, eps * 3);
   
   // test that accuracy also works at the borders:
   x[0] = min_;
   double fx = f(x);
   acc = f_accuracy(f, x, 0);
   double expected_acc = (fabs(fx) + 1.0) * eps;
   BOOST_CHECK_GT(acc, expected_acc / 3);
   BOOST_CHECK_LT(acc, expected_acc * 3);
   
   x[0] = max_;
   fx = f(x);
   acc = f_accuracy(f, x, 0);
   expected_acc = (fabs(fx) + 1.0) * eps;
   BOOST_CHECK_GT(acc, expected_acc / 3);
   BOOST_CHECK_LT(acc, expected_acc * 3);
   
   x[0] = max_ + 1.0;
   BOOST_CHECK_THROW(f_accuracy(f, x, 0), logic_error);
}

BOOST_AUTO_TEST_CASE(rangedthetaf1d){
    VarIdManager vm;
    vector<ParId> vpids = create_pids(vm, "p", 1);
    ParId pid = vpids[0];
    ParValues x0;
    x0.set(pid, 0.0);
    
    ParIds pids;
    pids.insert(pid);

    Matrix hesse(1,1);
    hesse(0,0) = 1.0;
    quadratic_function thetaf(pids, x0, hesse);
    
    ParValues pv_x, fixed, step;
    Ranges ranges;
    step.set(pid, 1.0);
    ranges.set(pid, make_pair(-inf, inf));
    RangedThetaFunction f(thetaf, fixed, step, ranges, false);
    
    vector<double> x(1);
    for(int i=0; i<10; ++i){
        x[0] = i-5.0;
        pv_x.set(pid, x[0]);
        double tf = thetaf(pv_x);
        double rf = f(x);
        BOOST_CHECK_EQUAL(tf, rf);
    }
    
    thetaf.reset_n_eval();
    x[0] = 1.0;
    vector<double> g;
    f.eval_with_derivative(x, g);
    BOOST_REQUIRE(g.size() == f.ndim());
    BOOST_CHECK_CLOSE(g[0], 1.0, 1e-6);
    BOOST_CHECK_EQUAL(thetaf.get_n_eval(), 1 + f.ndim());
}

BOOST_AUTO_TEST_CASE(rangedthetaf10d){
    VarIdManager vm;
    const size_t n = 10;
    vector<ParId> vpids = create_pids(vm, "p", n);
    ParValues x0;
    ParIds pids;
    for(size_t i=0; i<n; ++i){
        x0.set(vpids[i], 0.0);
        pids.insert(vpids[i]);
    }

    Matrix hesse(n,n);
    for(size_t i=0; i<n; ++i){
        hesse(i,i) = 1.0;
    }
    
    quadratic_function thetaf(pids, x0, hesse);
    
    ParValues pv_x, step;
    Ranges ranges;
    for(size_t i=0; i<n; ++i){
        step.set(vpids[i], 1.0);
        ranges.set(vpids[i], make_pair(-inf, inf));
    }
    // 1. everything free:
    {
        ParValues fixed;
        RangedThetaFunction f(thetaf, fixed, step, ranges, false);
        ParIds fixed_ids = f.get_fixed_parameters();
        BOOST_CHECK(fixed_ids.size() == 0);
        ParIds nonfixed = f.get_nonfixed_parameters();
        BOOST_CHECK(nonfixed.size() == n);
    }
    // 2. fix vpids[0] via step/range and vpids[n/2] via fixed:
    {
        ParValues fixed;
        step.set(vpids[0], 0.0);
        ranges.set(vpids[0], make_pair(0.0, 0.0));
        fixed.set(vpids[n/2], 0.0);
        RangedThetaFunction f(thetaf, fixed, step, ranges, false);
        ParIds fixed_ids = f.get_fixed_parameters();
        BOOST_REQUIRE(fixed_ids.size() == 2);
        ParIds::const_iterator it = fixed_ids.begin();
        BOOST_CHECK(*it == vpids[0]);
        ++it;
        BOOST_CHECK(*it == vpids[n/2]);
        ParIds nonfixed = f.get_nonfixed_parameters();
        BOOST_CHECK(nonfixed.size() == n-2);
    }
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(newton)

BOOST_AUTO_TEST_CASE(test1d){
    BOOST_CHECKPOINT("test1d entry");
    VarIdManager vm;
    ParId pid = vm.create_par_id("p0");
    ParValues x0;
    x0.set(pid, 0.0);
    ParIds pids;
    pids.insert(pid);

    Matrix hesse(1,1);
    hesse(0,0) = 1.0;

    BOOST_CHECKPOINT("test1d");
    quadratic_function f(pids, x0, hesse);
    newton_minimizer::options opts;
    //opts.debug = true;
    newton_minimizer min(opts);

    BOOST_CHECKPOINT("test1d");
    ParValues start, step;
    start.set(pid, 1.0);
    step.set(pid, 1.0);
    std::map<ParId, pair<double, double> > ranges;
    ranges[pid] = make_pair(-inf, inf);
    MinimizationResult minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(fabs(minres.values.get(pid)), opts.par_eps * step.get(pid));
    //cout << "eval = " << f.get_n_eval() << endl;
    BOOST_CHECK_LT(f.get_n_eval(), 100);
    //cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;

    BOOST_CHECKPOINT("test1d");
    // step size small:
    start.set(pid, 10.0);
    step.set(pid, 0.1);
    f.reset_n_eval();
    BOOST_REQUIRE_EQUAL(f.get_n_eval(), 0);
    minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(fabs(minres.values.get(pid)), 1.01 * opts.par_eps * step.get(pid));
    BOOST_CHECK_LT(f.get_n_eval(), 100);
    //cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;

    BOOST_CHECKPOINT("test1d");
    // step size large:
    start.set(pid, 0.1);
    step.set(pid, 100.0);
    f.reset_n_eval();
    minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(fabs(minres.values.get(pid)), 1.01 * opts.par_eps * step.get(pid));
    BOOST_CHECK_LT(f.get_n_eval(), 100);
    //cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;
}

BOOST_AUTO_TEST_CASE(test10d){
    VarIdManager vm;
    std::vector<ParId> vpids;
    ParIds pids;
    const size_t ndim = 10;
    ParValues x0, start, step, x;
    Matrix hesse(ndim,ndim);
    map<ParId, pair<double, double> > ranges;
    for(size_t i=0; i<ndim; ++i){
        stringstream ss;
        ss << "p" << i;
        ParId pid = vm.create_par_id(ss.str());
        vpids.push_back(pid);
        x0.set(pid, 1.0);
        start.set(pid, 1.0);
        step.set(pid, 1.0);
        x.set(pid, 0.0);
        hesse(i,i) = 1.0;
        pids.insert(pid);
        ranges[pid] = make_pair(-inf, inf);
    }
    quadratic_function f(pids, x0, hesse);

    newton_minimizer::options opts;
    newton_minimizer min(opts);

    /*cout << "f(x0) = " << f(x0) << endl;
    cout << "f(0) = " << f(x) << endl;
    cout << "f(start) = " << f(start) << endl;*/

    // starting at the minimum:
    f.reset_n_eval();
    MinimizationResult minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(minres.fval, pow(opts.par_eps * step.get(vpids[0]), 2) * ndim);
    BOOST_CHECK_LT(f.get_n_eval(), 500);
    /*cout << "minimum found at: ";
    dump(cout, minres.values, vm);
    cout << endl;
    cout << "f = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;*/

    // starting somewhat off, with too large stepsize:
    for(size_t i=0; i<ndim; ++i){
        start.set(vpids[i], 1.0 + i);
        step.set(vpids[i], 10.0);
    }
    f.reset_n_eval();
    minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(minres.fval, pow(opts.par_eps * step.get(vpids[0]), 2) * ndim);
    BOOST_CHECK_LT(f.get_n_eval(), 500);
    /*cout << "minimum found at: ";
    dump(cout, minres.values, vm);
    cout << endl;
    cout << "f = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;*/
}

// 30 dimensions with different "scales", not aligned to axes
BOOST_AUTO_TEST_CASE(test30d_skewed){
    const size_t ndim = 10;
    Matrix hesse(ndim, ndim), covariance(ndim, ndim);
    // make the covariance matrix a 30*30 matrix by transforming the diagonal matrix of eigenvalues lambda with
    // a householder transform trafo: cov = trafo * lambda * trafo
    Matrix lambda(ndim, ndim), trafo(ndim, ndim);
    vector<double> v(ndim);
    auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    double v_norm2 = 0.0;
    for(size_t i=0; i < ndim; ++i){
        lambda(i, i) = pow(1.0 + i, 2);
        v[i] = rnd.uniform() - 0.5;
        v_norm2 += v[i] * v[i];
    }
    // construct trafo = I - 2 v v^T:
    //printf("trafo:\n");
    for(size_t i=0; i<ndim; ++i){
        for(size_t j=0; j<ndim; ++j){
            trafo(i,j) = (i==j?1:0) - 2 * v[i] * v[j] / v_norm2;
            //printf(" %10.3g", trafo(i,j));
        }
        //printf("\n");
    }
    // covariance_ij = trafo_ik * lambda_kl * trafo_lj, but lambda is diagonal, so
    // covariance_ij = trafo_ik * lambda_kk * trafo_kj
    for(size_t i=0; i<ndim; ++i){
        for(size_t j=0; j<ndim; ++j){
            for(size_t k=0; k<ndim; ++k){
                covariance(i,j) += trafo(i,k) * lambda(k,k) * trafo(k,j);
            }
        }
    }

    // test:
    VarIdManager vm;
    std::vector<ParId> vpids;
    ParIds pids;
    ParValues x0, start, step;
    map<ParId, pair<double, double> > ranges;
    for(size_t i=0; i<ndim; ++i){
        stringstream ss;
        ss << "p" << i;
        ParId pid = vm.create_par_id(ss.str());
        vpids.push_back(pid);
        x0.set(pid, 0.0);
        start.set(pid, 1.0);
        step.set(pid, 1.0);
        pids.insert(pid);
        ranges[pid] = make_pair(-inf, inf);
    }
    hesse = covariance;
    hesse.invert_cholesky();

    /*printf("covariance:\n");
    for(size_t i=0; i<ndim; ++i){
            for(size_t j=0; j<ndim; ++j){
                    printf(" %10.4g", covariance(i,j));
            }
            printf("\n");
    }
    printf("true hesse:\n");
    for(size_t i=0; i<ndim; ++i){
            for(size_t j=0; j<ndim; ++j){
                    printf(" %10.4g", hesse(i,j));
            }
            printf("\n");
    }*/

    quadratic_function f(pids, x0, hesse);

    newton_minimizer::options opts;
    //opts.debug = true;
    newton_minimizer min(opts);

    f.reset_n_eval();
    MinimizationResult minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(minres.fval, 1e-6);
    BOOST_CHECK_LT(f.get_n_eval(), 400);

    // large step sizes seem to be much les problematic than too small ones, test that:
    f.reset_n_eval();
    for(size_t i=0; i<ndim; ++i){
        step.set( vpids[i], 100.0);
    }
    minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(minres.fval, 1e-6);
    BOOST_CHECK_LT(f.get_n_eval(), 800);

    // test whether it's faster with correct covariance:
    /*f.reset_n_eval();
    minres = min.minimize2(MinimizationProblem(f, start, covariance, ranges));
    BOOST_CHECK_LT(minres.fval, 1e-10);
    BOOST_CHECK_LT(f.get_n_eval(), 100);*/
}

BOOST_AUTO_TEST_CASE(range1d){
    VarIdManager vm;
    ParId pid = vm.create_par_id("p0");
    ParValues x0;
    x0.set(pid, 0.0);
    ParIds pids;
    pids.insert(pid);

    Matrix hesse(1,1);
    hesse(0,0) = 1.0;

    quadratic_function f(pids, x0, hesse);
    newton_minimizer::options opts;
    newton_minimizer min(opts);

    // minimize a 1d function with a range at the minimum:
    ParValues start, step;
    start.set(pid, 1.0);
    step.set(pid, 1.0);
    std::map<ParId, pair<double, double> > ranges;
    ranges[pid] = make_pair(0, inf);
    MinimizationResult minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(minres.fval, 1e-10);
    BOOST_CHECK_LT(f.get_n_eval(), 20);

    // set range at 0.2
    f.reset_n_eval();
    ranges[pid] = make_pair(0.2, inf);
    minres = min.minimize(f, start, step, ranges);
    //cout << minres.values.get(pid) << endl;
    BOOST_CHECK_LT(minres.fval, 0.2 * 0.2 * (1 + 1e-7));
    BOOST_CHECK_LT(f.get_n_eval(), 50);
}


BOOST_AUTO_TEST_CASE(range10d){
    VarIdManager vm;
    std::vector<ParId> vpids;
    ParIds pids;
    const size_t ndim = 10;
    ParValues x0, start, step, x;
    Matrix hesse(ndim,ndim);
    map<ParId, pair<double, double> > ranges;
    for(size_t i=0; i<ndim; ++i){
        stringstream ss;
        ss << "p" << i;
        ParId pid = vm.create_par_id(ss.str());
        vpids.push_back(pid);
        x0.set(pid, 1.0);
        start.set(pid, 1.0);
        step.set(pid, 1.0);
        x.set(pid, 0.0);
        hesse(i,i) = 1.0;
        pids.insert(pid);
        ranges[pid] = make_pair(i==0?1.3:-inf, inf);
    }
    quadratic_function f(pids, x0, hesse);

    newton_minimizer::options opts;
    newton_minimizer min(opts);

    /*cout << "f(x0) = " << f(x0) << endl;
    cout << "f(0) = " << f(x) << endl;
    cout << "f(start) = " << f(start) << endl;*/

    // starting close to the minimum:
    f.reset_n_eval();
    start.set(vpids[0], 2.0);
    MinimizationResult minres = min.minimize(f, start, step, ranges);
    BOOST_CHECK_LT(fabs(minres.values.get(vpids[0]) - 1.3), 1e-7);
    BOOST_CHECK_LT(minres.fval, 0.3 * 0.3 + 1e-3);
    BOOST_CHECK_LT(f.get_n_eval(), 500);
}

BOOST_AUTO_TEST_SUITE_END()
