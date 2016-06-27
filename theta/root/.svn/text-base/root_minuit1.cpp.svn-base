#include "root/root_minuit1.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/redirect_stdio.hpp"

#include "interface/exception.hpp"
#include "TMinuit.h"

#include <iomanip>

using namespace theta;
using namespace std;

namespace{
    
// for TMinuit, function evaluation is provided by deriving from TMinuit and overloading Eval.
class MyTMinuit: public TMinuit{
private:
    const theta::Function & f;
    size_t ndim;
    vector<double> min_;
    double f_at_min;
    
public:
    // note 2*n+1 for TMinuit
    explicit MyTMinuit(const Function & f_): TMinuit(2*f_.get_parameters().size()+1), f(f_), ndim(f.get_parameters().size()){
        min_.resize(ndim);
    }
    
    const vector<double> & get_min() const{
        return min_;
    }
    
    double get_f_at_min() const{
        return f_at_min;
    }
    
    // see http://root.cern.ch/root/html/TMinuit.html#TMinuit:Eval
    virtual Int_t Eval(Int_t npar, Double_t * grad, Double_t & fval, Double_t * par, Int_t flag){
        for(size_t i=0; i<ndim; ++i){
            if(std::isnan(par[i])){
                throw MinimizationException("minuit called likelihood function with NAN argument!");
            }
        }
        fval = f(par);
        /*out << " flag = " << flag;
        out << " x = ";
        for(size_t i=0; i<ndim; ++i){
            out << par[i] << " ";
        }
        out << "--> " << setprecision(10) << fval << endl;*/
        if(std::isinf(fval)){
            throw MinimizationException("function to minimize was infinity during minimization");
        }
        return 0;
    }
};

}


MinimizationResult root_minuit1::minimize(const theta::Function & f, const theta::ParValues & start,
        const theta::ParValues & steps, const Ranges & ranges){
    MyTMinuit min(f);
    Double_t args[10];
    Int_t res;
    args[0] = -1;
    min.mnexcm("SET PRINT", args, 1, res);
    theta_assert(res==0);
    min.mnexcm("SET NOW", args, 0, res);
    theta_assert(res==0);
    //1. setup parameters, limits and initial step sizes
    ParIds parameters = f.get_parameters();
    const size_t n = parameters.size();
    int ivar=0;
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++ivar){
        pair<double, double> range = ranges.get(*it);
        double def = start.get(*it);
        double step = steps.get(*it);
        stringstream ss;
        ss << "par" << ivar;
        std::string sname = ss.str();
        const char * name = sname.c_str();
        //use not the ranges directly, but a somewhat more narrow range (one permille of the respective border)
        // in order to avoid that the numerical evaluation of the numerical derivative at the boundaries pass these
        // boundaries ...
        if(step == 0.0){
            min.DefineParameter(ivar, name, def, 0.0, def, def);
            min.FixParameter(ivar);
        }
        else if(std::isinf(range.first)){
            if(std::isinf(range.second)){
                // non-bounded parameters have range (0.0, 0.0):
                min.DefineParameter(ivar, name, def, step, 0.0, 0.0);
            }
            else{
                min.DefineParameter(ivar, name, def, step, def - infinity * step, range.second - fabs(range.second) * 0.001);
            }
        }
        else{
            if(std::isinf(range.second)){
                min.DefineParameter(ivar, name, def, step, range.first + fabs(range.first) * 0.001, def + infinity * step);
            }
            else{ // both ends are finite
                theta_assert(range.first < range.second);
                min.DefineParameter(ivar, name, def, step, range.first + fabs(range.first) * 0.001, range.second - fabs(range.second) * 0.001);
            }
        }
    }
    
    // error definition is: change of 0.5 corresponds to 1 sigma; we have negative log-likelihood (not 2*nll / chi2):
    args[0] = 0.5;
    min.mnexcm("SET ERR", args, 1, res);
    theta_assert(res==0);
    
    if(strategy >= 0){
        args[0] = strategy;
        min.mnexcm("SET STR", args, 1, res);
        theta_assert(res==0);
    }
    min.mnexcm("MIG", args, 0, res);
    if(res!=0){
        stringstream s;
        s << "MIG returned " << res;
        throw MinimizationException(s.str());
    }
    if(hesse){
        min.mnexcm("HES", args, 0, res);
        if(res!=0){
            stringstream s;
            s << "HES returned " << res;
            throw MinimizationException(s.str());
        }
    }
    
    MinimizationResult result;
    double fedm, errdef;
    int idummy1, idummy2;
    min.mnstat(result.fval, fedm, errdef, idummy1, idummy2, res);
    //theta_assert(errdef == 0.5);
    vector<double> xmin(n);
    ivar = 0;
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++ivar){
        TString name;
        double val, err, ddummy1, ddummy2;
        min.mnpout(ivar, name, val, err, ddummy1, ddummy2, idummy1);
        result.values.set(*it, val);
        result.errors_plus.set(*it, err);
        result.errors_minus.set(*it, err);
    }
    
    // res is the accuracy of the covariance; only 3 means it's accurate.
    vector<double> emat(n*n, -1);
    if(res>=3){
        min.mnemat(&emat[0], n);
    }// else: leave matrix at -1
    result.covariance.reset(n, n);
    for(size_t i=0; i<parameters.size(); ++i){
        for(size_t j=0; j<parameters.size(); ++j){
            result.covariance(i,j) = emat[i*n+j];
        }
    }
    return result;
}

root_minuit1::root_minuit1(const Configuration & cfg): tolerance(-1), infinity(1e5), strategy(-1), n_retries(2), hesse(false) {
    if(cfg.setting.exists("n_retries")){
        n_retries = cfg.setting["n_retries"];
    }
    if(cfg.setting.exists("infinity")){
        infinity = cfg.setting["infinity"];
    }
    if(cfg.setting.exists("strategy")){
        strategy = cfg.setting["strategy"];
        if(strategy < 0 or strategy > 2){
            throw ConfigurationException("strategy must be between 0 and 2");
        }
    }
    if(cfg.setting.exists("hesse")){
        hesse = cfg.setting["hesse"];
    }
    if(cfg.setting.exists("tolerance")){
        tolerance = cfg.setting["tolerance"];
        if(tolerance <= 0){
            throw ConfigurationException("tolerance <= 0.0 not allowed");
        }
    }
}

REGISTER_PLUGIN(root_minuit1)

