#include "root/root_minuit.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"

#include "TError.h"

using namespace theta;
using namespace std;

// function adapter to be used by root minuit and
// "redirects" function calls to theta::Function
class RootMinuitFunctionAdapter: public ROOT::Math::IMultiGenFunction{
public:
    
    virtual ROOT::Math::IBaseFunctionMultiDim*	Clone() const{
        throw Exception("RootMinuitFunctionAdapter::Clone not implemented");
    }

    virtual unsigned int NDim() const{
        return ndim;
    }

    RootMinuitFunctionAdapter(const Function & f_): f(f_), ndim(f.get_parameters().size()){
    }

    virtual double DoEval(const double * x) const{
        for(size_t i=0; i<ndim; ++i){
            if(std::isnan(x[i])){
               throw MinimizationException("minuit called likelihood function with NAN argument!");
            }
        }
        double result = f(x);
        if(std::isinf(result)){
           throw MinimizationException("function to minimize was infinity during minimization");
        }
        return result;
    }

private:
    const theta::Function & f;
    const size_t ndim;
};


MinimizationResult root_minuit::minimize(const theta::Function & f, const theta::ParValues & start,
        const theta::ParValues & steps, const Ranges & ranges){
    //I would like to re-use min. However, it horribly fails after very few uses with
    // unsigned int ROOT::Minuit2::MnUserTransformation::IntOfExt(unsigned int) const: Assertion `!fParameters[ext].IsFixed()' failed.
    // when calling SetFixedVariable(...).
    //Using a "new" one every time seems very wastefull, but it seems to work ...
    std::auto_ptr<ROOT::Minuit2::Minuit2Minimizer> min(new ROOT::Minuit2::Minuit2Minimizer(type));
    //min->SetPrintLevel(0);
    if(max_function_calls > 0) min->SetMaxFunctionCalls(max_function_calls);
    if(max_iterations > 0) min->SetMaxIterations(max_iterations);
    MinimizationResult result;

    //1. setup parameters, limits and initial step sizes
    ParIds parameters = f.get_parameters();
    int ivar=0;
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++ivar){
        pair<double, double> range = ranges.get(*it);
        double def = start.get(*it);
        double step = steps.get(*it);
        stringstream ss;
        ss << "par" << ivar;
        string name = ss.str();
        //use not the ranges directly, but a somewhat more narrow range (one permille of the respective border)
        // in order to avoid that the numerical evaluation of the numerical derivative at the boundaries pass these
        // boundaries ...
        if(step == 0.0){
            min->SetFixedVariable(ivar, name, def);
        }
        else if(std::isinf(range.first)){
            if(std::isinf(range.second)){
                min->SetVariable(ivar, name, def, step);
            }
            else{
                min->SetUpperLimitedVariable(ivar, name, def, step, range.second - fabs(range.second) * 0.001);
            }
        }
        else{
            if(std::isinf(range.second)){
                min->SetLowerLimitedVariable(ivar, name, def, step, range.first + fabs(range.first) * 0.001);
            }
            else{ // both ends are finite
                if(range.first==range.second){
                    min->SetFixedVariable(ivar, name, range.first);
                }
                else{
                    min->SetLimitedVariable(ivar, name, def, step, range.first + fabs(range.first) * 0.001, range.second - fabs(range.second) * 0.001);
                }
            }
        }
    }

    //2. setup the function
    RootMinuitFunctionAdapter minuit_f(f);
    min->SetFunction(minuit_f);

    //3. setup tolerance
    double root_tol = min->Tolerance();
    if(root_tol == 0.01){
        root_tol *= 0.1;
    }
    min->SetTolerance(tolerance_factor * root_tol);
    //3.a. error definition. Unfortunately, SetErrorDef in ROOT is not documented, so I had to guess.
    // 0.5 seems to work somehow.
    min->SetErrorDef(0.5);
    
    //4. minimize. In case of failure, try harder. Discard all output generated in min->Minimize.
    bool success;
    success = min->Minimize();
    if(!success){
        for(int i=1; i<=n_retries; i++){
            success = min->Minimize();
            if(success) break;
        }
    }

    //5. do error handling
    if(not success){
        int status = min->Status();
        int status_1 = status % 10;
        //int status_2 = status / 10;
        stringstream s;
        s << "MINUIT returned status " << status;
        switch(status_1){
            case 1: s << " (Covariance was made pos defined)"; break;
            case 2: s << " (Hesse is invalid)"; break;
            case 3: s << " (Edm is above max)"; break;
            case 4: s << " (Reached call limit)"; break;
            case 5: s << " (Some other failure)"; break;
            default:
                s << " [unexpected status code]";
        }
        throw MinimizationException(s.str());
    }

    //6. convert result
    result.fval = min->MinValue();
    ivar = 0;
    const double * x = min->X();
    const double * errors = 0;
    bool have_errors = min->ProvidesError();
    if(have_errors) errors = min->Errors();
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++ivar){
        if(!std::isfinite(x[ivar])){
             throw MinimizationException("minuit: minimum is at non-finite point");
        }
        result.values.set(*it, x[ivar]);
        if(have_errors){
            result.errors_plus.set(*it, errors[ivar]);
            result.errors_minus.set(*it, errors[ivar]);
        }
        else{
            result.errors_plus.set(*it, -1);
            result.errors_minus.set(*it, -1);
        }
    }
    result.covariance.reset(parameters.size(), parameters.size());
    //I would use min->CovMatrixStatus here to check the validity of the covariance matrix,
    // if only it was documented ...
    if(min->ProvidesError()){
        for(size_t i=0; i<parameters.size(); ++i){
            for(size_t j=0; j<parameters.size(); ++j){
                result.covariance(i,j) = min->CovMatrix(i,j);
            }
        }
    }
    else{
        for(size_t i=0; i<parameters.size(); ++i){
            result.covariance(i,i) = -1;
        }
    }
    return result;
}

root_minuit::root_minuit(const Configuration & cfg): tolerance_factor(1),
        max_iterations(0), max_function_calls(0), n_retries(2) {
    gErrorIgnoreLevel = kFatal + 1;
    if(cfg.setting.exists("max_iterations")){
        max_iterations = cfg.setting["max_iterations"];
    }
    if(cfg.setting.exists("max_function_calls")){
        max_iterations = cfg.setting["max_function_calls"];
    }
    if(cfg.setting.exists("n_retries")){
        n_retries = cfg.setting["n_retries"];
    }
    string method = "migrad";
    if(cfg.setting.exists("method")){
        method = (string)cfg.setting["method"];
    }
    if(method=="migrad"){
        type = ROOT::Minuit2::kMigrad;
    }
    else if(method == "simplex"){
        type = ROOT::Minuit2::kSimplex;
    }
    else{
        stringstream s;
        s << "invalid method '" << method << "' (allowed are only 'migrad' and 'simplex')";
        throw ConfigurationException(s.str());
    }       
    if(cfg.setting.exists("tolerance_factor")){
        tolerance_factor = cfg.setting["tolerance_factor"];
        if(tolerance_factor <= 0){
            throw ConfigurationException("tolerance_factor is <= 0.0; this is not allowed");
        }
    }
}

REGISTER_PLUGIN(root_minuit)

