#include "plugins/mcmc_minimizer.hpp"
#include "interface/mcmc.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

using namespace theta;
using namespace std;

//the result class for the metropolisHastings routine.
class MCMCMinResult: public ResultMeanCov{
    public:
        //ipar_ is the parameter of interest
        MCMCMinResult(size_t npar, const ParIds & pids_): ResultMeanCov(npar), n(0), n_different(0), pids(pids_), min_nll(numeric_limits<double>::infinity()){
        }
                
        virtual void fill2(const double * x, double nll, size_t n_){
            ++n_different;
            n += n_;
            if(nll < min_nll){
                min_nll = nll;
                values_at_min.set(ParValues(x, pids));
            }
        }
        
        const ParValues & values_at_minimum() const{
            return values_at_min;
        }
        
        double get_min_nll() const{
            return min_nll;
        }
        
        double acceptance_rate() const{
            return 1.0 * n_different / n;
        }
        
    private:
        size_t n, n_different;
        ParIds pids;
        double min_nll;
        ParValues values_at_min;
};


MinimizationResult mcmc_minimizer::minimize(const Function & f, const ParValues & start, const ParValues & step, const Ranges & ranges){
    const ParIds & pids = f.get_parameters();
    size_t n = pids.size();
    Matrix sqrt_cov(n,n);
    options.startvalues.resize(n);
    ParIds::const_iterator it=pids.begin();
    for(size_t i=0; i<n; ++it, ++i){
        options.startvalues[i] = start.get(*it);
    }
    it=pids.begin();
    for(size_t i=0; i<n; ++it, ++i){
        sqrt_cov(i,i) = step.get(*it);
        theta_assert(std::isfinite(sqrt_cov(i,i)));
    }
    for (int i = 0; i < bootstrap_mcmcpars; i++) {
        ResultMeanCov res(n);
        metropolis_hastings_multigauss(f, options, res, *rnd_gen, sqrt_cov);
        if(res.getCount() < options.iterations / 2) {
            throw MinimizationException("during mcmcpars bootstrapping: more than 50% of the chain is infinite!");
        }
        options.startvalues = res.getMeans();
        Matrix cov = res.getCov();
        get_cholesky(cov, sqrt_cov);
    }
    
    MCMCMinResult result(n, pids);
    double first_nll = f(&options.startvalues[0]);
    if(!std::isfinite(first_nll)){
        throw MinimizationException("first nll value was not finite");
    }
    metropolis_hastings_multigauss(f, options, result, *rnd_gen, sqrt_cov);
    //now step 2: call the after_minimizer:
    const ParValues & start2 = result.values_at_minimum();
    ParValues step2;
    Matrix cov = result.getCov();
    it=pids.begin();
    for(size_t i=0; i<n; ++it, ++i){
        step2.set(*it, sqrt(cov(i, i)));
    }
    if(after_minimizer.get()){
        return after_minimizer->minimize(f, start2, step2, ranges);
    }
    else{
        MinimizationResult res;
        res.fval = result.get_min_nll();
        res.values.set(start2);
        res.errors_minus.set(step2);
        res.errors_plus.set(step2);
        res.covariance = cov;
        return res;
    }
}

mcmc_minimizer::mcmc_minimizer(const theta::Configuration & cfg): RandomConsumer(cfg, cfg.setting["name"]), name(cfg.setting["name"]), bootstrap_mcmcpars(0){
    Setting s = cfg.setting;
    options.ignore_inf_nll_start = true;
    options.iterations = s["iterations"];
    if(s.exists("burn-in")){
        options.burn_in = s["burn-in"];
    }
    else{
        options.burn_in = options.iterations / 10;
    }
    if(s.exists("after_minimizer")){
        after_minimizer = PluginManager<Minimizer>::build(Configuration(cfg, s["after_minimizer"]));
    }
    if(s.exists("stepsize_factor")){
        options.factor = s["stepsize_factor"];
    }
    if(s.exists("bootstrap_mcmcpars")){
        bootstrap_mcmcpars = s["bootstrap_mcmcpars"];
    }
}

REGISTER_PLUGIN(mcmc_minimizer)
