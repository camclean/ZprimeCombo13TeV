#include "plugins/mcmc_posterior_ratio.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

using namespace theta;
using namespace std;

//the result class for the metropolisHastings routine.
class MCMCPosteriorRatioResult: public MCMCResult{
    public:
        MCMCPosteriorRatioResult(size_t npar_): npar(npar_), min_nll_value(numeric_limits<double>::infinity()), n_total(0){}
        
        size_t getnpar() const{
            return npar;
        }
        
        void fill(const double *, double nll, size_t n_){
            nll_values.push_back(nll);
            if(nll < min_nll_value) min_nll_value = nll;
            n.push_back(n_);
            n_total += n_;
        }
        
        //return the negative logarithm of the average posterior
        double get_nl_average_posterior(){
            double posterior_sum = 0.0;
            //instead of calculating
            // - log (   1/N  * sum_{i=1}^N exp (-nll_i)    )
            // calculate
            // - log (  exp(-min_nll) * 1/N * sum_{i=1}^N  exp(-nll_i + min_nll)   ) = min_nll - log (  1/N * sum_{i=1}^N  exp (-nll_i + min_nll)   )
            // which is the same , but numerically much better
            for(size_t i=0; i<n.size(); ++i){
                posterior_sum += n[i] * exp(min_nll_value - nll_values[i]);
            }
            return min_nll_value - log(posterior_sum / n_total);
        }
        
    private:
        size_t npar;
        vector<double> nll_values;
        double min_nll_value;
        vector<size_t> n;
        size_t n_total;
};


void mcmc_posterior_ratio::produce(const theta::Data & data, const theta::Model & model) {
    if(!init){
        try{
            const ParIds & pids = model.get_parameters();
            ParValues sb_mode, b_mode;
            s_plus_b->mode(sb_mode);
            b_only->mode(b_mode);
            const size_t n = pids.size();
            ParIds::const_iterator it=pids.begin();
            options_sb.startvalues.resize(n);
            options_b.startvalues.resize(n);
            for(size_t i=0; i<n; ++it, ++i){
                options_sb.startvalues[i] = sb_mode.get(*it);
                options_b.startvalues[i] = b_mode.get(*it);
            }
            sqrt_cov_sb = get_sqrt_cov2(*rnd_gen, model, s_plus_b);
            sqrt_cov_b = get_sqrt_cov2(*rnd_gen, model, b_only);
            init = true;
        }catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            throw invalid_argument(ex.message);
        }
    }
    
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    
    //a. calculate s plus b:
    MCMCPosteriorRatioResult res_sb(nll->getnpar());
    metropolis_hastings_multigauss(*nll, options_sb, res_sb, *rnd_gen, sqrt_cov_sb);
    double nl_posterior_sb = res_sb.get_nl_average_posterior();

    //b. calculate b only:
    MCMCPosteriorRatioResult res_b(nll->getnpar());
    metropolis_hastings_multigauss(*nll, options_b, res_b, *rnd_gen, sqrt_cov_b);
    double nl_posterior_b = res_b.get_nl_average_posterior();

    if(std::isnan(nl_posterior_sb) || std::isnan(nl_posterior_b)){
        throw Exception("average posterior was NAN");
    }
    products_sink->set_product(c_nl_posterior_sb, nl_posterior_sb);
    products_sink->set_product(c_nl_posterior_b, nl_posterior_b);
}

mcmc_posterior_ratio::mcmc_posterior_ratio(const theta::Configuration & cfg): Producer(cfg), RandomConsumer(cfg, get_name()), init(false){
    Setting s = cfg.setting;
    s_plus_b = theta::PluginManager<Distribution>::build(theta::Configuration(cfg, s["signal-plus-background-distribution"]));
    b_only = theta::PluginManager<Distribution>::build(theta::Configuration(cfg, s["background-only-distribution"]));    
    options_sb.iterations = options_b.iterations = s["iterations"];
    if(s.exists("burn-in")){
        options_sb.burn_in = options_b.burn_in = s["burn-in"];
    }
    else{
        options_sb.burn_in = options_b.burn_in = options_sb.iterations / 10;
    }
    c_nl_posterior_sb = products_sink->declare_product(*this, "nl_posterior_sb", theta::typeDouble);
    c_nl_posterior_b =  products_sink->declare_product(*this, "nl_posterior_b",  theta::typeDouble);
}

REGISTER_PLUGIN(mcmc_posterior_ratio)
