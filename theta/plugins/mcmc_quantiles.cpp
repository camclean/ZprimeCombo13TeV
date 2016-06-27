#include "plugins/mcmc_quantiles.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>
#include <iomanip>

using namespace theta;
using namespace std;

//the result class for the metropolisHastings routine.
class MCMCPosteriorQuantilesResult: public MCMCResult{
    public:
        //ipar_ is the parameter of interest
        MCMCPosteriorQuantilesResult(size_t npar_, size_t ipar_, size_t iterations): npar(npar_), ipar(ipar_), n_iterations(iterations),  n_iterations_total(0), n_iterations_different(0){
            par_values.reserve(iterations);
        }
        
        virtual size_t getnpar() const{
            return npar;
        }
        
        //just save the parameter value we are interested in
        virtual void fill(const double * x, double, size_t n_){
            n_iterations_total += n_;
            ++n_iterations_different;
            for(size_t i=0; i<n_; ++i){
               par_values.push_back(x[ipar]);
            }
        }

        double get_accrate() const{
            return 1.0 * n_iterations_different / n_iterations_total;
        }
        
        //return the quantile q
        double get_quantile(double q){
            if(par_values.size()!=n_iterations){
                throw invalid_argument("MCMCPosteriorQuantilesResult: called get_quantile before chain has finished!");
            }
            int index = static_cast<int>(q * n_iterations);
            if(index >= static_cast<int>(n_iterations)) index = n_iterations-1;
            std::nth_element(par_values.begin(), par_values.begin() + index, par_values.end());
            return par_values[index];
        }
        
    private:
        size_t npar;
        size_t ipar;
        const size_t n_iterations;
        size_t n_iterations_total;
        size_t n_iterations_different;
        vector<double> par_values;
};

void mcmc_quantiles::produce(const Data & data, const Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(!init || (re_init > 0 && itoy % re_init == 0)){
        try{
            mcmc_strategy->init(model, override_parameter_distribution);
            //find the number of the parameter of interest:
            const ParIds & pars = nll->get_parameters();
            size_t i = 0;
            for(ParIds::const_iterator it=pars.begin(); it!=pars.end(); ++it, ++i){
                if(*it == par_id){
                    ipar = i;
                }
            }
            init = true;
        }
        catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            throw invalid_argument(ex.message);
        }
    }
    ++itoy;
    MCMCPosteriorQuantilesResult result(nll->getnpar(), ipar, mcmc_strategy->get_options().iterations);
    mcmc_strategy->run_mcmc(*nll, result);
    for(size_t i=0; i<quantiles.size(); ++i){
        products_sink->set_product(columns[i], result.get_quantile(quantiles[i]));
    }
    products_sink->set_product(c_accrate, result.get_accrate());
}


mcmc_quantiles::mcmc_quantiles(const theta::Configuration & cfg): Producer(cfg),
   init(false), par_id(get_parameter(cfg, "parameter")), re_init(0), itoy(0) {
    Setting s = cfg.setting;
    string parameter = s["parameter"];
    size_t n = s["quantiles"].size();
    if(n==0){
        throw ConfigurationException("mcmc_quantiles: list of requested quantiles is empty");
    }
    quantiles.reserve(n);
    for(size_t i=0; i<n; ++i){
        quantiles.push_back(s["quantiles"][i]);
        if(quantiles[i]<=0.0 || quantiles[i]>=1.0){
            throw ConfigurationException("mcmc_quantiles: quantiles out of range (must be strictly between 0 and 1)");
        }
    }
    if(s.exists("re-init")){
        re_init = s["re-init"];
    }
    mcmc_strategy = construct_mcmc_strategy(cfg);
    for(size_t i=0; i<quantiles.size(); ++i){
        stringstream ss;
        ss << "quant" << setw(5) << setfill('0') << static_cast<int>(quantiles[i] * 10000 + 0.5);
        columns.push_back(products_sink->declare_product(*this, ss.str(), theta::typeDouble));
    }
    c_accrate = products_sink->declare_product(*this, "accrate", theta::typeDouble);
}

REGISTER_PLUGIN(mcmc_quantiles)

