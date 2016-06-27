#include "plugins/mcmc_mean_prediction.hpp"
#include "interface/mcmc.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>

using namespace theta;
using namespace std;

//the result class for the metropolisHastings routine.
class MCMCMeanPredictionResult: public MCMCResult{
    public:
        MCMCMeanPredictionResult(const Model & model_, const ObsIds & observables, size_t npar_): model(model_), npar(npar_), obs_ids(observables), n(0){
            par_ids = model.get_parameters();
            theta_assert(par_ids.size() == npar);
            nll_min = std::numeric_limits<double>::infinity();
        }
        
        size_t getnpar() const{
            return npar;
        }
        
        //just save the parameter value we are interested in
        void fill(const double * x, double nll, size_t n_){
            n += n_;
            ParValues values(x, par_ids);
            //get the prediction using this values from the model:
            Data pred;
            model.get_prediction(pred, values);
            if(nll < nll_min){
               nll_min = nll;
               best = pred;
            }
            //add the prediction to sum, squaresum. Note
            // That if this is the first time, use "=" instead of "+=".
            for(ObsIds::const_iterator it=obs_ids.begin(); it!=obs_ids.end(); ++it){
               Histogram1D & h_pred = pred[*it];
               h_pred *= n_;
               Histogram1D & h_sum = sum[*it];
               Histogram1D & h_squaresum = squaresum[*it];
               if(h_sum.get_nbins()==0){
                   h_sum = h_pred;
                   h_pred *= h_pred;
                   h_pred *= 1.0 / n_;
                   h_squaresum = h_pred;
               }
               else{
                   h_sum += h_pred;
                   h_pred *= h_pred;
                   h_squaresum.add_with_coeff(1.0 / n_, h_pred);
               }
            }
        }

        void get_mean_width(Histogram1D & mean, Histogram1D & width, const ObsId & oid) const{
           const Histogram1D & h_sum = sum[oid];
           const Histogram1D & h_squaresum = squaresum[oid];
           mean = h_sum;
           mean *= 1.0 / n;
           width = Histogram1D(h_sum.get_nbins(), h_sum.get_xmin(), h_sum.get_xmax());
           for(size_t i=0; i<h_sum.get_nbins(); ++i){
              double ssum = h_squaresum.get(i);
              double sum = h_sum.get(i);
              double v = 1.0 / (n-1) * (ssum - 1.0 / n * sum * sum);
              width.set(i, sqrt(v));
           }
        }
        
        Histogram1D get_best(const ObsId & oid) const {
            return best[oid];
        }
        
    private:
        const Model & model;
        size_t npar;
        ObsIds obs_ids;
        //the total number of iterations seen:
        size_t n;
        ParIds par_ids;
        Data sum;
        Data squaresum;
        //information for the "best" point:
        double nll_min;
        Data best;
};


void mcmc_mean_prediction::produce(const Data & data, const Model & model) {
    if(!init){
        try{
            mcmc_strategy->init(model, override_parameter_distribution);
            init = true;
        }
        catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            invalid_argument(ex.message);
        }
    }
    
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    MCMCMeanPredictionResult result(model, observables, nll->getnpar());
    mcmc_strategy->run_mcmc(*nll, result);
    
    size_t i=0;
    for(ObsIds::const_iterator it=observables.begin(); it!=observables.end(); ++it, ++i){
        Histogram1D mean, width;
        result.get_mean_width(mean, width, *it);
        products_sink->set_product(c_mean[i], mean);
        products_sink->set_product(c_width[i], width);
        products_sink->set_product(c_best[i], result.get_best(*it));
    }
}

mcmc_mean_prediction::mcmc_mean_prediction(const theta::Configuration & cfg): Producer(cfg), init(false){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    Setting s = cfg.setting;
    size_t n = s["observables"].size();
    if(n==0){
        throw ConfigurationException("list 'observables' is empty");
    }
    for(size_t i=0; i<n; ++i){
        observables.insert(vm->get_obs_id(s["observables"][i]));
    }
    mcmc_strategy = construct_mcmc_strategy(cfg);
    for(ObsIds::const_iterator it=observables.begin(); it!=observables.end(); ++it){
        c_mean.push_back(products_sink->declare_product(*this, vm->get_name(*it) + "_mean", theta::typeHisto));
        c_width.push_back(products_sink->declare_product(*this, vm->get_name(*it) + "_width", theta::typeHisto));
        c_best.push_back(products_sink->declare_product(*this, vm->get_name(*it) + "_best", theta::typeHisto));
    }
}

REGISTER_PLUGIN(mcmc_mean_prediction)
