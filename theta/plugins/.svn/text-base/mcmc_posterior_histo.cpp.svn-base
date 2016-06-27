#include "plugins/mcmc_posterior_histo.hpp"
#include "interface/mcmc.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

using namespace theta;
using namespace std;

namespace{

//the result class for the metropolisHastings routine, saving the histograms
class MCMCPosteriorHistoResult: public MCMCResult{
    public:
        //ipar_ is the parameter of interest
        MCMCPosteriorHistoResult(const vector<size_t> & ipars_, size_t npar_, const vector<size_t> & nbins,
                                 const vector<double> & lower, const vector<double> & upper):
                           npar(npar_), ipars(ipars_){
            histos.resize(ipars_.size());
            for(size_t i=0; i<ipars_.size(); ++i){
                histos[i] = Histogram1D(nbins[i], lower[i], upper[i]);
            }
        }
        
        size_t getnpar() const{
            return npar;
        }
        
        //fill the parameters we are interested in into their histograms ...
        void fill(const double * x, double, size_t n){
            for(size_t i=0; i<ipars.size(); ++i){
                histos[i].fill(x[ipars[i]], n);
            }
        }
        
        const theta::Histogram1D & get_histo(size_t i) const{
            return histos[i];
        }
        
    private:
        size_t npar;
        vector<size_t> ipars;
        vector<theta::Histogram1D> histos;
};

// the result class for the smoothed version:
class MCMCPosteriorHistoResultSmoothed: public MCMCResult{
public:
    MCMCPosteriorHistoResultSmoothed(const vector<size_t> & ipars_, const vector<size_t> & nbins,
                                 const vector<double> & lower, const vector<double> & upper,
                                 const NLLikelihood & nll_): nll(nll_), npar(nll.getnpar()), ipars(ipars_){
            histos.resize(ipars.size());
            histos_tmp.resize(ipars.size());
            theta_assert(ipars.size()==nbins.size());
            theta_assert(ipars.size()==lower.size());
            theta_assert(ipars.size()==upper.size());
            for(size_t i=0; i<ipars.size(); ++i){
                histos[i] = Histogram1D(nbins[i], lower[i], upper[i]);
                histos_tmp[i] = Histogram1D(nbins[i], lower[i], upper[i]);
            }
    }
    
    void fill(const double * x, double nll0, size_t n){
        vector<double> myvalues(npar);
        if(std::isnan(nll0) || (std::isinf(nll0) && nll0 < 0)){
            throw range_error("nll0 is nan/-inf in mcmc_posterior_histo");
        }
        for(size_t ih=0; ih<histos.size(); ++ih){
            copy(x, x + npar, myvalues.begin());
            const double xmin = histos[ih].get_xmin();
            const double x_binwidth = (histos[ih].get_xmax() - histos[ih].get_xmin()) / histos[ih].get_nbins();
            for(size_t i=0; i<histos[ih].get_nbins(); ++i){
                myvalues[ipars[ih]] = xmin + (i + 0.5) * x_binwidth;
                double nll_value = nll(&myvalues[0]);
                if(std::isnan(nll_value) || (std::isinf(nll_value) && nll_value < 0)){
                    throw range_error("nll value is nan/-inf in mcmc_posterior_histo");
                }
                histos_tmp[ih].set(i, exp(-nll_value + nll0));
            }
            histos[ih].add_with_coeff(n / histos_tmp[ih].get_sum(), histos_tmp[ih]);
        }
    }
    
    size_t getnpar() const{
        return npar;
    }
    
    const Histogram1D & get_histo(size_t i) const{
        return histos[i];
    }
private:
    const NLLikelihood & nll;
    size_t npar;
    vector<size_t> ipars;
    vector<theta::Histogram1D> histos, histos_tmp;
};


}



void mcmc_posterior_histo::produce(const Data & data, const Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(!init){
        try{
            mcmc_strategy->init(model, override_parameter_distribution);
            //find ipars:
            ParIds nll_pars = nll->get_parameters();
            ipars.resize(parameters.size());
            for(size_t i=0; i<parameters.size(); ++i){
                for(ParIds::const_iterator it=nll_pars.begin(); it!=nll_pars.end(); ++it, ++ipars[i]){
                    if(*it == parameters[i]) break;
                }
            }
            init = true;
        }
        catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            throw invalid_argument(ex.message);
        }
    }
    
    if(!smooth){
        MCMCPosteriorHistoResult result(ipars, nll->getnpar(), nbins, lower, upper);
        mcmc_strategy->run_mcmc(*nll, result);
        for(size_t i=0; i<parameters.size(); ++i){
            products_sink->set_product(columns[i], result.get_histo(i));
        }
    }
    else{
        MCMCPosteriorHistoResultSmoothed result(ipars, nbins, lower, upper, *nll);
        mcmc_strategy->run_mcmc(*nll, result);
        for(size_t i=0; i<parameters.size(); ++i){
            products_sink->set_product(columns[i], result.get_histo(i));
        }
    }
}

mcmc_posterior_histo::mcmc_posterior_histo(const theta::Configuration & cfg): Producer(cfg), init(false), smooth(false){
    Setting s = cfg.setting;
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    size_t n = s["parameters"].size();
    for(size_t i=0; i<n; ++i){
        string parameter_name = s["parameters"][i];
        parameter_names.push_back(parameter_name);
        parameters.push_back(vm->get_par_id(parameter_name));
        nbins.push_back(static_cast<unsigned int>(s["histo_" + parameter_name]["nbins"]));
        lower.push_back(s["histo_" + parameter_name]["range"][0]);
        upper.push_back(s["histo_" + parameter_name]["range"][1]);
    }
    mcmc_strategy = construct_mcmc_strategy(cfg);
    if(s.exists("smooth")){
        smooth = s["smooth"];
    }
    for(size_t i=0; i<parameters.size(); ++i){
        columns.push_back(products_sink->declare_product(*this, "posterior_" + parameter_names[i], theta::typeHisto));
    }
}

REGISTER_PLUGIN(mcmc_posterior_histo)
