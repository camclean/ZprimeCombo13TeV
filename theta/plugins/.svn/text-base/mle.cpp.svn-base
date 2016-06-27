#include "plugins/mle.hpp"
#include "interface/log2_dot.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/minimizer.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>
#include <limits>

using namespace theta;
using namespace std;

namespace{

int get_index(const ParId & pid, const ParIds & pids){
   int result = 0;
   for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it, ++result){
      if((*it)==pid) return result;
   }
   throw invalid_argument("get_index: no such parameter");
}

}

void mle::produce(const theta::Data & data, const theta::Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(not start_step_ranges_init){
        const Distribution & d = nll->get_parameter_distribution();
        ranges.set_from(d);
        d.mode(start);
        try{
           step.set(asimov_likelihood_widths(model, override_parameter_distribution));
        }
        catch(const Exception & ex){
           throw logic_error("mle: could not initialize step size with asimov likelihood widths: " + ex.message);
        }
        start_step_ranges_init = true;
    }
    MinimizationResult minres = minimizer->minimize(*nll, start, step, ranges);
    products_sink->set_product(c_nll, minres.fval);
    for(size_t i=0; i<save_ids.size(); ++i){
        products_sink->set_product(parameter_columns[i], minres.values.get(save_ids[i]));
        if(minres.errors_plus.contains(save_ids[i])){
            products_sink->set_product(error_columns[i], 0.5 * (minres.errors_plus.get(save_ids[i]) + minres.errors_minus.get(save_ids[i])) );
        }
        else{
            products_sink->set_product(error_columns[i], NAN);
        }
    }
    if(write_prediction){
        Data pred;
        model.get_prediction(pred, minres.values);
        const ObsIds & obs = model.get_observables();
        for(ObsIds::const_iterator it=obs.begin(); it!=obs.end(); ++it){
            products_sink->set_product(c_pred[*it], pred[*it]);
        }
    }
    if(write_covariance){
       const size_t N = save_ids.size();
       Histogram1D h(N*N, 0, N*N);
       const ParIds & pars = nll->get_parameters();
       for(size_t i=0; i<N; ++i){
           int index_i = get_index(save_ids[i], pars);
           for(size_t j=0; j<N; ++j){
               int index_j = get_index(save_ids[j], pars);
               h.set(i*N + j, minres.covariance(index_i,index_j));
           }
       }
       products_sink->set_product(c_covariance, h);
    }
    if(write_ks_ts){
        const ObsIds & obs = model.get_observables();
        DataWithUncertainties pred;
        model.get_prediction(pred, minres.values);
        double ks_ts = 0.0;
        for(ObsIds::const_iterator it=obs.begin(); it!=obs.end(); ++it){
            const Histogram1D & data_o = data[*it];
            const Histogram1DWithUncertainties & pred_o = pred[*it];
            //data_o.check_compatibility(pred_o);
            double sum_d=0, sum_p=0;
            for(size_t i=0; i<data_o.get_nbins(); ++i){
                sum_d += data_o.get(i);
                sum_p += pred_o.get_value(i);
                ks_ts = max(ks_ts, fabs(sum_d - sum_p));
            }
        }
        products_sink->set_product(c_ks_ts, ks_ts);
    }
    if(write_pchi2){
        const ObsIds & obs = model.get_observables();
        Data pred;
        model.get_prediction(pred, minres.values);
        double pchi2 = 0.0;
        for(ObsIds::const_iterator it=obs.begin(); it!=obs.end(); ++it){
            const Histogram1D & data_o = data[*it];
            const Histogram1D & pred_o = pred[*it];
            pchi2 += template_pchisquare(data_o.get_data(), pred_o.get_data(), data_o.size());
        }
        products_sink->set_product(c_pchi2, pchi2);
    }
}

mle::mle(const theta::Configuration & cfg): Producer(cfg), start_step_ranges_init(false), write_covariance(false), write_ks_ts(false), write_pchi2(false), write_prediction(false){
    Setting s = cfg.setting;
    minimizer = PluginManager<Minimizer>::build(Configuration(cfg, s["minimizer"]));
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    size_t n_parameters = s["parameters"].size();
    for (size_t i = 0; i < n_parameters; i++) {
        string par_name = s["parameters"][i];
        save_ids.push_back(vm->get_par_id(par_name));
        parameter_names.push_back(par_name);
    }
    if(s.exists("write_covariance")){
        write_covariance = s["write_covariance"];
    }
    if(s.exists("write_ks_ts")){
        write_ks_ts = s["write_ks_ts"];
    }
    if(s.exists("write_pchi2")){
        write_pchi2 = s["write_pchi2"];
    }
    if(s.exists("write_prediction")){
        write_prediction = s["write_prediction"];
    }
    c_nll = products_sink->declare_product(*this, "nll", theta::typeDouble);
    for(size_t i=0; i<save_ids.size(); ++i){
        parameter_columns.push_back(products_sink->declare_product(*this, parameter_names[i], theta::typeDouble));
        error_columns.push_back(products_sink->declare_product(*this, parameter_names[i] + "_error", theta::typeDouble));
    }
    if(write_covariance){
       c_covariance = products_sink->declare_product(*this, "covariance", theta::typeHisto);
    }
    if(write_ks_ts){
       c_ks_ts = products_sink->declare_product(*this, "ks_ts", theta::typeDouble);
    }
    if(write_pchi2){
       c_pchi2 = products_sink->declare_product(*this, "pchi2", theta::typeDouble);
    }
    if(write_prediction){
        ObsIds oids = vm->get_all_observables();
        for(ObsIds::const_iterator oit=oids.begin(); oit!=oids.end(); ++oit){
            c_pred[*oit] = products_sink->declare_product(*this, "pred_" + vm->get_name(*oit), theta::typeHisto);
        }
    }
}

REGISTER_PLUGIN(mle)
