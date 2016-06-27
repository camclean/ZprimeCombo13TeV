#include "interface/model.hpp"
#include "interface/histogram-function.hpp"
#include "interface/distribution.hpp"
#include "interface/plugin.tcc"
#include "interface/log2_dot.hpp"

#include <limits>

using namespace std;
using namespace theta;

REGISTER_PLUGIN_BASETYPE(Model);

const ParIds & Model::get_parameters() const{
    return parameters;
}

const ParIds & Model::get_rvobservables() const{
    return rvobservables;
}

const ObsIds & Model::get_observables() const{
    return observables;
}


atomic_int theta::n_nll_eval, theta::n_nll_eval_with_derivative;

namespace{

struct n_nll_eval_reset{
    n_nll_eval_reset(){
        atomic_set(&n_nll_eval, 0);
        atomic_set(&n_nll_eval_with_derivative, 0);
    }
} resetter;

    
/* default_model */
class default_model_nll: public NLLikelihood{
friend class theta::default_model;
public:
    using Function::operator();
    virtual double operator()(const ParValues & values) const;
    
    virtual void set_override_distribution(const boost::shared_ptr<Distribution> & d);
    virtual const Distribution & get_parameter_distribution() const{
        if(override_distribution) return *override_distribution;
        else return model.get_parameter_distribution();
    }
    virtual double eval_with_derivative(const ParValues & v, ParValues & der) const;
        
protected:
    const default_model & model;
    const Data & data;
    bool robust_nll;
    boost::shared_ptr<Distribution> override_distribution;
    
    default_model_nll(const default_model & m, const Data & data);
    
protected:
    //cached predictions:
    mutable Data predictions;
    mutable std::map<ParId, Histogram1D> der_pred;
};

// includes additive Barlow-Beeston uncertainties, where the extra nuisance parameters of this method (1 per bin) have been "profiled out".
class default_model_bbadd_nll: public default_model_nll {
friend class theta::default_model;
public:
    using Function::operator();
    virtual double operator()(const ParValues & values) const;
    virtual double eval_with_derivative(const ParValues & values, ParValues & der) const;
    
private:
    default_model_bbadd_nll(const default_model & m, const Data & data);
    mutable DataWithUncertainties predictions_wu;
    mutable std::map<ParId, Histogram1DWithUncertainties> der_pred_wu;
};

} // anon. namespace


void default_model::set_prediction(const ObsId & obs_id, boost::ptr_vector<Function> & coeffs_, boost::ptr_vector<HistogramFunction> & histos_){
    observables.insert(obs_id);
    const size_t n = coeffs_.size();
    theta_assert(n > 0);
    theta_assert(n==histos_.size());
    for(std::vector<std::pair<ObsId, size_t> >::const_iterator it=last_indices.begin(); it!=last_indices.end(); ++it){
        theta_assert(it->first != obs_id);
    }
    const size_t imin = hfs.size();
    coeffs.transfer(coeffs.end(), coeffs_.begin(), coeffs_.end(), coeffs_);
    hfs.transfer(hfs.end(), histos_.begin(), histos_.end(), histos_);
    last_indices.push_back(make_pair(obs_id, hfs.size()));
    const size_t imax = hfs.size();
    size_t nbins = 0;
    double xmin = NAN, xmax = NAN;
    bool first = true;
    ParIds obs_pids;
    for(size_t i = imin; i < imax; ++i){
        if(first){
            hfs[i].get_histogram_dimensions(nbins, xmin, xmax);
            histo_dimensions.push_back(hdim());
            histo_dimensions.back().xmin = xmin;
            histo_dimensions.back().xmax = xmax;
            histo_dimensions.back().nbins = nbins;
            first = false;
        }
        else{
            size_t nbins_tmp = 0;
            double xmin_tmp = NAN, xmax_tmp = NAN;
            hfs[i].get_histogram_dimensions(nbins_tmp, xmin_tmp, xmax_tmp);
            if(nbins!=nbins_tmp || xmin!=xmin_tmp || xmax!=xmax_tmp){
                throw invalid_argument("default_model::set_prediction: histogram dimensions mismatch");
            }
        }
        parameters.insert_all(hfs[i].get_parameters());
        parameters.insert_all(coeffs[i].get_parameters());
        obs_pids.insert_all(hfs[i].get_parameters());
        obs_pids.insert_all(coeffs[i].get_parameters());
    }
    obs_parameters.push_back(obs_pids);
    theta_assert(histo_dimensions.size() == last_indices.size());
    theta_assert(histo_dimensions.size() == obs_parameters.size());
}

template<typename HT>
void default_model::get_prediction_impl(DataT<HT> & result, const ParValues & parameters) const{
    size_t imin = 0;
    size_t nobs = last_indices.size();
    for(size_t i=0; i<nobs; ++i){
        const ObsId & oid = last_indices[i].first;
        result[oid].reset(histo_dimensions[i].nbins, histo_dimensions[i].xmin, histo_dimensions[i].xmax);
        for(; imin < last_indices[i].second; ++imin){
            hfs[imin].add_with_coeff_to(result[oid], coeffs[imin](parameters), parameters);
        }
    }
}

void default_model::get_prediction(DataWithUncertainties & result, const ParValues & parameters) const {
    get_prediction_impl<Histogram1DWithUncertainties>(result, parameters);
}

void default_model::get_prediction(Data & result, const ParValues & parameters) const {
    get_prediction_impl<Histogram1D>(result, parameters);
}

void default_model::get_prediction_with_derivative(const ObsId & oid, Histogram1D & h, std::map<ParId, Histogram1D> & der, const ParValues & parameters) const{
    size_t iobs = 0;
    while(last_indices[iobs].first != oid) ++iobs;
    get_prediction_with_derivative(iobs, h, der, parameters);
}

void default_model::get_prediction_with_derivative(size_t iobs, Histogram1D & h, std::map<ParId, Histogram1D> & der, const ParValues & parameters) const{
    const size_t imin = iobs == 0? 0 : last_indices[iobs-1].second;
    const size_t imax = last_indices[iobs].second;
    h.reset(histo_dimensions[iobs].nbins, histo_dimensions[iobs].xmin, histo_dimensions[iobs].xmax);
    for(ParIds::const_iterator pit=this->parameters.begin(); pit!=this->parameters.end(); ++pit){
        der[*pit].reset(histo_dimensions[iobs].nbins, histo_dimensions[iobs].xmin, histo_dimensions[iobs].xmax);
    }
    ParValues coeff_ders(parameters);
    for(size_t i=imin; i<imax; ++i){
        // the model prediction for the observable is sum_i  coeff_i * hf_i, so its derivative w.r.t. p
        // is sum_i dcoeff_i/dp * hf_i + coeff_i * dhf_i / dp
        Histogram1D hi(histo_dimensions[iobs].nbins, histo_dimensions[iobs].xmin, histo_dimensions[iobs].xmax);
        coeff_ders.clear();
        double coeff = coeffs[i].eval_with_derivative(parameters, coeff_ders);
        
        // coeff_i * dhf_i / dp part of derivative:
        hfs[i].eval_and_add_derivatives(hi, der, coeff, parameters);
        // update result itself:
        h.add_with_coeff(coeff, hi);
        
        // add the dcoeff_i/dp * hf_i part of derivative:
        const ParIds & coeff_pids = coeffs[i].get_parameters();
        for(ParIds::const_iterator pit=coeff_pids.begin(); pit!=coeff_pids.end(); ++pit){
            der[*pit].add_with_coeff(coeff_ders.get_unchecked(*pit), hi);
        }
    }
}

void default_model::get_prediction_with_derivative(const ObsId & oid, Histogram1DWithUncertainties & h, std::map<ParId, Histogram1DWithUncertainties> & der, const ParValues & parameters) const{
    size_t iobs = 0;
    while(last_indices[iobs].first != oid) ++iobs;
    get_prediction_with_derivative(iobs, h, der, parameters);
}

void default_model::get_prediction_with_derivative(size_t iobs, Histogram1DWithUncertainties & h,
                                            std::map<ParId, Histogram1DWithUncertainties > & der, const ParValues & parameters) const{
    const size_t imin = iobs == 0? 0 : last_indices[iobs-1].second;
    const size_t imax = last_indices[iobs].second;
    h.reset(histo_dimensions[iobs].nbins, histo_dimensions[iobs].xmin, histo_dimensions[iobs].xmax);
    for(ParIds::const_iterator pit=this->parameters.begin(); pit!=this->parameters.end(); ++pit){
        der[*pit].reset(histo_dimensions[iobs].nbins, histo_dimensions[iobs].xmin, histo_dimensions[iobs].xmax);
    }
    ParValues coeff_ders(parameters);
    for(size_t i=imin; i<imax; ++i){
        // the model prediction for the observable is sum_i  coeff_i * hf_i, so its derivative w.r.t. p
        // is sum_i dcoeff_i/dp * hf_i + coeff_i * dhf_i / dp.
        //
        // For the i-th component squared uncertainties, we have
        // d unc2_i / dp = 2 * coeff_i * dcoeff_i / dp  * unc2_i    +    coeff_i^2 * d / dp unc2_i
        Histogram1DWithUncertainties hi(histo_dimensions[iobs].nbins, histo_dimensions[iobs].xmin, histo_dimensions[iobs].xmax);
        coeff_ders.clear();
        double coeff = coeffs[i].eval_with_derivative(parameters, coeff_ders);
        
        // coeff_i * dhf_i / dp part of derivative, and the "coeff_i^2 * d / dp unc2_i" part for the derivatives:
        hfs[i].eval_and_add_derivatives(hi, der, coeff, parameters);
        // update result itself:
        h.add_with_coeff(coeff, hi);
        
        // add the dcoeff_i/dp * hf_i part of the derivative values and the   "2 * coeff_i * dcoeff_i / dp  * unc2_i" for the squared uncertainties:
        const ParIds & coeff_pids = coeffs[i].get_parameters();
        for(ParIds::const_iterator pit=coeff_pids.begin(); pit!=coeff_pids.end(); ++pit){
            der[*pit].add_with_coeff_values(coeff_ders.get_unchecked(*pit), hi);
            der[*pit].add_with_coeff_unc2(2.0 * coeff_ders.get_unchecked(*pit) * coeff, hi);
        }
    }
}

std::auto_ptr<NLLikelihood> default_model::get_nllikelihood(const Data & data) const{
    if(not(data.get_observables() == observables)){
        throw invalid_argument("default_model::get_nllikelihood: observables of model and data mismatch!");
    }
    if(not(data.get_rvobs_values().contains_all(rvobservables))){
        throw invalid_argument("default_model::get_nllikelihood: real-values observables of model and data mismatch!");
    }
    default_model_nll * result;
    if(bb_uncertainties){
        result = new default_model_bbadd_nll(*this, data);
    }
    else{
        result = new default_model_nll(*this, data);
    }
    result->robust_nll = robust_nll;
    return std::auto_ptr<NLLikelihood>(result);
}

default_model::default_model(const Configuration & ctx): bb_uncertainties(false), robust_nll(false){
    Setting s = ctx.setting;
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    if(s.exists("bb_uncertainties")){
        bb_uncertainties =  s["bb_uncertainties"];
    }
    if(s.exists("robust_nll")){
        robust_nll =  s["robust_nll"];
    }
    //go through observables to find the template definition for each of them:
    ObsIds observables = vm->get_all_observables();
    for (ObsIds::const_iterator obsit = observables.begin(); obsit != observables.end(); obsit++) {
        string obs_name = vm->get_name(*obsit);
        if(not s.exists(obs_name)) continue;
        Setting obs_setting = s[obs_name];
        if(obs_setting.size()==0) throw ConfigurationException("observable '" + vm->get_name(*obsit) + "' is empty");
        boost::ptr_vector<HistogramFunction> histos;
        boost::ptr_vector<Function> coeffs;
        for (size_t i = 0; i < obs_setting.size(); i++) {
            auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(ctx, obs_setting[i]["histogram"]));
            auto_ptr<Function> coeff_function = PluginManager<Function>::build(Configuration(ctx, obs_setting[i]["coefficient-function"]));
            coeffs.push_back(coeff_function);
            histos.push_back(hf);
        }
        set_prediction(*obsit, coeffs, histos);
    }
    if(ctx.setting.exists("rvobs-distribution")){
        rvobservable_distribution = PluginManager<Distribution>::build(Configuration(ctx, ctx.setting["rvobs-distribution"]));
        rvobservables = rvobservable_distribution->get_parameters();
        // add parameters:
        parameters.insert_all(rvobservable_distribution->get_distribution_parameters());
    }
    // type checking for rvobs ParIds vs. parameter ParIds:
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it){
        if(vm->get_type(*it) != "par"){
            throw ConfigurationException("Type error: parameter '" + vm->get_name(*it) + "' is used as model parameter, but was not declared as such.");
        }
    }
    for(ParIds::const_iterator it=rvobservables.begin(); it!=rvobservables.end(); ++it){
        if(vm->get_type(*it) != "rvobs"){
            throw ConfigurationException("Type error: parameter '" + vm->get_name(*it) + "' is used as real-valued observable, but was not declared as such.");
        }
    }
    
    // additional_nll_term:
    if(ctx.setting.exists("additional_nll_term")){
        additional_nll_term = PluginManager<Function>::build(Configuration(ctx, ctx.setting["additional_nll_term"]));
        parameters.insert_all(additional_nll_term->get_parameters());
    }
    
    // parameter distribution:
    if(parameters.size() == 0){
        parameter_distribution.reset(new EmptyDistribution());
    }
    parameter_distribution = PluginManager<Distribution>::build(Configuration(ctx, s["parameter-distribution"]));
    if(not (parameter_distribution->get_parameters() == parameters)){
        stringstream ss;
        ss << "'parameter-distribution' has to define the same set of parameters the model depends on. However";
        ParIds dist_pars = parameter_distribution->get_parameters();
        ParIds all_pars = parameters;
        all_pars.insert_all(dist_pars);
        for(ParIds::const_iterator p_it=all_pars.begin(); p_it!=all_pars.end(); ++p_it){
            if(parameters.contains(*p_it) && dist_pars.contains(*p_it)) continue;
            if(parameters.contains(*p_it)){
               ss << ", the model depends on '"<< vm->get_name(*p_it) << "' which the parameter distribution does not include";
            }
            else ss << ", the parameter distribution depends on '" << vm->get_name(*p_it) << "' which the model does not depend on";
        }
        throw ConfigurationException(ss.str());
    }
}

default_model::~default_model(){
}

/* default_model_nll */
default_model_nll::default_model_nll(const default_model & m, const Data & dat): model(m),  data(dat) {
    par_ids = model.get_parameters();
}

void default_model_nll::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
}


double default_model_nll::operator()(const ParValues & values) const{
    atomic_inc(&n_nll_eval);
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->eval_nl(values);
    }
    else{
        result += model.get_parameter_distribution().eval_nl(values);
    }
    if(std::isinf(result)) return result;
    //2. get the prediction of the model:
    model.get_prediction(predictions, values);
    //3. the template likelihood
    const std::vector<std::pair<ObsId, size_t> > & li = model.get_last_indices();
    if(robust_nll){
        for(std::vector<std::pair<ObsId, size_t> >::const_iterator li_it=li.begin(); li_it!=li.end(); ++li_it){
            const ObsId & oid = li_it->first;
            const double * pred_data = predictions[oid].get_data();
            const double * data_data = data[oid].get_data();
            result += template_nllikelihood_robust(data_data, pred_data, data[oid].get_nbins());
        }
    }
    else{
        for(std::vector<std::pair<ObsId, size_t> >::const_iterator li_it=li.begin(); li_it!=li.end(); ++li_it){
            const ObsId & oid = li_it->first;
            const double * pred_data = predictions[oid].get_data();
            const double * data_data = data[oid].get_data();
            result += template_nllikelihood(data_data, pred_data, data[oid].get_nbins());
        }
    }
    if(std::isinf(result)) return result;
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.get_rvobs_values());
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    const Function * additional_term = model.get_additional_nll_term();
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}


double default_model_nll::eval_with_derivative(const ParValues & values, ParValues & der) const{
    atomic_inc(&n_nll_eval_with_derivative);
    double result = 0.0;
    der.set_zero(model.get_parameters());
    //1. the model prior first
    if(override_distribution){
        result = override_distribution->eval_nl_with_derivative(values, der);
    }
    else{
        result = model.get_parameter_distribution().eval_nl_with_derivative(values, der);
    }
    //2., 3. get the [prediction and add the template likelihood
    const std::vector<std::pair<ObsId, size_t> > & li = model.get_last_indices();
    for(std::vector<std::pair<ObsId, size_t> >::const_iterator li_it=li.begin(); li_it!=li.end(); ++li_it){
        const ObsId & oid = li_it->first;
        Histogram1D & hpred = predictions[oid];
        model.get_prediction_with_derivative(oid, hpred, der_pred, values);
        // The negative Poisson log is
        //  -ln L = p - d * log(p)
        // where d is the data and p the predicted yield.
        // Taking the derivative w.r.t. the model parameter theta is
        // d / dtheta -ln L = dp / dtheta * (1 - d / p)
        const double * pred_data = hpred.get_data();
        const double * data_data = data[oid].get_data();
        const size_t nbins = hpred.get_nbins();
        for(size_t i = 0; i<nbins; ++i){
            const double d = data_data[i];
            const double pred = pred_data[i];
            result += pred;
            if(d > 0.0){
                if(pred <= 0.0){
                    result = numeric_limits<double>::infinity();
                }
                result -= d * utils::log(pred);
            }
        }
        for(map<ParId, Histogram1D>::const_iterator der_it=der_pred.begin(); der_it!=der_pred.end(); ++der_it){
            const double * dpred_dpid = der_it->second.get_data();
            for(size_t i = 0; i<nbins; ++i){
                const double d = data_data[i];
                const double pred = pred_data[i];
                if(pred > 0.0){
                    der.add_unchecked(der_it->first, dpred_dpid[i] * (1 - d / pred));
                }
                // otherwise: if d==0.0 we have nothing to add; if d >0 the result is infinity and the derivative does not
                // really make sense anyway ...
            }
        }
    }
    
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.get_rvobs_values());
        ParValues der_tmp(model.get_parameters());
        result += rvobs_dist->eval_nl_with_derivative(all_values, der_tmp);
        der.add(der_tmp);
    }
    //5. The additional likelihood terms, if set:
    const Function * additional_term = model.get_additional_nll_term();
    if(additional_term){
       ParValues der_tmp(model.get_parameters());
       result += additional_term->eval_with_derivative(values, der_tmp);
       der.add(der_tmp);
    }
    return result;    
}



// bbadd
default_model_bbadd_nll::default_model_bbadd_nll(const default_model & m, const Data & dat):
     default_model_nll(m, dat){
}

double default_model_bbadd_nll::operator()(const ParValues & values) const{
    atomic_inc(&n_nll_eval);
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->eval_nl(values);
    }
    else{
        result += model.get_parameter_distribution().eval_nl(values);
    }
    //2. get the prediction of the model, with uncertainties:
    model.get_prediction(predictions_wu, values);
    //3. the template likelihood. This is the only thing different w.r.t. the "non-bb" version ...
    const ObsIds & obs_ids = model.get_observables();
    for(ObsIds::const_iterator obsit=obs_ids.begin(); obsit!=obs_ids.end(); obsit++){
        const Histogram1DWithUncertainties & pred_obs = predictions_wu[*obsit];
        const Histogram1D & data_obs = data[*obsit];
        const size_t nbins = data_obs.get_nbins();
        theta_assert(nbins == pred_obs.get_nbins());
        for(size_t ibin=0; ibin < nbins; ++ibin){
            const double p = pred_obs.get_value(ibin);
            const double d = data_obs.get(ibin);
            const double p_unc2 = pred_obs.get_uncertainty2(ibin);
            double beta = 0.0;
            if(p_unc2 > 0.0){
                double dummy;
                theta::utils::roots_quad(dummy, beta, p + p_unc2, p_unc2 * (p - d));
                result += 0.5 * beta * beta / p_unc2;
            }
            const double new_pred = beta + p;
            // As special case, new_pred == 0.0 can happen (for p == p_unc2 and data == 0). In this case,
            // the log-term can be skipped as it has vanishing contribution to the nll.
            result += new_pred;
            if(d > 0.0){
                if(new_pred <= 0.0){
                    return numeric_limits<double>::infinity();
                }
                result -= d * utils::log(new_pred);
            }
        }
        
    }
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.get_rvobs_values());
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    const Function * additional_term = model.get_additional_nll_term();
    if(additional_term){
        result += (*additional_term)(values);
    }
    return result;
}


double default_model_bbadd_nll::eval_with_derivative(const ParValues & values, ParValues & der) const{
    atomic_inc(&n_nll_eval_with_derivative);
    double result = 0.0;
    der.set_zero(model.get_parameters());
    //1. the model prior first
    if(override_distribution){
        result = override_distribution->eval_nl_with_derivative(values, der);
    }
    else{
        result = model.get_parameter_distribution().eval_nl_with_derivative(values, der);
    }
    //2., 3. get the prediction and add the template likelihood
    model.get_prediction(predictions_wu, values);
    const std::vector<std::pair<ObsId, size_t> > & li = model.get_last_indices();
    size_t iobs = 0;
    for(std::vector<std::pair<ObsId, size_t> >::const_iterator li_it=li.begin(); li_it!=li.end(); ++li_it, ++iobs){
        const ObsId & oid = li_it->first;
        Histogram1DWithUncertainties & hpred_wu = predictions_wu[oid];
        model.get_prediction_with_derivative(iobs, hpred_wu, der_pred_wu, values);
        // The negative Poisson log is
        //  -ln L = p - d * log(p)
        // where d is the data and p the predicted yield.
        // Taking the derivative w.r.t. the model parameter theta is
        // d / dtheta -ln L = dp / dtheta * (1 - d / p)
        size_t ider = 0; // ider: count parameters: update the result only for ider==0; update derivatives in "der" always
        const ParIds & obs_parameters = model.get_parameters(iobs);
        map<ParId, Histogram1DWithUncertainties>::const_iterator der_it=der_pred_wu.begin();
        for(ParIds::const_iterator pit=obs_parameters.begin(); pit!=obs_parameters.end(); ++pit, ++ider){
            // obs_parameters only contains a subset of all the model parameters, and der_pred_wu contains *all* model parameters
            // So we can skip through der_pred_wu we find pit:
            while(der_it->first < *pit) ++der_it;
            //theta_assert(der_it != der_pred_wu.end());
            
            const double * data_data = data[oid].get_data();
            const double * dpred_dpid = der_it->second.get_data();
            const size_t nbins = hpred_wu.get_nbins();
            for(size_t ibin = 0; ibin<nbins; ++ibin){
                const double p = hpred_wu.get(ibin);
                const double d = data_data[ibin];
                const double p_unc2 = hpred_wu.get_uncertainty2(ibin);
                // sqrt-expression appearing in the solution of the quadratic eq. for beta. We
                // need this a couple of times, so save it for later:
                const double S = sqrt((p - p_unc2)*(p - p_unc2) + 4 * d * p_unc2); 
                double beta = 0.0;
                if(p_unc2 > 0.0){
                    beta = 0.5 * (-p -p_unc2 + S);
                    if(ider == 0){
                        result += 0.5 * beta * beta / p_unc2;
                    }
                }
                const double new_pred = beta + p;
                // As special case, new_pred == 0.0 can happen (for p == p_unc2 and data == 0). In this case,
                // the log-term can be skipped as it has vanishing contribution to the nll.
                if(ider == 0){
                    result += new_pred;
                    if(d > 0.0){
                        if(new_pred <= 0.0){
                            return numeric_limits<double>::infinity();
                        }
                        result -= d * utils::log(new_pred);
                    }
                }
                if(new_pred > 0.0){
                    // calculate the derivative of new_pred w.r.t. pid:
                    double dnewpred_dpid = 0.5 * dpred_dpid[ibin] * (1 + (p - p_unc2) / S);
                    der.add_unchecked(der_it->first, dnewpred_dpid * (1 - d / new_pred));
                    if(p_unc2 > 0.0){
                        double dbeta_dpid = dnewpred_dpid - dpred_dpid[ibin];
                        der.add_unchecked(der_it->first, beta / p_unc2 * dbeta_dpid - 0.5 * beta * beta / (p_unc2 * p_unc2) * der_it->second.get_uncertainty2(ibin));
                    }
                }
            }
        }
    }
    
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.get_rvobs_values());
        ParValues der_tmp(model.get_parameters());
        result += rvobs_dist->eval_nl_with_derivative(all_values, der_tmp);
        der.add(der_tmp);
    }
    //5. The additional likelihood terms, if set:
    const Function * additional_term = model.get_additional_nll_term();
    if(additional_term){
       ParValues der_tmp(model.get_parameters());
       result += additional_term->eval_with_derivative(values, der_tmp);
       der.add(der_tmp);
    }
    return result;    
}

REGISTER_PLUGIN_DEFAULT(default_model)
