#include "llvm/llvm_model.hpp"
#include "interface/log2_dot.hpp"
#include "interface/plugin.hpp"
#include "interface/utils.hpp"
#include "interface/distribution.hpp"

#include <iostream>

#include "llvm/Module.h"
#include "llvm/DerivedTypes.h"
#include "llvm/LLVMContext.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"

#include "llvm/Support/IRBuilder.h"

using namespace std;
using namespace theta;
using namespace llvm;
using namespace theta::utils;

namespace {
    // conversion utilities from theta to llvm data structures: the first three are for
    // converting data / histograms
    size_t get_total_nbins(const ObsIds & oids, const VarIdManager & vm){
        size_t total_nbins = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            total_nbins += vm.get_nbins(*it);
            if(total_nbins % 2) ++total_nbins;
        }
        return total_nbins;
    }
    
    // make the concatenated histogram, based on the observables in oids and binning info in vm.
    // assumes that h_concat has been allocated correctly with get_total_nbins(!!)
    void get_concatenated_from_data(Histogram1D & h_concat, const Data & dat, const ObsIds & oids){
        size_t offset = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            size_t nbins = dat[*it].get_nbins();
            utils::copy_fast(h_concat.get_data() + offset, dat[*it].get_data(), nbins);
            if(nbins % 2){
                theta_assert(*(h_concat.get_data() + offset + nbins) == 0);
            }
            offset += nbins;
            if(offset % 2) ++offset;
        }
    }
    
    void get_data_from_concatenated(Data & d, const Histogram1D & h_concat, const ObsIds & oids, const VarIdManager & vm){
        theta_assert(get_total_nbins(oids, vm) == h_concat.size());
        size_t offset = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            size_t nbins = vm.get_nbins(*it);
            const std::pair<double, double> & range = vm.get_range(*it);
            d[*it] = Histogram1D(nbins, range.first, range.second);
            utils::copy_fast(d[*it].get_data(), h_concat.get_data() + offset, nbins);
            if(nbins % 2){
                theta_assert(*(h_concat.get_data() + offset + nbins) == 0);
            }
            offset += nbins;
            if(offset % 2) ++offset;
        }
    }
    
    void get_data_from_concatenated(DataWithUncertainties & d, const Histogram1D & h_concat, const Histogram1D & h_unc2_concat, const ObsIds & oids, const VarIdManager & vm){
        theta_assert(get_total_nbins(oids, vm) == h_concat.size());
        size_t offset = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            size_t nbins = vm.get_nbins(*it);
            const std::pair<double, double> & range = vm.get_range(*it);
            Histogram1DWithUncertainties & h = d[*it] = Histogram1DWithUncertainties(nbins, range.first, range.second);
            const double * pred_data = h_concat.get_data() + offset;
            const double * pred_unc2_data = h_unc2_concat.get_data() + offset;
            for(size_t i=0; i< nbins; ++i){
            	h.set_unc2(i, pred_data[i], pred_unc2_data[i]);
            }
            if(nbins % 2){
                theta_assert(*(h_concat.get_data() + offset + nbins) == 0);
            }
            offset += nbins;
            if(offset % 2) ++offset;
        }
    }

    // convert ParValues into double * ( use &par_values[0])
    void get_par_values(vector<double> & par_values, const ParValues & values, const ParIds & par_ids){
        par_values.resize(par_ids.size());
        size_t i=0;
        for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
            par_values[i] = values.get_unchecked(*it);
        }
    }
    
    
    // par_values must have correct size!
    void get_par_values(vector<double> & par_values, const ParValues & values, const vector<ParId> & par_ids){
        //par_values.resize(par_ids.size());
        size_t i=0;
        for(vector<ParId>::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
            par_values[i] = values.get_unchecked(*it);
        }
    }
}


class LLVMPluginBuilderF: public PluginBuilder<theta::Function>{
public:
    virtual std::auto_ptr<theta::Function> build(const Configuration & cfg, const std::string & type){
        string new_type = type;
        if(type=="exp_function" || type == "multiply") new_type = "llvm_" + type;
        return PluginManager<theta::Function>::build_type(cfg, new_type);
    }
};

class LLVMPluginBuilderHF: public PluginBuilder<theta::HistogramFunction>{
public:
    virtual std::auto_ptr<theta::HistogramFunction> build(const Configuration & cfg, const std::string & type){
        string new_type = type;
        if(type == "cubiclinear_histomorph") new_type = "llvm_" + type;
        return PluginManager<theta::HistogramFunction>::build_type(cfg, new_type);
    }
};



struct llvm_pb_sentinel{
   llvm_pb_sentinel(){
       std::auto_ptr<PluginBuilder<theta::Function> > pbf(new LLVMPluginBuilderF());
       PluginManager<theta::Function>::set_plugin_builder(pbf);
       std::auto_ptr<PluginBuilder<theta::HistogramFunction> > pbhf(new LLVMPluginBuilderHF());
       PluginManager<theta::HistogramFunction>::set_plugin_builder(pbhf);

   }
   ~llvm_pb_sentinel(){
       PluginManager<theta::Function>::reset_plugin_builder();
       PluginManager<theta::HistogramFunction>::reset_plugin_builder();
   }
};


llvm_model::~llvm_model(){}

void llvm_model::set_prediction(const ObsId & obs_id, boost::ptr_vector<theta::Function> & coeffs_, boost::ptr_vector<HistogramFunction> & histos_){
    observables.insert(obs_id);
    const size_t n = coeffs_.size();
    if(n!=coeffs_.size()) throw invalid_argument("Model::setPrediction: number of histograms and coefficients do not match");
    if(histos[obs_id].size()>0 || coeffs[obs_id].size()>0)
        throw invalid_argument("Model::setPrediction: prediction already set for this observable");
    coeffs[obs_id].transfer(coeffs[obs_id].end(), coeffs_.begin(), coeffs_.end(), coeffs_);
    histos[obs_id].transfer(histos[obs_id].end(), histos_.begin(), histos_.end(), histos_);
    for(boost::ptr_vector<theta::Function>::const_iterator it=coeffs[obs_id].begin(); it!=coeffs[obs_id].end(); ++it){
        ParIds pids = (*it).get_parameters();
        parameters.insert_all(pids);
    }
    size_t nbins = 0;
    double xmin = NAN, xmax = NAN;
    bool first = true;
    for(boost::ptr_vector<HistogramFunction>::const_iterator it=histos[obs_id].begin(); it!=histos[obs_id].end(); ++it){
        if(first){
            it->get_histogram_dimensions(nbins, xmin, xmax);
            first = false;
        }
        else{
            size_t nbins_tmp = 0;
            double xmin_tmp = NAN, xmax_tmp = NAN;
            it->get_histogram_dimensions(nbins_tmp, xmin_tmp, xmax_tmp);
            if(nbins!=nbins_tmp || xmin!=xmin_tmp || xmax!=xmax_tmp){
                throw invalid_argument("llvm_model::set_prediction: histogram dimensions mismatch");
            }
        }
        const ParIds & pids = (*it).get_parameters();
        parameters.insert_all(pids);
    }
}


template<typename HT>
void llvm_model::get_prediction_impl(DataT<HT> & result, const ParValues & parameters) const{
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    for(; h_it != histos.end(); ++h_it, ++c_it){
        const ObsId & oid = h_it->first;
        histos_type::const_mapped_reference hfs = *(h_it->second);
        coeffs_type::const_mapped_reference h_coeffs = *(c_it->second);
        //theta_assert(hfs.size() > 0 && hfs.size() == h_coeffs.size());
        // overwrite result[oid] with first term:
        hfs[0].apply_functor(copy_to<HT>(result[oid]), parameters);
        result[oid] *= h_coeffs[0](parameters);
        // add the rest:
        for (size_t i = 1; i < hfs.size(); i++) {
            hfs[i].apply_functor(add_with_coeff_to<HT>(result[oid], h_coeffs[i](parameters)), parameters);
        }
    }
}

void llvm_model::get_prediction(Data & result, const ParValues & parameters) const {
    if(llvm_always){
        if(model_evaluate==0) generate_llvm(false);
        theta_assert(model_evaluate);
        const size_t n_total_bins = get_total_nbins(observables, vm);
        Histogram1D pred_concat(n_total_bins);
        vector<double> parameter_values(this->parameters.size());
        get_par_values(parameter_values, parameters, this->parameters);
      	model_evaluate(&parameter_values[0], pred_concat.get_data());
        get_data_from_concatenated(result, pred_concat, observables, vm);
    }
    else{
    	get_prediction_impl<Histogram1D>(result, parameters);
    }
}

void llvm_model::get_prediction(DataWithUncertainties & result, const theta::ParValues & parameters) const{
	if(llvm_always){
		if(model_evaluate_unc==0) generate_llvm(true);
		theta_assert(model_evaluate_unc);
		const size_t n_total_bins = get_total_nbins(observables, vm);
		Histogram1D pred_concat(n_total_bins);
		vector<double> parameter_values(this->parameters.size());
		get_par_values(parameter_values, parameters, this->parameters);
		Histogram1D pred_unc2_concat(n_total_bins);
		model_evaluate_unc(&parameter_values[0], pred_concat.get_data(), pred_unc2_concat.get_data());
		get_data_from_concatenated(result, pred_concat, pred_unc2_concat, observables, vm);
	}
	else{
		get_prediction_impl<Histogram1DWithUncertainties>(result, parameters);
	}
}

namespace{
Histogram1D get_uncertainty2_histo(const theta::HistogramFunction & f_){
	const llvm_enabled_histogram_function * f = dynamic_cast<const llvm_enabled_histogram_function *>(&f_);
	if(f!=0){
		return f->get_uncertainty2_histogram();
	}
	// handle constant histograms:
	if(f_.get_parameters().size()==0){
		Histogram1DWithUncertainties h_wu;
		f_.apply_functor(copy_to<Histogram1DWithUncertainties>(h_wu), ParValues());
		return h_wu.get_uncertainty2_histogram();
	}
	else{
		string hf_name = demangle(typeid(f_).name());
		throw invalid_argument("HistogramFunction '"+ hf_name +"' does not support llvm uncertainty.");
	}
}

/*void dump_module(llvm::Module & module, const std::string & fname){
	PassManager pm;
	std::string err;
	raw_ostream * os = new raw_fd_ostream(fname.c_str(), err);
	if(!err.empty()) throw invalid_argument("could not open module file " + fname);
	pm.add(createPrintModulePass(os, true));
	pm.run(module);
}*/
}

// the generated function is
//    void model_evaluate(const double * par_values, double * data);  -- for with_uncertainties = false
//  or
//    void model_evaluate_unc(const double * par_values, double * data, double * data_unc2);  -- for with_uncertainties = true
// where data and data_unc2 must both be set to zero prior to calling this function.
void llvm_model::generate_llvm(bool with_uncertainties) const {
    if(!module.get()) module.reset(new llvm_module(parameters));
    LLVMContext & context = module->module->getContext();
    Type * double_t = Type::getDoubleTy(context);
    Type * i32_t = Type::getInt32Ty(context);
    Type * void_t = Type::getVoidTy(context);
    std::vector<Type*> arg_types(3, double_t->getPointerTo());
    std::string suffix;
    if(with_uncertainties){
    	suffix = "_unc";
    }
    else{
    	arg_types.resize(2);
    }
    FunctionType * FT = FunctionType::get(void_t, arg_types, false);
    // note: has to be linked externally to prevent the optimizer from optimizing it away completely ...
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "model_evaluate" + suffix, module->module);
    llvm::Function::arg_iterator iter = F->arg_begin();
    // the function parameters:
    Value * par_values = iter++;
    Value * data = iter++;
    Value * data_unc2 = 0;
    if(with_uncertainties) data_unc2 = iter;
    IRBuilder<> Builder(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    // generate code for all Functions and HistorgamFunctions:
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    size_t iobs = 0;
    size_t obs_offset = 0;
    for(; h_it != histos.end(); ++h_it, ++c_it, ++iobs){
        theta_assert(h_it->first == c_it->first);
        histos_type::const_mapped_reference h_functions = *(h_it->second);
        coeffs_type::const_mapped_reference h_coeffs = *(c_it->second);
        theta_assert(h_functions.size() == h_coeffs.size());
        // data_iobs = data + obs_offset:
        Value * data_iobs = Builder.CreateGEP(data, ConstantInt::get(i32_t, obs_offset), "data_iobs");
        Value * data_unc2_iobs = 0;
        if (with_uncertainties) data_unc2_iobs = Builder.CreateGEP(data_unc2, ConstantInt::get(i32_t, obs_offset), "data_unc2_iobs");
        for(size_t i=0; i<h_functions.size(); ++i){
            stringstream ss_prefix;
            ss_prefix << "c" << iobs << "_p" << i;
            string prefix = ss_prefix.str();
            llvm::Function * coeff_function = create_llvm_function(&h_coeffs[i], *module, prefix);
            Value * coeff = Builder.CreateCall(coeff_function, par_values, "coeff");
            llvm::Function * histo_function = create_llvm_histogram_function(&h_functions[i], *module, prefix);
            llvm::Value * c = Builder.CreateCall3(histo_function, coeff, par_values, data_iobs);
            if(with_uncertainties){
            	Value * coeff2 = Builder.CreateFMul(coeff, coeff, "c2");
            	Value * coeff3 = Builder.CreateFMul(coeff2, c, "c3");
            	Value * coeff4 = Builder.CreateFMul(coeff3, c, "c4");
            	BB = module->emit_add_with_coeff(BB, coeff4, data_unc2_iobs, get_uncertainty2_histo(h_functions[i]));
            	Builder.SetInsertPoint(BB);
            }
        }
        obs_offset += vm.get_nbins(h_it->first);
        if(obs_offset % 2) ++obs_offset;
    }
    Builder.CreateRetVoid();
    module->optimize();
    if(with_uncertainties){
    	model_evaluate_unc = reinterpret_cast<t_model_evaluate_unc>(module->getFunctionPointer(F));
    }
    else{
    	model_evaluate = reinterpret_cast<t_model_evaluate>(module->getFunctionPointer(F));
    }
}

std::auto_ptr<NLLikelihood> llvm_model::get_nllikelihood(const Data & data) const{
    if(not(data.get_observables()==observables)){
        throw invalid_argument("Model::createNLLikelihood: observables of model and data mismatch!");
    }
    if(bb_uncertainties){
    	if(model_evaluate_unc == 0){
			generate_llvm(true);
		}
    	return std::auto_ptr<NLLikelihood>(new llvm_model_nll_bb(*this, data, model_evaluate_unc, nbins_total));
    }
    else{
    	if(model_evaluate == 0){
			generate_llvm(false);
		}
    	return std::auto_ptr<NLLikelihood>(new llvm_model_nll(*this, data, model_evaluate, nbins_total));
    }
}


llvm_model::llvm_model(const Configuration & ctx, bool llvm_always_): vm(*ctx.pm->get<VarIdManager>()),
 llvm_always(llvm_always_), model_evaluate(0), model_evaluate_unc(0), bb_uncertainties(false){
    Setting s = ctx.setting;
    if(s.exists("bb_uncertainties")){
        bb_uncertainties =  s["bb_uncertainties"];
    }
    ObsIds observables = vm.get_all_observables();
    llvm_pb_sentinel b;
    //go through observables to find the template definition for each of them:
    for (ObsIds::const_iterator obsit = observables.begin(); obsit != observables.end(); obsit++) {
        string obs_name = vm.get_name(*obsit);
        if(not s.exists(obs_name)) continue;
        Setting obs_setting = s[obs_name];
        boost::ptr_vector<HistogramFunction> histos;
        boost::ptr_vector<theta::Function> coeffs;
        for (size_t i = 0; i < obs_setting.size(); i++) {
            Configuration context(ctx, obs_setting[i]["histogram"]);
            auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(context);
            Configuration cfg(ctx, obs_setting[i]["coefficient-function"]);
            auto_ptr<theta::Function> coeff_function = PluginManager<theta::Function>::build(cfg);
            coeffs.push_back(coeff_function);
            theta_assert(coeff_function.get()==0);
            histos.push_back(hf);
        }
        set_prediction(*obsit, coeffs, histos);
    }

    if(ctx.setting.exists("llvm_always")){
        llvm_always = ctx.setting["llvm_always"];
    }
    if(ctx.setting.exists("rvobs-distribution")){
        rvobservable_distribution = PluginManager<Distribution>::build(Configuration(ctx, ctx.setting["rvobs-distribution"]));
        rvobservables = rvobservable_distribution->get_parameters();
        parameters.insert_all(rvobservable_distribution->get_distribution_parameters());
    }
    // type checking:
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it){
        if(vm.get_type(*it) != "par"){
            throw ConfigurationException("Type error: parameter '" + vm.get_name(*it) + "' is used as model parameter, but was not declared as such.");
        }
    }
    for(ParIds::const_iterator it=rvobservables.begin(); it!=rvobservables.end(); ++it){
        if(vm.get_type(*it) != "rvobs"){
            throw ConfigurationException("Type error: parameter '" + vm.get_name(*it) + "' is used as real-valued observable, but was not declared as such.");
        }
    }
    parameter_distribution = PluginManager<Distribution>::build(Configuration(ctx, s["parameter-distribution"]));
	if(not (parameter_distribution->get_parameters() == get_parameters())){
		stringstream ss;
		ss << "'parameter-distribution' has to define the same set of parameters the model depends on. However";
		ParIds dist_pars = parameter_distribution->get_parameters();
		const ParIds & m_pars = parameters;
		ParIds all_pars = parameters;
		all_pars.insert_all(dist_pars);
		for(ParIds::const_iterator p_it=all_pars.begin(); p_it!=all_pars.end(); ++p_it){
			if(m_pars.contains(*p_it) && dist_pars.contains(*p_it)) continue;
			if(m_pars.contains(*p_it)){
			   ss << ", the model depends on '"<< vm.get_name(*p_it) << "' which the parameter distribution does not include";
			}
			else ss << ", the parameter distribution depends on '" << vm.get_name(*p_it) << "' which the model does not depend on";
		}
		throw ConfigurationException(ss.str());
	}
    nbins_total = get_total_nbins(observables, vm);
}


/* llvm_model_nll */
void llvm_model_nll::fill_par_ids_vec(){
    par_ids_vec.clear();
    for(ParIds::const_iterator it = par_ids.begin(); it!=par_ids.end(); ++it){
        par_ids_vec.push_back(*it);
    }
    parameter_values.resize(par_ids_vec.size());
}


llvm_model_nll::llvm_model_nll(const llvm_model & m, const Data & dat, t_model_evaluate model_evaluate_, size_t nbins_total): model(m),
  data_concatenated(nbins_total, 0, 1), model_evaluate(model_evaluate_){
    theta::Function::par_ids = model.get_parameters();
    fill_par_ids_vec();
    get_concatenated_from_data(data_concatenated, dat, model.get_observables());
    pred_concatenated = Histogram1D(data_concatenated.size(), 0, 1);
    rvobs_values = dat.get_rvobs_values();
}

void llvm_model_nll::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
}

double llvm_model_nll::operator()(const ParValues & values) const{
    //0. convert values to vector:
    get_par_values(parameter_values, values, par_ids_vec);
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->eval_nl(values);
    }
    else{
        result += model.get_parameter_distribution().eval_nl(values);
    }
    //2. get the prediction of the model:
    pred_concatenated.set_all_values(0.0);
    model_evaluate(&parameter_values[0], pred_concatenated.get_data());
    //3. the template likelihood
    result += template_nllikelihood(data_concatenated.get_data(), pred_concatenated.get_data(), data_concatenated.size());
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(rvobs_values);
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    const Function * additional_term = model.get_additional_nll_term();
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}

/* llvm_model_nll_bb */
void llvm_model_nll_bb::fill_par_ids_vec(){
    par_ids_vec.clear();
    for(ParIds::const_iterator it = par_ids.begin(); it!=par_ids.end(); ++it){
        par_ids_vec.push_back(*it);
    }
    parameter_values.resize(par_ids_vec.size());
}


llvm_model_nll_bb::llvm_model_nll_bb(const llvm_model & m, const Data & dat, t_model_evaluate_unc model_evaluate_unc_, size_t nbins_total): model(m),
  data_concatenated(nbins_total, 0, 1), model_evaluate_unc(model_evaluate_unc_){
    theta::Function::par_ids = model.get_parameters();
    fill_par_ids_vec();
    get_concatenated_from_data(data_concatenated, dat, model.get_observables());
    pred_concatenated = Histogram1D(data_concatenated.size(), pred_concatenated.get_xmin(), pred_concatenated.get_xmax());
    pred_unc2_concatenated = pred_concatenated;
    rvobs_values = dat.get_rvobs_values();
}


void llvm_model_nll_bb::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
}

double llvm_model_nll_bb::operator()(const ParValues & values) const{
    //0. convert values to vector:
    get_par_values(parameter_values, values, par_ids_vec);
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->eval_nl(values);
    }
    else{
        result += model.get_parameter_distribution().eval_nl(values);
    }
    //2. get the prediction of the model:
    pred_concatenated.set_all_values(0.0);
    pred_unc2_concatenated.set_all_values(0.0);
    model_evaluate_unc(&parameter_values[0], pred_concatenated.get_data(), pred_unc2_concatenated.get_data());
    //3. the template likelihood
    const size_t nbins = pred_concatenated.size();
    for(size_t ibin=0; ibin < nbins; ++ibin){
		const double p = pred_concatenated.get(ibin);
		const double d = data_concatenated.get(ibin);
		const double p_unc2 = pred_unc2_concatenated.get(ibin);
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
    //result += template_nllikelihood(data_concatenated.get_data(), pred_concatenated.get_data(), data_concatenated.size());
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(rvobs_values);
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    const Function * additional_term = model.get_additional_nll_term();
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}


REGISTER_PLUGIN(llvm_model)

