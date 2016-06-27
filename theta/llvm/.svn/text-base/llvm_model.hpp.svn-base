#ifndef LLVM_MODEL_HPP
#define LLVM_MODEL_HPP

#include "interface/decls.hpp"
#include "interface/model.hpp"
#include "llvm/llvm_interface.hpp"

#include <vector>
#include <string>
#include <set>
#include <map>


typedef double (*t_model_evaluate)(const double *, double *);
typedef double (*t_model_evaluate_unc)(const double *, double *, double*);

class llvm_model_nll;


// like default_model, with same configuration.
// One additional parameter: llvm_always = true;  (default is false) will use the
// compiled prediction function not only for likelihood construction but also for get_prediction(Data).
class llvm_model: public theta::Model{
friend class llvm_model_nll;
private:
    typedef boost::ptr_map<theta::ObsId, boost::ptr_vector<theta::HistogramFunction> > histos_type;
    typedef boost::ptr_map<theta::ObsId, boost::ptr_vector<theta::Function> > coeffs_type;
    theta::VarIdManager vm;
    histos_type histos;
    coeffs_type coeffs;
    std::auto_ptr<theta::Distribution> parameter_distribution;
    bool llvm_always;
    std::auto_ptr<theta::Distribution> rvobservable_distribution;
    std::auto_ptr<theta::Function> additional_nll_term;
    
    mutable t_model_evaluate model_evaluate;
    mutable t_model_evaluate_unc model_evaluate_unc;
    mutable std::auto_ptr<llvm_module> module;
    
    size_t nbins_total;
    
    bool bb_uncertainties;

    void set_prediction(const theta::ObsId & obs_id, boost::ptr_vector<theta::Function> & coeffs, boost::ptr_vector<theta::HistogramFunction> & histos);
    // generate and compile llvm, fill mode_get_prediction function pointer.
    void generate_llvm(bool with_uncertainties) const;

    template<typename HT>
    void get_prediction_impl(theta::DataT<HT> & result, const theta::ParValues & parameters) const;
    
 public:
    llvm_model(const theta::Configuration & cfg, bool llvm_always = false);
    //the pure virtual functions:
    virtual void get_prediction(theta::Data & result, const theta::ParValues & parameters) const;
    virtual void get_prediction(theta::DataWithUncertainties & result, const theta::ParValues & parameters) const;
    virtual std::auto_ptr<theta::NLLikelihood> get_nllikelihood(const theta::Data & data) const;
    


    virtual const theta::Distribution & get_parameter_distribution() const {
       return *parameter_distribution;
    }
    virtual ~llvm_model();
    
    virtual const theta::Function * get_additional_nll_term() const{
        return additional_nll_term.get();
    }

    virtual const theta::Distribution * get_rvobservable_distribution() const{
          return rvobservable_distribution.get();
    }
    
};
 
class llvm_model_nll: public theta::NLLikelihood{
friend class llvm_model;
public:
    using theta::Function::operator();
    virtual double operator()(const theta::ParValues & values) const;
    
    virtual void set_override_distribution(const boost::shared_ptr<theta::Distribution> & d);
    virtual const theta::Distribution & get_parameter_distribution() const{
        if(override_distribution) return *override_distribution;
        else return model.get_parameter_distribution();
    }
private:
    const llvm_model & model;
    
    void fill_par_ids_vec();
    
    std::vector<theta::ParId> par_ids_vec;
    
    theta::ParValues rvobs_values;
    boost::shared_ptr<theta::Distribution> override_distribution;

    theta::Histogram1D data_concatenated;
    t_model_evaluate model_evaluate;
    //cached predictions:
    mutable theta::Histogram1D pred_concatenated;
    mutable std::vector<double> parameter_values;
    
    llvm_model_nll(const llvm_model & m, const theta::Data & data, t_model_evaluate model_evaluate, size_t nbins_total);
};

class llvm_model_nll_bb: public theta::NLLikelihood{
friend class llvm_model;
public:
    using theta::Function::operator();
    virtual double operator()(const theta::ParValues & values) const;

    virtual void set_override_distribution(const boost::shared_ptr<theta::Distribution> & d);
    virtual const theta::Distribution & get_parameter_distribution() const{
        if(override_distribution) return *override_distribution;
        else return model.get_parameter_distribution();
    }
private:
    const llvm_model & model;

    void fill_par_ids_vec();

    std::vector<theta::ParId> par_ids_vec;

    theta::ParValues rvobs_values;
    boost::shared_ptr<theta::Distribution> override_distribution;

    theta::Histogram1D data_concatenated;
    t_model_evaluate_unc model_evaluate_unc;

    //cached predictions:
    mutable theta::Histogram1D pred_concatenated;
    mutable theta::Histogram1D pred_unc2_concatenated;
    mutable std::vector<double> parameter_values;

    llvm_model_nll_bb(const llvm_model & m, const theta::Data & data, t_model_evaluate_unc model_evaluate_unc, size_t nbins_total);
};

#endif


