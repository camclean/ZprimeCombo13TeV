#include "plugins/multiply.hpp"
#include "interface/plugin.hpp"
#include "interface/redirect_stdio.hpp"

using namespace theta;
using namespace std;

multiply::multiply(const Configuration & cfg): literal_factor(1.0){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    string type = cfg.setting["type"];
    if(type == "mult"){
        theta::out << "Warning: function plugin with type='mult' is obsolete. Use type='multiply' instead and adapt configuration accordingly (see documentation; in particular use 'factors' setting instead of 'parameters')." << endl;
        //compatibility mode: search for "parameters", instead of "factors"
        size_t n = cfg.setting["parameters"].size();
        for(size_t i=0; i<n; ++i){
            ParId pid = vm->get_par_id(cfg.setting["parameters"][i]);
            v_pids.push_back(pid);
            par_ids.insert(pid);
        }
        return;
    }
    multiple_f_per_par = false;
    size_t n = cfg.setting["factors"].size();
    for(size_t i=0; i<n; ++i){
        size_t n_pars_before = par_ids.size();
        size_t n_pars_i = 0;
        Setting::Type t = cfg.setting["factors"][i].get_type();
        if(t==Setting::TypeFloat){
            literal_factor *= static_cast<double>(cfg.setting["factors"][i]);
        }
        else if(t==Setting::TypeString){
           ParId pid = vm->get_par_id(cfg.setting["factors"][i]);
           v_pids.push_back(pid);
           par_ids.insert(pid);
           n_pars_i = 1;
        }
        else if(t==Setting::TypeGroup){
            std::auto_ptr<Function> f = PluginManager<Function>::build(Configuration(cfg, cfg.setting["factors"][i]));
            par_ids.insert_all(f->get_parameters());
            n_pars_i = f->get_parameters().size();
            functions.push_back(f);
        }
        else{
           std::stringstream ss;
           ss << "Invalid config setting type in setting 'factors' at index " << i;
           throw ConfigurationException(ss.str());
        }
        if(par_ids.size() != n_pars_before + n_pars_i){
            multiple_f_per_par = true;
        }
    }
}

double multiply::operator()(const ParValues & v) const{
    double result = literal_factor;
    std::vector<theta::ParId>::const_iterator it_end = v_pids.end();
    for(std::vector<theta::ParId>::const_iterator it = v_pids.begin(); it!= it_end; ++it){
        result *= v.get_unchecked(*it);
    }
    boost::ptr_vector<theta::Function>::const_iterator fit_end = functions.end();
    for(boost::ptr_vector<theta::Function>::const_iterator fit = functions.begin(); fit != fit_end; ++fit){
        result *= (*fit)(v);
    }
    return result;
}

double multiply::eval_with_derivative(const theta::ParValues & v, theta::ParValues & der) const{
    if(multiple_f_per_par){
        throw invalid_argument("multiply: two or more factors share a parameter. For this case, the derivative is not implemented.");
    }
    double result = literal_factor;
    std::vector<theta::ParId>::const_iterator it_end = v_pids.end();
    for(std::vector<theta::ParId>::const_iterator it = v_pids.begin(); it!= it_end; ++it){
        result *= v.get_unchecked(*it);
        der.set(*it, 1.0);
    }
    boost::ptr_vector<theta::Function>::const_iterator fit_end = functions.end();
    for(boost::ptr_vector<theta::Function>::const_iterator fit = functions.begin(); fit != fit_end; ++fit){
        result *= fit->eval_with_derivative(v, der);
    }
    return result;
}


REGISTER_PLUGIN(multiply)
REGISTER_PLUGIN_NAME(multiply,mult)
