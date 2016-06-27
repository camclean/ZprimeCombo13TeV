#include "plugins/add.hpp"
#include "interface/plugin.hpp"

using namespace theta;
using namespace std;

add::add(const Configuration & cfg): literal_addend(0.0){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    size_t n = cfg.setting["addends"].size();
    for(size_t i=0; i<n; ++i){
        Setting::Type t = cfg.setting["addends"][i].get_type();
        if(t==Setting::TypeFloat){
            literal_addend += static_cast<double>(cfg.setting["addends"][i]);
        }
        else if(t==Setting::TypeString){
           ParId pid = vm->get_par_id(cfg.setting["addends"][i]);
           v_pids.push_back(pid);
           par_ids.insert(pid);
        }
        else if(t==Setting::TypeGroup){
            std::auto_ptr<Function> f = PluginManager<Function>::build(Configuration(cfg, cfg.setting["addends"][i]));
            par_ids.insert_all(f->get_parameters());
            functions.push_back(f);
        }
        else{
           std::stringstream ss;
           ss << "Invalid config setting type in setting 'addends' at index " << i;
           throw ConfigurationException(ss.str());
        }
    }
}

double add::operator()(const ParValues & v) const{
    double result = literal_addend;
    for(size_t i=0; i<v_pids.size(); ++i){
        result += v.get_unchecked(v_pids[i]);
    }
    for(size_t i=0; i<functions.size(); ++i){
        result += functions[i](v);
    }
    return result;
}


REGISTER_PLUGIN(add)

