#include "interface/variables.hpp"
#include "interface/variables-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/exception.hpp"
#include "interface/utils.hpp"

#include <sstream>
#include <cmath>

using namespace theta;
using namespace theta::utils;
using namespace std;

ParId VarIdManager::create_par_id(const std::string & name, const std::string & type) {
    if (name_to_pid.find(name) != name_to_pid.end()) {
        stringstream ss;
        ss << "VarIdManager::createParId: parameter '"<< name <<"' defined twice";
        throw invalid_argument(ss.str());
    }
    ParId result(next_pid_id);
    ++next_pid_id;
    pid_to_name[result] = name;
    pid_to_type[result] = type;
    name_to_pid.insert(make_pair(name, result));
    return result;
}

ObsId VarIdManager::create_obs_id(const std::string & name, size_t nbins, double min, double max) {
    if (name_to_oid.find(name) != name_to_oid.end()) {
        stringstream ss;
        ss << "VarIdManager::createObsId: observable '" << name << "' defined twice";
        throw invalid_argument(ss.str());
    }
    if (min >= max) {
        stringstream ss;
        ss << "Observable " << name << " has min >= max, i.e., empty range";
        throw invalid_argument(ss.str());
    }
    if(nbins==0){
        stringstream ss;
        ss << "Observable '" << name << "' has no bins";
        throw invalid_argument(ss.str());
    }
    
    ObsId result(next_oid_id);
    ++next_oid_id;
    oid_to_name[result] = name;
    name_to_oid.insert(make_pair(name, result));
    oid_to_range[result] = make_pair(min, max);
    oid_to_nbins[result] = nbins;
    return result;
}


std::string VarIdManager::get_name(const ParId & id) const {
    std::map<ParId, std::string>::const_iterator it = pid_to_name.find(id);
    if (it == pid_to_name.end()) {
        throw invalid_argument("VarIdManager::get_name: did not find given ParId.");
    }
    return it->second;
}

std::string VarIdManager::get_name(const ObsId & id) const {
    std::map<ObsId, std::string>::const_iterator it = oid_to_name.find(id);
    if (it == oid_to_name.end()) {
        throw invalid_argument("VarIdManager::get_name: did not find given ObsId.");
    }
    return it->second;
}

std::string VarIdManager::get_type(const ParId & id) const {
    std::map<ParId, std::string>::const_iterator it = pid_to_type.find(id);
    if (it == pid_to_type.end()) {
        throw invalid_argument("VarIdManager::getType: did not find given ParId.");
    }
    return it->second;
}

ParId VarIdManager::get_par_id(const std::string & name) const {
    std::map<std::string, ParId>::const_iterator it = name_to_pid.find(name);
    if (it == name_to_pid.end()) {
        stringstream ss;
        ss << __FUNCTION__ << ": did not find variable '" << name << "'";
        throw invalid_argument(ss.str());
    }
    return ParId(it->second);
}

ObsId VarIdManager::get_obs_id(const std::string & name) const {
    std::map<std::string, ObsId>::const_iterator it = name_to_oid.find(name);
    if (it == name_to_oid.end()) {
        stringstream ss;
        ss << __FUNCTION__ << ": did not find variable '" << name << "'";
        throw invalid_argument(ss.str());
    }
    return ObsId(it->second);
}

size_t VarIdManager::get_nbins(const ObsId & id) const{
    std::map<ObsId, size_t>::const_iterator it = oid_to_nbins.find(id);
    if (it == oid_to_nbins.end()) {
        throw invalid_argument("VarIdManager::get_nbins: did not find given variable id.");
    }
    return it->second;
}

const pair<double, double> & VarIdManager::get_range(const ObsId & id) const{
    std::map<ObsId, pair<double, double> >::const_iterator it = oid_to_range.find(id);
    if (it == oid_to_range.end()) {
        throw invalid_argument("VarIdManager::get_range: did not find given variable id.");
    }
    return it->second;
}

ObsIds VarIdManager::get_all_observables() const{
    std::map<ObsId, pair<double, double> >::const_iterator it = oid_to_range.begin();
    ObsIds result;
    for(; it!= oid_to_range.end(); ++it){
       result.insert(it->first);
    }
    return result;
}

ParIds VarIdManager::get_all_parameters() const{
    std::map<ParId, string>::const_iterator it = pid_to_name.begin();
    ParIds result;
    for(; it!= pid_to_name.end(); ++it){
       result.insert(it->first);
    }
    return result;
}

/* ParValues */
void ParValues::fail_get(const ParId & pid) const{
    std::stringstream ss;
    ss << "ParValues::get: given VarId " << pid.id << " not found";
    throw invalid_argument(ss.str());
}

std::ostream & theta::operator<<(std::ostream & out, const ParIds & pids){
    for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
        out << *it << " ";
    }
    return out;
}


void theta::apply_vm_settings(Configuration & ctx){
    Setting s = ctx.setting;
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    if(s.exists("observables")){
        size_t nobs = s["observables"].size();
        for (size_t i = 0; i < nobs; ++i) {
            string obs_name = s["observables"][i].get_name();
            double min = s["observables"][i]["range"][0].get_double_or_inf();
            double max = s["observables"][i]["range"][1].get_double_or_inf();
            int nbins = s["observables"][i]["nbins"];
            if(nbins <= 0){
                throw ConfigurationException("Observable '" + obs_name + "' has nbins <= 0 which is not allowed.");
            }
            vm->create_obs_id(obs_name, static_cast<size_t>(nbins), min, max);
        }
    }
    if(s.exists("parameters")){
        //get parameters:
        size_t npar = s["parameters"].size();
        for (size_t i = 0; i < npar; i++) {
            string par_name = s["parameters"][i];
            vm->create_par_id(par_name);
        }
    }
    if(s.exists("rvobservables")){
        size_t npar = s["rvobservables"].size();
        for (size_t i = 0; i < npar; i++) {
            string par_name = s["rvobservables"][i];
            vm->create_par_id(par_name, "rvobs");
        }
    }
}

