#include "interface/distribution.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.tcc"

using namespace theta;
using namespace std;

REGISTER_PLUGIN_BASETYPE(theta::Distribution);

Distribution::~Distribution(){}

double Distribution::eval_nl_with_derivative(const ParValues & values, ParValues & derivative) const{
    throw Exception("not implemented for " + demangle(typeid(*this).name()));
}

EmptyDistribution::~EmptyDistribution(){}

const std::pair<double, double> & EmptyDistribution::support(const ParId & p) const{
    throw invalid_argument("called EmptyDistribution::support");
}


Ranges::Ranges(const std::map<ParId, std::pair<double, double> > & ranges_): ranges(ranges_){
}

Ranges::Ranges(const theta::Distribution & dist){
    set_from(dist);
}

void Ranges::set_from(const theta::Distribution & dist){
    const ParIds & pars = dist.get_parameters();
    for(ParIds::const_iterator it=pars.begin(); it!=pars.end(); ++it){
        ranges[*it] = dist.support(*it);
        theta_assert(ranges[*it].first <= ranges[*it].second);
   }
}

void Ranges::set(const ParId & pid, const std::pair<double, double> & range){
    if(range.first > range.second){
        throw invalid_argument("invalid range");
    }
    ranges[pid] = range;
}

const pair<double, double> & Ranges::get(const ParId & pid) const{
    map<ParId, pair<double, double> >::const_iterator it = ranges.find(pid);
    if(it==ranges.end()){
        throw invalid_argument("no range set for this parid");
    }
    return it->second;
}

bool Ranges::fixed(const ParId & pid) const{
    map<ParId, pair<double, double> >::const_iterator it = ranges.find(pid);
    if(it==ranges.end()){
        throw invalid_argument("no range set for this parid");
    }
    return it->second.first == it->second.second;
}

// truncate the values to the ranges.
void Ranges::trunc(ParValues & values) const{
    for(map<ParId, pair<double, double> >::const_iterator it=ranges.begin(); it!=ranges.end(); ++it){
        if(!values.contains(it->first)) continue;
        double val = values.get_unchecked(it->first);
        if(val < it->second.first) val = it->second.first;
        else if(val > it->second.second) val = it->second.second;
        values.set(it->first, val);
    }
}

