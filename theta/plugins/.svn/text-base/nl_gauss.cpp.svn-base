#include "plugins/nl_gauss.hpp"

#include "interface/exception.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

nl_gauss::nl_gauss(const theta::Configuration & cfg) {
    const size_t n = cfg.setting["rows"].size();
    if(n==0) throw ConfigurationException("empty 'rows'");
    if(n != cfg.setting["mu"].size() || n != cfg.setting["covariance"].size()) throw ConfigurationException("inconsistent length of 'rows'/'mu'/'covariance'");
    rows.reserve(n);
    mu.resize(n);
    row_values.resize(n);
    Matrix cov_tmp(n,n);
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    for(size_t i=0; i<n; ++i){
        if(cfg.setting["rows"][i].get_type() == Setting::TypeString){
            ParId pid = vm->get_par_id(cfg.setting["rows"][i]);
            rows.push_back(new IdFunction(pid));
            par_ids.insert(pid);
        }
        else{
            rows.push_back(PluginManager<Function>::build(Configuration(cfg, cfg.setting["rows"][i])));
            par_ids.insert_all(rows.back().get_parameters());
        }
        mu[i] = cfg.setting["mu"][i];
        Setting s_row = cfg.setting["covariance"][i];
        if(s_row.size() != n){
            throw ConfigurationException("invalid cov matrix size");
        }
        for(size_t j=0; j<n; ++j){
             cov_tmp(i,j) = s_row[j];
        }
    }
    inv_cov = cov_tmp;
    inv_cov.invert_cholesky();
    theta_assert(mu.size() == n and inv_cov.get_n_rows() == n and inv_cov.get_n_cols()==n and rows.size()==n and row_values.size()==n);
}

double nl_gauss::operator()(const theta::ParValues & values) const{
    const size_t n = rows.size();
    for(size_t i=0; i<n; ++i){
        row_values[i] = rows[i](values);
    }
    double result = 0.0;
    for(size_t i=0; i<n; ++i){
        for(size_t j=0; j<n; ++j){
            result += (row_values[i] - mu[i]) * (row_values[j] - mu[j]) * inv_cov(i,j);
        }
    }
    return result / 2;
}

REGISTER_PLUGIN(nl_gauss)
