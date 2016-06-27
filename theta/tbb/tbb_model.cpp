#include "tbb/tbb_model.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"

#include <tbb/tbb.h>

using namespace theta;
using namespace tbb;

namespace{
int get_int(const Configuration & cfg, const std::string & key, int def){
    if(cfg.setting.exists(key)){
        return cfg.setting[key];
    }
    else{
        return def;
    }
}
}

template<typename HT>
void tbb_model::get_prediction_impl(DataT<HT> & result, const ParValues & parameters) const{
    // allocate all Histograms in "result" before the multi-threaded part:
    size_t nobs = last_indices.size();
    for(size_t i=0; i<nobs; ++i){
        const ObsId & oid = last_indices[i].first;
        result[oid].reset(histo_dimensions[i].nbins, histo_dimensions[i].xmin, histo_dimensions[i].xmax);
    }
    parallel_for(blocked_range<size_t>(0, last_indices.size()), [&](const blocked_range<size_t> & r){
        size_t i = r.begin() == 0 ? 0 : last_indices[r.begin()-1].second; // index into coeffs / hfs to start with
        for(size_t k=r.begin(); k!=r.end(); ++k) {
            const std::pair<ObsId, size_t> & obs_li = last_indices[k];
            for(; i<obs_li.second; ++i){
                hfs[i].add_with_coeff_to(result[obs_li.first], coeffs[i](parameters), parameters);
            }
        }
    });
}

void tbb_model::get_prediction(DataWithUncertainties & result, const ParValues & parameters) const {
    get_prediction_impl<Histogram1DWithUncertainties>(result, parameters);
}

void tbb_model::get_prediction(Data & result, const ParValues & parameters) const {
    get_prediction_impl<Histogram1D>(result, parameters);
}

tbb_model::tbb_model(const Configuration & cfg): default_model(cfg), n_threads(get_int(cfg, "n_threads", task_scheduler_init::automatic)), tbb_init(n_threads){
}

REGISTER_PLUGIN(tbb_model)
