#include "interface/decls.hpp"
#include "interface/model.hpp"
#include "interface/data.hpp"

#include "tbb/task_scheduler_init.h"


/** \brief Model using Intel's Threading Building Blocks for parallelizing the evaluation on channel-level
 * 
 * It inherits from the \c default_model, see documentation there. The only additional
 * (optional) configuration setting is an integer \c n_threads which is the number of threads to use. If not given,
 * the default from the tbb library is used (which uses a machine-dependent value).
 * 
 * The parallelization is done on the level of channels, i.e., the prediction is evaluated in parallel
 * for several observables / channels. Therefore, using this parallelization is beneficial especially in cases
 * of many "complicated" observables (with many processes / uncertainties / bins).
 */
class tbb_model: public theta::default_model{
private:
    template<typename HT>
    void get_prediction_impl(theta::DataT<HT> & result, const theta::ParValues & parameters) const;
    
    int n_threads;
    tbb::task_scheduler_init tbb_init;
    
public:
    explicit tbb_model(const theta::Configuration & cfg);
    virtual void get_prediction(theta::DataWithUncertainties & result, const theta::ParValues & parameters) const;
    virtual void get_prediction(theta::Data & result, const theta::ParValues & parameters) const;  
};

