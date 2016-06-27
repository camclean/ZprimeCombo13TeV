#include "interface/phys.hpp"
#include "interface/histogram.hpp"
#include "interface/random.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"
#include "interface/model.hpp"

#include "test/utils.hpp"

#include "libconfig/libconfig.h++"

#include <iostream>
#include <iomanip>

#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

namespace{
Histogram1DWithUncertainties apply(const HistogramFunction & hf, const ParValues & values){
    size_t nbins;
    double xmin, xmax;
    hf.get_histogram_dimensions(nbins, xmin, xmax);
    Histogram1DWithUncertainties result(nbins, xmin, xmax);
    hf.add_with_coeff_to(result, 1.0, values);
    return result;
}


// if the relative or(!) absolute difference is < eps
bool close(const Histogram1DWithUncertainties & lhs, const Histogram1DWithUncertainties & rhs, bool verbose = false, double eps = 1e-7){
    if(lhs.get_nbins() != rhs.get_nbins()){
        if(verbose){
            cout << "nbins different" << endl;
        }
        return false;
    }
    for(size_t i=0; i<lhs.get_nbins(); ++i){
        double diff = fabs(lhs.get(i) - rhs.get(i));
        if(diff < eps) continue;
        double reldiff = diff / max(fabs(lhs.get(i)), fabs(rhs.get(i)));
        if(reldiff < eps) continue;
        if(verbose){
            cout << "bin " << i << " has large difference in value: lhs=" << setprecision(20) << lhs.get(i) << " rhs=" << setprecision(20) << rhs.get(i) << endl;
        }
        return false;
    }
    for(size_t i=0; i<lhs.get_nbins(); ++i){
        double diff = fabs(lhs.get_uncertainty2(i) - rhs.get_uncertainty2(i));
        if(diff < eps) continue;
        double reldiff = diff / max(fabs(lhs.get_uncertainty2(i)), fabs(rhs.get_uncertainty2(i)));
        if(reldiff < eps) continue;
        if(verbose){
            cout << "bin " << i << " has large difference in uncertainty2: lhs=" << setprecision(20) << lhs.get_uncertainty2(i) << " rhs=" << setprecision(20) << rhs.get_uncertainty2(i) << endl;
        }
        return false;
    }
    return true;
}
    
    
}


BOOST_AUTO_TEST_SUITE(model_tests)

BOOST_AUTO_TEST_CASE(model0){
    BOOST_CHECKPOINT("model0 entry");
    load_core_plugins();
    utils::fill_theta_dir(0);
    
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParIds pars;
    ObsIds obs;
    const size_t nbins = 100;
    ParId beta1 = vm->create_par_id("beta1");
    ParId beta2 = vm->create_par_id("beta2");
    ObsId obs0 = vm->create_obs_id("obs0", nbins, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("flat-histo = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "gauss-histo = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.5; normalize_to = 1.0;};\n"
            "c1 = {type = \"multiply\"; factors=(\"beta1\");};\n"
            "c2 = {type = \"multiply\"; factors=(\"beta2\");};\n"
            "dist-flat = {\n"
            "       type = \"flat_distribution\";\n"
            "       beta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; }; \n"
            "       beta2 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            " };\n"
            "m = {\n"
            "  obs0 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c1\";\n"
            "          histogram = \"@gauss-histo\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c2\";\n"
            "           histogram = \"@flat-histo\";\n"
            "       };\n"
            "   };\n"
            "  parameter-distribution = \"@dist-flat\";\n"
            "};\n"
            , vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model");
    std::auto_ptr<Model> m;
    try{
        m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(Exception & ex){
        std::cerr << ex.message << endl;
    }
    catch(libconfig::SettingNotFoundException & ex){
        std::cerr << ex.getPath() << " not found" << endl;
    }
    
    BOOST_REQUIRE(m.get()!=0);
    
    BOOST_CHECKPOINT("building signal histo");
    std::auto_ptr<HistogramFunction> f_signal_histo = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["gauss-histo"]));
    BOOST_CHECKPOINT("building bkg histo");
    std::auto_ptr<HistogramFunction> f_bkg_histo = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["flat-histo"]));
    
    ParValues values;
    Histogram1DWithUncertainties signal = apply(*f_signal_histo,values);
    Histogram1DWithUncertainties background = apply(*f_bkg_histo, values);
    
    values.set(beta1, 1.0);
    values.set(beta2, 0.0);
    Histogram1DWithUncertainties s;
    DataWithUncertainties pred;
    BOOST_CHECKPOINT("");
    m->get_prediction(pred, values);
    BOOST_CHECKPOINT("");
    s = pred[obs0];
    BOOST_CHECKPOINT("");
    //s should be signal only:
    for(size_t i = 0; i<nbins; i++){
        BOOST_REQUIRE(signal.get_value(i)==s.get_value(i));
    }
    //background only:
    values.set(beta1, 0.0);
    values.set(beta2, 1.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    for(size_t i = 0; i<nbins; i++){
        BOOST_REQUIRE(background.get_value(i)==s.get_value(i));
    }
    //zero prediction:
    values.set(beta1, 0.0);
    values.set(beta2, 0.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    for(size_t i = 0; i<nbins; i++){
        BOOST_REQUIRE(0.0==s.get_value(i));
    }

    //The likelihood, take double background. Use average as data:
    values.set(beta1, 1.0);
    values.set(beta2, 2.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    Data data;
    data[obs0] = s.get_values_histogram();
    std::auto_ptr<NLLikelihood> nll = m->get_nllikelihood(data);
    double x[2];
    x[0] = 0.9;
    x[1] = 1.9;
    double nll09 = (*nll)(x);
    x[0] = 1.0;
    x[1] = 2.0;
    double nll10 = (*nll)(x);
    x[0] = 1.1;
    x[1] = 2.1;
    double nll11 = (*nll)(x);
    BOOST_CHECKPOINT("check4");
    BOOST_CHECK(nll10 < nll11);
    BOOST_CHECK(nll10 < nll09);
    
    //test the model derivative:
    values.set(beta1, 1.0).set(beta2, 1.0);
    map<ParId, Histogram1D> ders;
    Histogram1D h1d;
    m->get_prediction_with_derivative(obs0, h1d, ders, values);
    BOOST_REQUIRE_EQUAL(h1d.get_nbins(), nbins);
    BOOST_REQUIRE_EQUAL(ders[beta1].get_nbins(), nbins);
    BOOST_REQUIRE_EQUAL(ders[beta2].get_nbins(), nbins);
    // the derivative w.r.t. beta should be the signal and background histos:
    for(size_t i=0; i<nbins; ++i){
        BOOST_CHECK_EQUAL(ders[beta1].get(i), signal.get(i));
        BOOST_CHECK_EQUAL(ders[beta2].get(i), background.get(i));
        BOOST_CHECK_EQUAL(h1d.get(i), signal.get(i) + background.get(i));
    }
    
    m->get_prediction(data, values);
    nll = m->get_nllikelihood(data);
    ParValues der;
    double nll0 = nll->eval_with_derivative(values, der);
    BOOST_CHECK(der.contains(beta1));
    BOOST_CHECK(der.contains(beta2));
    double nll1 = (*nll)(values);
    BOOST_CHECK_EQUAL(nll0, nll1);
    BOOST_CHECK_EQUAL(der.get(beta1), 0.0);
    BOOST_CHECK_EQUAL(der.get(beta2), 0.0);
    values.set(beta1, 1.1);
    nll->eval_with_derivative(values, der);
    BOOST_CHECK(der.get(beta1)!=0.0);
    BOOST_CHECK(der.get(beta2)!=0.0);
}

BOOST_AUTO_TEST_CASE(model_unc){
    load_core_plugins();
    utils::fill_theta_dir(0);
    
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParIds pars;
    ObsIds obs;
    const size_t nbins = 4;
    ParId beta1 = vm->create_par_id("beta1");
    ParId beta2 = vm->create_par_id("beta2");
    ParId beta3 = vm->create_par_id("beta3");
    ObsId obs0 = vm->create_obs_id("obs0", nbins, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc("flat-histo = {type = \"direct_data_histo\"; nbins = 4; range = (-1.0, 1.0); data = (100.0, 100.0, 50.0, 50.0); uncertainties = (10.0, 10.0, 5.0, 5.0); };\n"
            "c1 = {type = \"multiply\"; factors=(\"beta1\");};\n"
            "c2 = {type = \"multiply\"; factors=(\"beta2\");};\n"
            "c3 = {type = \"multiply\"; factors=(\"beta3\");};\n"
            "dist-flat = {\n"
            "       type = \"flat_distribution\";\n"
            "       beta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       beta2 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       beta3 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            " };\n"
            "m = {\n"
            "  obs0 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c1\";\n"
            "          histogram = \"@flat-histo\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c2\";\n"
            "           histogram = \"@flat-histo\";\n"
            "       };\n"
            "       background2 = {\n"
            "           coefficient-function = \"@c3\";\n"
            "           histogram = \"@flat-histo\";\n"
            "       };\n"
            "   };\n"
            "  parameter-distribution = \"@dist-flat\";\n"
            "  bb_uncertainties = true;\n"
            "};\n"
            , vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model");
    std::auto_ptr<Model> m;
    try{
        m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(Exception & ex){
        std::cerr << ex.message << endl;
    }
    catch(libconfig::SettingNotFoundException & ex){
        std::cerr << ex.getPath() << " not found" << endl;
    }
    
    ParValues values;
    
    values.set(beta1, 1.0);
    values.set(beta2, 1.0);
    values.set(beta3, 1.0);
    Histogram1DWithUncertainties s;
    DataWithUncertainties pred;
    m->get_prediction(pred, values);
    s = pred[obs0];
    //s should be signal only:
    BOOST_CHECK(close_to_relative(s.get(0), 300.0));
    BOOST_CHECK(close_to_relative(s.get(1), 300.0));
    BOOST_CHECK(close_to_relative(s.get(2), 150.0));
    BOOST_CHECK(close_to_relative(s.get(3), 150.0));
    
    BOOST_CHECK(close_to_relative(s.get_uncertainty(0), sqrt(3) * 10.0));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(1), sqrt(3) * 10.0));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(2), sqrt(3) * 5.0));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(3), sqrt(3) * 5.0));
    
    values.set(beta1, 2.0);
    m->get_prediction(pred, values);
    s = pred[obs0];
    BOOST_CHECK(close_to_relative(s.get(0), 400.0));
    BOOST_CHECK(close_to_relative(s.get(1), 400.0));
    BOOST_CHECK(close_to_relative(s.get(2), 200.0));
    BOOST_CHECK(close_to_relative(s.get(3), 200.0));
    
    BOOST_CHECK(close_to_relative(s.get_uncertainty(0), sqrt(20*20 + 2 * 10*10)));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(1), sqrt(20*20 + 2 * 10*10)));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(2), sqrt(10*10 + 2 * 5*5)));
    BOOST_CHECK(close_to_relative(s.get_uncertainty(3), sqrt(10*10 + 2 * 5*5)));
    
    // check model derivative:
    values.set(beta1, 1.0).set(beta2, 1.0).set(beta3, 1.0);
    map<ParId, Histogram1DWithUncertainties> ders;
    m->get_prediction_with_derivative(obs0, s, ders, values);
    BOOST_REQUIRE_EQUAL(ders.size(), 3);
    // The derivative values should be eual to the flat-histo, the derivative uncertainties-squared should
    // be 2*flat_histo.unc2
    const Histogram1DWithUncertainties flat_histo = get_constant_histogram(Configuration(cfg, cfg.setting["flat-histo"]));
    double beta_value = 1.0;
    Histogram1DWithUncertainties der_expected(flat_histo);
    for(size_t i=0; i<nbins; ++i){
        der_expected.set_unc2(i, flat_histo.get(i), flat_histo.get_uncertainty2(i) * 2.0 * beta_value);
    }
    BOOST_CHECK(close(der_expected, ders[beta1]));
    BOOST_CHECK(close(der_expected, ders[beta2]));
    BOOST_CHECK(close(der_expected, ders[beta3]));
    
    //some non-trivial point:
    beta_value = 1.7;
    values.set(beta1, beta_value);
    m->get_prediction_with_derivative(obs0, s, ders, values);
    der_expected = flat_histo;
    for(size_t i=0; i<nbins; ++i){
        der_expected.set_unc2(i, flat_histo.get(i), flat_histo.get_uncertainty2(i) * 2.0 * beta_value);
    }
    BOOST_CHECK(close(der_expected, ders[beta1]));
    
    // check nll derivative:
    values.set(beta1, 1.0).set(beta2, 1.0).set(beta3, 1.0);
    Data adata;
    m->get_prediction(adata, values);
    std::auto_ptr<NLLikelihood> nll = m->get_nllikelihood(adata);
    double nll0, nllp, h;
    nll0 = (*nll)(values);
    ParValues der;
    double nll0_2 = nll0 = nll->eval_with_derivative(values, der);
    BOOST_CHECK_EQUAL(nll0, nll0_2);
    BOOST_REQUIRE(der.contains(beta1) && der.contains(beta2) && der.contains(beta3));
    
    const double s_eps = 1e-6;
    values.set(beta1, 1.0 + s_eps);
    h = (1.0 + s_eps) - 1.0;
    nllp = (*nll)(values);
    double numder = (nllp - nll0) / h;
    cout << numder << " " << der.get(beta1) << endl;
    // derivative at the minimum should be very close to 0.0:
    BOOST_CHECK_SMALL(der.get(beta1), 1e-8);
    BOOST_CHECK_SMALL(der.get(beta2), 1e-8);
    BOOST_CHECK_SMALL(der.get(beta3), 1e-8);
    BOOST_CHECK_SMALL(numder - der.get(beta1), 1e-4);
    
    values.set(beta1, 1.0);
    values.set(beta2, 1.0 + s_eps);
    nllp = (*nll)(values);
    numder = (nllp - nll0) / h;
    cout << numder << " " << der.get(beta2) << endl;
    BOOST_CHECK_SMALL(numder - der.get(beta2), 1e-4);
    
    // go away from the minimum:
    double v0 = 1.5;
    values.set(beta1, v0);
    nll0 = (*nll)(values);
    values.set(beta1, v0 + s_eps);
    nllp = (*nll)(values);
    h = (v0 + s_eps) - v0;
    numder = (nllp - nll0) / h;
    nll->eval_with_derivative(values, der);
    BOOST_CHECK_CLOSE(der.get(beta1), numder, 1e-4);
    cout << numder << " " << der.get(beta1) << endl;

}

BOOST_AUTO_TEST_SUITE_END()
