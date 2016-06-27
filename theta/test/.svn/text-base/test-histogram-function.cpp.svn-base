#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables.hpp"
#include "test/utils.hpp"
#include <vector>
#include <iostream>

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
}

BOOST_AUTO_TEST_SUITE(histogram_function_tests)

BOOST_AUTO_TEST_CASE(root_histogram_range){
    bool loaded =  load_root_plugins();
    if(!loaded){
       std::cout << "In test root_histogram_range: root plugin not loaded, not executing test" << endl;
       return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ConfigCreator cc(
            "root-histo1 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\";};\n" // no range
            "root-histo2 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-4.0, 20.0); };\n" // full range
            "root-histo3 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-2.0, 2.0); };\n" // sub range
            "root-histo4 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-4.1, 20.1); };\n" // range excluding underflow / overflow
            "root-histo5 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\"; range = (-3.9, 19.9); };\n" // invalid range
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building hf");
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo1"]));
    ParValues pv;
    Histogram1DWithUncertainties h1 = apply(*hf,pv);
    BOOST_REQUIRE(h1.get_nbins()==24);
    BOOST_REQUIRE(h1.get_xmin()==-4);
    BOOST_REQUIRE(h1.get_xmax()==20);
    hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo2"]));
    Histogram1DWithUncertainties h2 = apply(*hf,pv);
    BOOST_REQUIRE(h2.get_nbins()==24);
    BOOST_REQUIRE(h2.get_xmin()==-4);
    BOOST_REQUIRE(h2.get_xmax()==20);
    for(int i=0; i<24; ++i){
        BOOST_CHECK(h1.get_value(i) == i+13.0);
        BOOST_CHECK(h1.get_value(i) == h2.get_value(i));
    }
    hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo3"]));
    Histogram1DWithUncertainties h3 = apply(*hf,pv);
    BOOST_REQUIRE(h3.get_nbins()==4);
    BOOST_REQUIRE(h3.get_xmin()==-2);
    BOOST_REQUIRE(h3.get_xmax()==2);
    for(int i=0; i<4; ++i){
       BOOST_CHECK(h3.get_value(i) == i+15.0);
    }
    BOOST_CHECKPOINT("building h4");
    hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo4"]));
    Histogram1DWithUncertainties h4 = apply(*hf,pv);
    BOOST_REQUIRE(h4.get_nbins()==26);
    BOOST_REQUIRE(h4.get_xmin()==-5);
    BOOST_REQUIRE(h4.get_xmax()==21);
    for(int i=0; i<26; ++i){
       BOOST_CHECK(h4.get_value(i) == i+12);
    }
    
    bool except = false;
    try{
        hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo5"]));
    }
    catch(ConfigurationException & ex){
       except = true;
    }
    BOOST_CHECK(except);
}


BOOST_AUTO_TEST_CASE(root_histogram){
    bool loaded =  load_root_plugins();
    if(!loaded){
       std::cout << "In test root_histogram: root plugin not loaded, not executing test" << endl;
       return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ConfigCreator cc(
            "root-histo1 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo1d\";};\n"
            "root-histo2 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo2d\";};\n"
            "root-histo3 = {type = \"root_histogram\"; filename=\"testhistos.root\"; histoname = \"histo3d\";};\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building hf");
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo1"]));
    ParValues pv;
    Histogram1DWithUncertainties h = apply(*hf,pv);
    BOOST_REQUIRE(h.get_nbins()==24);
    for(size_t i=0; i<h.get_nbins(); ++i){
        BOOST_ASSERT(h.get_value(i)==i+13);
    }
    //2D histogram:
    hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo2"]));
    h = apply(*hf,pv);
    BOOST_REQUIRE(h.get_nbins()==10*11);
    //calculate the expected integral (excluding overflow / underflow!!);
    // if this matches, we are satisfied and believe the rest is Ok as well:
    double expected_integral = 0.0;
    for(int i=1; i<=10; ++i){
       for(int j=1; j<=11; ++j){
          expected_integral += (i + 0.78) * (j + 3.02);
       }
    }
    //BOOST_ASSERT(close_to_relative(h.get_sum(),expected_integral));
    //3D histogram, same as 2D:
    hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["root-histo3"]));
    h = apply(*hf,pv);
    BOOST_REQUIRE(h.get_nbins()==10*11*12);
    //calculate the expected integral (excluding overflow / underflow!!);
    // if this matches, we are satisfied and believe the rest is Ok as well:
    expected_integral = 0.0;
    for(int i=1; i<=10; ++i){
       for(int j=1; j<=11; ++j){
          for(int k=1; k<=12; ++k){
             expected_integral += (i + 0.12) * (j + 1.34) * (k + 5.67);
          }
       }
    }
    //BOOST_ASSERT(close_to_relative(h.get_sum(),expected_integral));
}


BOOST_AUTO_TEST_CASE(cubiclinear_histomorph){
    load_core_plugins();
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    
    const size_t nbins = 1;
    ParId delta = vm->create_par_id("delta");
    vm->create_obs_id("obs", nbins, -1, 1);
    BOOST_CHECKPOINT("parsing config");
    ConfigCreator cc("flat-histo0 = {type = \"direct_data_histo\"; nbins = 1; range = [-1.0, 1.0]; data = [1.0]; uncertainties = [0.1]; };\n"
            "flat-histo1 = {type = \"direct_data_histo\"; nbins = 1; range = [-1.0, 1.0]; data = [1.12]; };\n"
            "flat-histo-1 = {type = \"direct_data_histo\"; nbins = 1; range = [-1.0, 1.0]; data = [0.83]; };\n"
            "histo = { type = \"cubiclinear_histomorph\"; parameters = (\"delta\"); nominal-histogram = \"@flat-histo0\";\n"
            "      delta-plus-histogram = \"@flat-histo1\"; delta-minus-histogram = \"@flat-histo-1\"; };\n"
            "histo2 = { type = \"cubiclinear_histomorph\"; parameters = (\"delta\"); nominal-histogram = \"@flat-histo0\";\n"
            "      delta-plus-histogram = \"@flat-histo1\"; delta-minus-histogram = \"@flat-histo-1\"; normalize_to_nominal = true; };\n"
            , vm);
            
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building hf");
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["histo"]));
    BOOST_CHECKPOINT("hf built");
    std::auto_ptr<HistogramFunction> hf_nominal = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["flat-histo0"]));
    ParValues pv;
    Histogram1DWithUncertainties flat0 = apply(*hf_nominal, pv);
    BOOST_CHECK(flat0.get_nbins()==1);
    BOOST_CHECK(flat0.get_value(0)==1.0);
    //check central and +- 1 sigma values:
    pv.set(delta, 0.0);
    Histogram1DWithUncertainties h = apply(*hf, pv);
    BOOST_REQUIRE(h.get_nbins()==1);
    BOOST_REQUIRE(h.get_xmin()==-1);
    BOOST_REQUIRE(h.get_xmax()==1);
    BOOST_CHECK(close_to_relative(h.get_value(0), 1.0));
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1));
    pv.set(delta, 1.0);
    h = apply(*hf,pv);
    BOOST_CHECK(close_to_relative(h.get_value(0), 1.12));
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1));
    pv.set(delta, -1.0);
    h = apply(*hf,pv);
    BOOST_CHECK(close_to_relative(h.get_value(0), 0.83));
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1));
    //+- 2 sigma values, should be interpolated linearly:
    pv.set(delta, 2.0);
    h = apply(*hf,pv);
    BOOST_CHECK(close_to_relative(h.get_value(0), 1.24));
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1));
    pv.set(delta, -2.0);
    h = apply(*hf,pv);
    BOOST_CHECK(close_to_relative(h.get_value(0), 0.66));
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1));
    //cutoff at zero:
    pv.set(delta, -10.0);
    h = apply(*hf,pv);
    BOOST_ASSERT(h.get_value(0) == 0.0);
    //derivative at zero should be smooth:
    pv.set(delta, 1e-8);
    h = apply(*hf,pv);
    double eps = h.get_value(0);
    pv.set(delta, -1e-8);
    h = apply(*hf,pv);
    double eps_minus = h.get_value(0);
    BOOST_ASSERT(close_to(eps - 1, 1 - eps_minus, 1000.0));
    //derivative at zero should be (0.12 + 0.17) / 2.
    double der = (eps - eps_minus) / (2e-8);
    BOOST_ASSERT(fabs(der - (0.12 + 0.17)/2) < 1e-8);
    
    // check derivative:
    Histogram1D h1;
    map<ParId, Histogram1D> ders;
    pv.set(delta, 0.0);
    ders[delta].reset(nbins, -1, 1);
    hf->eval_and_add_derivatives(h1, ders, 1.0, pv);
    BOOST_CHECK_CLOSE(h1.get(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(ders[delta].get(0), (0.12 + 0.17)/2, 1e-10);
    // non-trivial coeff:
    ders[delta].reset();
    double coeff = M_PI;
    hf->eval_and_add_derivatives(h1, ders, coeff, pv);
    BOOST_CHECK_CLOSE(h1.get(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(ders[delta].get(0), (0.12 + 0.17)/2 * coeff, 1e-10);
    // non-trivial content before calling:
    ders[delta].reset();
    ders[delta].set(0, coeff);
    h1.set(0, 0.0);
    hf->eval_and_add_derivatives(h1, ders, coeff, pv);
    BOOST_CHECK_CLOSE(h1.get(0), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(ders[delta].get(0), coeff + (0.12 + 0.17)/2 * coeff, 1e-10);
    
    // non-trivial deltas:
    double dvals[] = {-10.0, -2.0, -1.1, -0.7, -0.2, 0.2, 0.8, 1.7, 5.0};
    const double s_eps = sqrt(numeric_limits<double>::epsilon());
    for(size_t i=0; i<sizeof(dvals) / sizeof(double); ++i){
        double dval = dvals[i];
        cout << "dval = " << dval << endl;
        pv.set(delta, dval);
        ders[delta].reset();
        hf->eval_and_add_derivatives(h1, ders, 1.0, pv);
        // get numerical derivative:
        double h0 = h1.get(0);
        double dval_plus = dval + s_eps;
        pv.set(delta, dval_plus);
        double hplus = apply(*hf, pv).get(0);
        der = (hplus - h0) / (dval_plus - dval);
        BOOST_CHECK_CLOSE(der, ders[delta].get(0), 1e-5);
    }
        
    // h2:
    std::auto_ptr<HistogramFunction> hf2 = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["histo2"]));
    pv.set(delta, 0.0);
    h = apply(*hf2,pv);
    BOOST_CHECK(h.get_value(0) == 1.0);
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1));
    
    pv.set(delta, 1.0); // should interpolate to 1.12, so uncertainty after re-normalization should shrink by that factor ...
    h = apply(*hf2,pv);
    BOOST_CHECK(h.get_value(0) == 1.0);
    BOOST_CHECK(close_to_relative(h.get_uncertainty(0), 0.1 / 1.12));
}



BOOST_AUTO_TEST_SUITE_END()
