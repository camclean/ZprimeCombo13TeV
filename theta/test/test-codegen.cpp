#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/random.hpp"
#include "interface/histogram-function.hpp"
#include "test/utils.hpp"

#include "libconfig/libconfig.h++"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace theta;

namespace{
    
bool operator==(const Histogram1DWithUncertainties& h0, const Histogram1DWithUncertainties & h1){
    if(h0.get_nbins()!=h1.get_nbins()) return false;
    if(h0.get_xmin()!=h1.get_xmin()) return false;
    if(h0.get_xmax()!=h1.get_xmax()) return false;
    for(size_t i=0; i<h0.get_nbins(); ++i){
        if(h0.get_value(i)!=h1.get_value(i)) return false;
        if(h0.get_uncertainty(i)!=h1.get_uncertainty(i)) return false;
    }
    return true;
}

bool operator==(const Histogram1D& h0, const Histogram1D & h1){
    if(h0.get_nbins()!=h1.get_nbins()) return false;
    if(h0.get_xmin()!=h1.get_xmin()) return false;
    if(h0.get_xmax()!=h1.get_xmax()) return false;
    for(size_t i=0; i<h0.get_nbins(); ++i){
        if(h0.get(i)!=h1.get(i)) return false;
    }
    return true;
}

using ::close_to;
bool close_to(const Histogram1D& h0, const Histogram1D & h1, double scale = 1.0, bool print = false){
	h0.check_compatibility(h1);
	for(size_t i=0; i<h0.get_nbins(); ++i){
		if(!close_to(h0.get(i), h1.get(i), scale)){
			if(print){
				int k = i;
				printf("bin %d different: h0[%d]=%.20g;  h1[%d]=%.20g;  diff=%.5g\n", k, k, h0.get(i), k, h1.get(i), h0.get(i) - h1.get(i));
			}
			return false;
		}
	}
	return true;
}

bool close_to(const Histogram1DWithUncertainties& h0, const Histogram1DWithUncertainties & h1, double scale = 1.0, bool print = false){
	return close_to(h0.get_values_histogram(), h1.get_values_histogram(), scale, print) &&
		   close_to(h0.get_uncertainty2_histogram(), h1.get_uncertainty2_histogram(), scale, print);
}

bool close_to(const DataWithUncertainties & d1, const DataWithUncertainties & d2, double scale = 1.0, bool print = false){
	if(d1.get_observables() != d2.get_observables()){
		if(print){
			printf("observables are not the same!\n");
		}
		return false;
	}
	ObsIds obs = d1.get_observables();
	for(ObsIds::const_iterator it = obs.begin(); it!=obs.end(); ++it){
		if(!close_to(d1[*it], d2[*it], scale, print)) return false;
	}
	return true;
}

void dump(const Histogram1DWithUncertainties & h, const string & name){
	cout << "histogram " << name << "(" << h.get_nbins() << ", " << h.get_xmin() << ", " << h.get_xmax() << ") ";
	for(size_t i=0; i<h.get_nbins(); ++i){
		int k = i;
		printf("%d: %6.4g +- %6.4g\n", k, h.get(i), h.get_uncertainty(i));
	}
}

double cubic_interpolation_delta(double nom, double minus, double plus, double delta){
	double hplus = plus - nom;
	double hminus = minus - nom;
	double sum = hplus + hminus;
	double diff = hplus - hminus;
	double result = 0.5 * delta * diff;
	result += (delta * delta - 0.5 * pow(fabs(delta), 3)) * sum;
	return result;
}




}

BOOST_AUTO_TEST_SUITE(codegen)

BOOST_AUTO_TEST_CASE(llvm_function_wrapper){
	load_core_plugins();
    if(!load_llvm_plugins()){
        std::cout << "skipping llvm test" << endl;
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    vm->create_par_id("p2");
    ConfigCreator cc(
            "f0 = {type = \"exp_function\"; parameters = (\"p1\"); lambdas_plus = (0.2); lambdas_minus = (0.1); };\n"
            "llvm_f0 = {type = \"llvm_enable_function\"; function = \"@f0\";};\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function");
    std::auto_ptr<Function> f0 = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f0"]));
    std::auto_ptr<Function> llvm_f0;
    try{
         llvm_f0 = PluginManager<Function>::build(Configuration(cfg, cfg.setting["llvm_f0"]));
    }
    catch(Exception & ex){
        std::cout << ex.message << endl;
        throw;
    }
    ParValues vals;
    vals.set(p1, 1.7);
    BOOST_CHECK((*f0)(vals) == (*llvm_f0)(vals));
    vals.set(p1, -1.8);
    BOOST_CHECK((*f0)(vals) == (*llvm_f0)(vals));
}


BOOST_AUTO_TEST_CASE(llvm_multiply){
	load_core_plugins();
    if(!load_llvm_plugins()){
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ConfigCreator cc(
            "f0 = {type = \"exp_function\"; parameters = (\"p1\"); lambdas_plus = (0.2); lambdas_minus = (0.1); };\n"
            "mul = {type = \"multiply\"; factors = (\"@f0\", \"p2\", 1.8); };\n"
            "llvm_mul = {type = \"llvm_multiply\"; factors = \"@mul.factors\";};\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function");
    std::auto_ptr<Function> mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["mul"]));
    std::auto_ptr<Function> llvm_mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["llvm_mul"]));
    ParValues vals;
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);//no setSeed to keep it the same every time ...
    
    for(size_t i=0; i<30; ++i){
        vals.set(p1, rnd.uniform() * 10 - 5);
        vals.set(p2, rnd.uniform() * 10 - 5);
        BOOST_CHECK((*mul)(vals) == (*llvm_mul)(vals));
    }
}


BOOST_AUTO_TEST_CASE(llvm_exp_function){
	load_core_plugins();
    if(!load_llvm_plugins()){
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ConfigCreator cc(
            "f = {type = \"exp_function\"; parameters = (\"p1\", \"p2\"); lambdas_plus = (0.1, 0.2); lambdas_minus = (0.1, 0.18); };\n"
            "llvm_f = {type = \"llvm_exp_function\"; parameters = \"@f.parameters\"; lambdas_plus = \"@f.lambdas_plus\"; lambdas_minus = \"@f.lambdas_minus\"; };\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building function mul");
    std::auto_ptr<Function> mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f"]));
    BOOST_CHECKPOINT("building function llvm_mul");
    std::auto_ptr<Function> llvm_mul = PluginManager<Function>::build(Configuration(cfg, cfg.setting["llvm_f"]));
    ParValues vals;
    
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);//no setSeed to keep it the same every time ...
    
    for(size_t i=0; i<30; ++i){
        vals.set(p1, rnd.uniform() * 10 - 5);
        vals.set(p2, rnd.uniform() * 10 - 5);
        BOOST_CHECK((*mul)(vals) == (*llvm_mul)(vals));
    }
}


BOOST_AUTO_TEST_CASE(constant_histo){
	load_core_plugins();
    if(!load_llvm_plugins()){
        return;
    }
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    utils::fill_theta_dir(0);
    ConfigCreator cc(
            "histo = {type = \"direct_data_histo\"; range = (0.0, 3.0); nbins = 3; data = (17.0, 23.0, 47.7); };\n"
            "llvm_h = {type = \"llvm_enable_histogram_function\"; histogram_function = \"@histo\"; debug = true; };\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building histo");
    std::auto_ptr<HistogramFunction> histo = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["histo"]));
    BOOST_CHECKPOINT("building llvm_h");
    std::auto_ptr<HistogramFunction> llvm_h = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["llvm_h"]));
    BOOST_CHECKPOINT("evaluating");
    Histogram1DWithUncertainties h;
    histo->apply_functor(copy_to<Histogram1DWithUncertainties>(h), ParValues());
    Histogram1DWithUncertainties h_l;
    llvm_h->apply_functor(copy_to<Histogram1DWithUncertainties>(h_l), ParValues());
    BOOST_REQUIRE(h.get_nbins() == h_l.get_nbins());
    for(size_t i=0; i<h_l.get_nbins(); ++i){
        BOOST_REQUIRE(h.get_value(i) == h_l.get_value(i));
    }
}

BOOST_AUTO_TEST_CASE(cubiclinear_histo){
	BOOST_CHECKPOINT("cubiclinear enter");
	load_core_plugins();
	if(!load_llvm_plugins()){
		return;
	}
	boost::shared_ptr<VarIdManager> vm(new VarIdManager);
	ParId p = vm->create_par_id("p");
	utils::fill_theta_dir(0);
	ConfigCreator cc(
			"histo = {type = \"direct_data_histo\"; range = (0.0, 3.0); nbins = 3; data = (10.0, 10.1, 10.2); };\n"
			"histop = {type = \"direct_data_histo\"; range = (0.0, 3.0); nbins = 3; data = (12.0, 12.1, 15.2); };\n"
			"histom = {type = \"direct_data_histo\"; range = (0.0, 3.0); nbins = 3; data = (7.0, 8.1, 9.2); };\n"
			"llvm_hf = {type = \"llvm_cubiclinear_histomorph\"; nominal-histogram = \"@histo\"; parameters = (\"p\");\n"
			"    p-plus-histogram = \"@histop\"; p-minus-histogram = \"@histom\";};\n"
			"llvm_hf2 = {type = \"llvm_enable_histogram_function\"; histogram_function = \"@llvm_hf\";};\n"
			"hf = {type = \"cubiclinear_histomorph\"; nominal-histogram = \"@histo\"; parameters = (\"p\");\n"
			"    p-plus-histogram = \"@histop\"; p-minus-histogram = \"@histom\";};\n"
			, vm);
	const theta::Configuration & cfg = cc.get();
	BOOST_CHECKPOINT("building hf");
	std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["llvm_hf2"]));
	BOOST_CHECKPOINT("hf built");
	std::auto_ptr<HistogramFunction> hf_reference = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting["hf"]));

	Histogram1D hnominal = get_constant_histogram(Configuration(cfg, cfg.setting["histo"])).get_values_histogram();

	Histogram1D result, expected;
	ParValues values;
	values.set(p, 0.0);
	hf->apply_functor(copy_to<Histogram1D>(result), values);
	BOOST_CHECK(result == hnominal);
	result.set_all_values(-1.0);

	// linear extrapolation:
	values.set(p, 1.1);
	hf->apply_functor(copy_to<Histogram1D>(result), values);
	hf_reference->apply_functor(copy_to<Histogram1D>(expected), values);
	BOOST_CHECK(close_to(result, expected));

	values.set(p, 3.1);
	hf->apply_functor(copy_to<Histogram1D>(result), values);
	hf_reference->apply_functor(copy_to<Histogram1D>(expected), values);
	BOOST_CHECK(close_to(result, expected));

	double delta = 1.7;
	values.set(p, -delta);
	hf->apply_functor(copy_to<Histogram1D>(result), values);
	hf_reference->apply_functor(copy_to<Histogram1D>(expected), values);
	BOOST_CHECK(close_to(result, expected));

	//cubic interpolation
	for(delta=-0.33; delta < 0.9; delta += .1){
		values.set(p, delta);
		hf->apply_functor(copy_to<Histogram1D>(result), values);
		hf_reference->apply_functor(copy_to<Histogram1D>(expected), values);
		BOOST_CHECK(close_to(result, expected));
	}

	// * truncation at 0.0
	for(delta = -4.0; delta > -18.; delta -= 6.0){
		values.set(p, delta);
		hf->apply_functor(copy_to<Histogram1D>(result), values);
		hf_reference->apply_functor(copy_to<Histogram1D>(expected), values);
		BOOST_CHECK(close_to(result, expected));
	}

	// * parameter_factors
	// * more than 1 syst

	// * FIXME: normalize_to_nominal
}


// compare model prediction of llvm model to classical model; including likelihood.
BOOST_AUTO_TEST_CASE(model){
	load_core_plugins();
    if(!load_llvm_plugins()){
        return;
    }
    utils::fill_theta_dir(0);
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParIds pars;
    ObsIds obs;
    const size_t nbins0 = 101;
    const size_t nbins1 = 50;
    ParId beta1 = vm->create_par_id("beta1");
    ParId beta2 = vm->create_par_id("beta2");
    ParId delta0 = vm->create_par_id("delta0");
    ParId delta1 = vm->create_par_id("delta1");
    ObsId obs0 = vm->create_obs_id("obs0", nbins0, -1, 1);
    ObsId obs1 = vm->create_obs_id("obs1", nbins1, -1, 1);
    pars.insert(beta1);
    pars.insert(beta2);
    obs.insert(obs0);
    
    //boost::ptr_vector<Function> coeffs;
    BOOST_CHECKPOINT("parsing config");
    
    ConfigCreator cc(
            "flat-histo0 = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0;};\n"
            "gauss-histo0 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.5; normalize_to = 1.0;};\n"
            "flat-histo1 = {type = \"fixed_poly\"; observable=\"obs1\"; coefficients = [1.0]; normalize_to = 1.7;};\n"
            "gauss-histo1 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.23; mean = 0.43; normalize_to = 1.12;};\n"
            "c-gauss = {type = \"multiply\"; factors=(1.1, \"beta1\", {type = \"exp_function\"; parameters = (\"delta0\", \"delta1\"); lambdas_plus = (0.1, 0.12); lambdas_minus = (0.1, 0.13);});};\n"
            "c-flat = {type = \"multiply\"; factors=(\"beta2\");};\n"
            "dist-flat = {\n"
            "       type = \"flat_distribution\";\n"
            "       beta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       beta2 = { range = (\"-inf\", \"inf\"); fix-sample-value = 1.0; };\n"
            "       delta0 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
            "       delta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
            " };\n"
            "m = {\n"
            "  obs0 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c-gauss\";\n"
            "          histogram = \"@gauss-histo0\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c-flat\";\n"
            "           histogram = \"@flat-histo0\";\n"
            "       };\n"
            "  };"
            "  obs1 = {\n"
            "       signal = {\n"
            "          coefficient-function = \"@c-gauss\";\n"
            "          histogram = \"@gauss-histo1\";\n"
            "       };\n"
            "       background = {\n"
            "           coefficient-function = \"@c-flat\";\n"
            "           histogram = \"@flat-histo1\";\n"
            "       };\n"
            "   };\n"
            "  parameter-distribution = \"@dist-flat\";\n"
            "};\n"
            "\n"
            "llvm_m = {\n"
            "   type = \"llvm_model\";\n"
            "   llvm_always = true;\n"
            "   obs0 = \"@m.obs0\";\n"
            "   obs1 = \"@m.obs1\";\n"
            "   parameter-distribution = \"@m.parameter-distribution\";\n"
            "};"
            , vm);
    
    BOOST_CHECKPOINT("config parsed");
    
    const theta::Configuration & cfg = cc.get();
    BOOST_CHECKPOINT("building model m");
    std::auto_ptr<Model> m;
    try{
       m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
    }
    catch(libconfig::SettingNotFoundException & ex){
        std::cout << ex.getPath() << endl;
        throw;
    }
    BOOST_CHECKPOINT("building model llvm_m");
    std::auto_ptr<Model> llvm_m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["llvm_m"]));
    
    BOOST_REQUIRE(m->get_parameters() == llvm_m->get_parameters());
    BOOST_REQUIRE(m->get_observables() == llvm_m->get_observables());
        
    Data d;
    Data llvm_d;
    
    // calculate some predictions:
    ParValues vals;
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);//no setSeed to keep it the same every time ...
    
    BOOST_CHECKPOINT("starting with random tests");
    for(size_t i=0; i<50; ++i){
        vals.set(beta1, rnd.uniform() * 5);
        vals.set(beta2, rnd.uniform() * 5);
        vals.set(delta0, rnd.uniform() * 10 - 5);
        vals.set(delta1, rnd.uniform() * 10 - 5);
        m->get_prediction(d, vals);
        llvm_m->get_prediction(llvm_d, vals);
        BOOST_REQUIRE(d.get_observables().contains(obs0) && d.get_observables().contains(obs1));
        BOOST_REQUIRE(d.get_observables() == llvm_d.get_observables());
        BOOST_CHECK(d[obs0] == llvm_d[obs0]);
        BOOST_CHECK(d[obs1] == llvm_d[obs1]);
    }
    
    
    // calculate some nll values:
    vals.set(beta1, 1.17);
    vals.set(beta2, 0.83);
    vals.set(delta0, 0.12);
    vals.set(delta1, -0.18);    
    std::auto_ptr<NLLikelihood> nll = m->get_nllikelihood(d);
    std::auto_ptr<NLLikelihood> llvm_nll = llvm_m->get_nllikelihood(d);
    for(size_t i=0; i<50; ++i){
        vals.set(beta1, rnd.uniform() * 5);
        vals.set(beta2, rnd.uniform() * 5);
        vals.set(delta0, rnd.uniform() * 10 - 5);
        vals.set(delta1, rnd.uniform() * 10 - 5);
        double nll_value = (*nll)(vals);
        double llvm_nll_value = (*llvm_nll)(vals);
        //cout << i << " " << nll_value << " " << llvm_nll_value << " diff = " << fabs(nll_value - llvm_nll_value) << endl;
        // note: the values are not completely equal, as llvm_nll calculates the template likelihood over ALL bins in one
        // go and nll calculates it for each observable seperately. The differences can only be at most N_obs roundings, though
        BOOST_CHECK(close_to_relative(nll_value, llvm_nll_value));
    }
}


// compare model prediction of llvm model to classical model with barlow-beeston uncertainties. Also compare likelihood values.
BOOST_AUTO_TEST_CASE(model_bb_unc){
	load_core_plugins();
    if(!load_llvm_plugins()){
        return;
    }
    utils::fill_theta_dir(0);
    // use a model with 2 observables and 2 real-valued observables.
    // obs0 has two components: signal and background; signal uses template morphing with 2 parameters (delta0, delta1),
    //    no stat. unc., normalize_to_nominal = true., parameter-factors = (1.0, 0.5)
    // obs1 has two components: signal and background; signal uses template morphing with 2 parameters (delta0, delta1), with stat. unc.,
    //     normalize_to_nominal = false; parameter-factors = (0.5, -0.5)
    //
    // the real-values observables rvob0 and rvobs1 have gaussian distribution around rv_delta0, rv_delta1

    // use different number of bins as the implementation is nbin-dependent. Make sure to use even and odd, <=8 and >8 and ==1.
    const int N = 5;
    size_t vnbins0[N] = {1, 8, 16, 27, 7};
    size_t vnbins1[N] = {50, 31, 6, 5, 1};

    for(int itest=0; itest<1; ++itest){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
	ParIds pars;
	ObsIds obs;
	const size_t nbins0 = vnbins0[itest];
	const size_t nbins1 = vnbins1[itest];
	// the template morphing parameters for "signal" component:
	ParId delta0 = vm->create_par_id("delta0");
	ParId delta1 = vm->create_par_id("delta1");
	ParId rv_delta0 = vm->create_par_id("rv_delta0");
	ParId rv_delta1 = vm->create_par_id("rv_delta1");
	ParId rvobs0 = vm->create_par_id("rvobs0", "rvobs");
	ParId rvobs1 = vm->create_par_id("rvobs1", "rvobs");
	ObsId obs0 = vm->create_obs_id("obs0", nbins0, -1, 1);
	ObsId obs1 = vm->create_obs_id("obs1", nbins1, -1, 1);

	BOOST_CHECKPOINT("parsing config");
	// FIXME: use uncertainties!

	ConfigCreator cc(
			    // obs0:
	            "flat-histo0 = {type = \"fixed_poly\"; observable=\"obs0\"; coefficients = [1.0]; normalize_to = 1.0; };\n"
	            "gauss-histo0 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.0; normalize_to = 1.0; };\n"
				"gauss-histo0-p0 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.6; mean = 0.0; normalize_to = 1.1; };\n"
				"gauss-histo0-m0 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.4; mean = 0.0; normalize_to = 0.9; };\n"
				"gauss-histo0-p1 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = 0.1; normalize_to = 1.0; };\n"
				"gauss-histo0-m1 = {type = \"fixed_gauss\"; observable=\"obs0\"; width = 0.5; mean = -0.1; normalize_to = 1.0; };\n"
				""
				// obs1:
	            "flat-histo1 = {type = \"fixed_poly\"; observable=\"obs1\"; coefficients = [1.0]; normalize_to = 1.7; relative_bb_uncertainty = 0.13;};\n"
	            "gauss-histo1 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.23; mean = 0.43; normalize_to = 1.12;  relative_bb_uncertainty = 0.2;};\n"
				"gauss-histo1-p0 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.25; mean = 0.43; normalize_to = 1.3;};\n"
				"gauss-histo1-m0 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.20; mean = 0.43; normalize_to = 1.0;};\n"
				"gauss-histo1-p1 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.23; mean = 0.47; normalize_to = 1.12;};\n"
				"gauss-histo1-m1 = {type = \"fixed_gauss\"; observable=\"obs1\"; width = 0.23; mean = 0.38; normalize_to = 1.12;};\n"
				""
				""
				"obs0-signal-histo = { type = \"cubiclinear_histomorph\"; parameters = (\"delta0\", \"delta1\"); parameter-factors = (1.0, 0.5);\n"
				"  nominal-histogram = \"@gauss-histo0\"; delta0-plus-histogram = \"@gauss-histo0-p0\"; delta0-minus-histogram = \"@gauss-histo0-m0\";\n"
				"  delta1-plus-histogram = \"@gauss-histo0-p1\"; delta1-minus-histogram = \"@gauss-histo0-m1\"; normalize_to_nominal = true;};\n"
				"obs1-signal-histo = { type = \"cubiclinear_histomorph\"; parameters = (\"delta0\", \"delta1\"); parameter-factors = (0.5, -0.5);\n"
							"  nominal-histogram = \"@gauss-histo1\"; delta0-plus-histogram = \"@gauss-histo1-p0\"; delta0-minus-histogram = \"@gauss-histo1-m0\";\n"
							"  delta1-plus-histogram = \"@gauss-histo1-p1\"; delta1-minus-histogram = \"@gauss-histo1-m1\";};\n"
				""
	            "c-gauss = {type = \"multiply\"; factors=(1000.);};\n"
	            "c-flat = {type = \"multiply\"; factors=(1000.);};\n"
	            "dist-flat = {\n"
	            "       type = \"flat_distribution\";\n"
	            "       delta0 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
	            "       delta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
				"       rv_delta0 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
				"       rv_delta1 = { range = (\"-inf\", \"inf\"); fix-sample-value = 0.0; };\n"
	            " };\n"
	            "m = {\n"
	            "  obs0 = {\n"
	            "       signal = {\n"
	            "          coefficient-function = \"@c-gauss\";\n"
	            "          histogram = \"@obs0-signal-histo\";\n"
	            "       };\n"
	            "       background = {\n"
	            "           coefficient-function = \"@c-flat\";\n"
	            "           histogram = \"@flat-histo0\";\n"
	            "       };\n"
	            "  };"
	            "  obs1 = {\n"
	            "       signal = {\n"
	            "          coefficient-function = \"@c-gauss\";\n"
	            "          histogram = \"@obs1-signal-histo\";\n"
	            "       };\n"
	            "       background = {\n"
	            "           coefficient-function = \"@c-flat\";\n"
	            "           histogram = \"@flat-histo1\";\n"
	            "       };\n"
	            "   };\n"
	            "  parameter-distribution = \"@dist-flat\";\n"
				"  rvobs-distribution = { type = \"product_distribution\";"
				"       distributions = ({type = \"gauss1d\"; mean = \"rv_delta0\"; parameter = \"rvobs0\"; width = 1.0; range = [\"-inf\", \"inf\"];},"
				"						 {type = \"gauss1d\"; mean = \"rv_delta1\"; parameter = \"rvobs1\"; width = 1.0; range = [\"-inf\", \"inf\"];});};"
				"   bb_uncertainties = true;"
	            "};\n"
	            "\n"
	            "llvm_m = {\n"
	            "   type = \"llvm_model\";\n"
	            "   llvm_always = true;\n"
	            "   obs0 = \"@m.obs0\";\n"
	            "   obs1 = \"@m.obs1\";\n"
	            "   parameter-distribution = \"@m.parameter-distribution\";\n"
				"   rvobs-distribution = \"@m.rvobs-distribution\";\n"
				"   bb_uncertainties = true;"
	            "};"
	            , vm);
	BOOST_CHECKPOINT("config parsed");
	const theta::Configuration & cfg = cc.get();

	BOOST_CHECKPOINT("building model");
	std::auto_ptr<Model> m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["m"]));
	BOOST_CHECKPOINT("building llvm model");
	std::auto_ptr<Model> llvm_m = PluginManager<Model>::build(Configuration(cfg, cfg.setting["llvm_m"]));
	BOOST_CHECKPOINT("llvm model built");

	// calculate some predictions:
	ParValues values;
	values.set(delta0, 0.0);
	values.set(delta1, 0.0);
	values.set(rv_delta0, 0.0);
	values.set(rv_delta1, 0.0);
	ParIds rvobs;
	rvobs.insert(rvobs0);
	rvobs.insert(rvobs1);
	DataWithUncertainties d_wu, llvm_d_wu;
	m->get_prediction(d_wu, values);
	llvm_m->get_prediction(llvm_d_wu, values);
	BOOST_CHECK(d_wu[obs0].get_uncertainty(nbins0 / 2) == 0.0);
	BOOST_CHECK(d_wu[obs1].get_uncertainty(nbins1 / 2) > 0.0);
	//double scale = d_wu[obs0].get_values().get_sum() / nbins0;
	double scale = 0;
	BOOST_CHECK(close_to(d_wu[obs0], llvm_d_wu[obs0], scale));
	BOOST_CHECK(close_to(d_wu[obs1], llvm_d_wu[obs1], scale));
	// again, to verify "overwriting" is working:
	m->get_prediction(d_wu, values);
	llvm_m->get_prediction(llvm_d_wu, values);
	BOOST_CHECK(close_to(d_wu, llvm_d_wu));


	// check normalize_to_nominal of obs0:
	const double norm_obs0 = d_wu[obs0].get_values_histogram().get_sum();
	values.set(delta0, 1.1);
	m->get_prediction(d_wu, values);
	llvm_m->get_prediction(llvm_d_wu, values);
	const double new_norm_obs0 = d_wu[obs0].get_values_histogram().get_sum();
	BOOST_CHECK(close_to_relative(norm_obs0, new_norm_obs0));
	bool res = close_to(d_wu, llvm_d_wu, scale, true);
	BOOST_CHECK(res);
	if(!res){
		dump(d_wu[obs0], "d_obs0");
		dump(llvm_d_wu[obs0], "llvm_d_obs0");
		dump(d_wu[obs1], "d_obs1");
		dump(llvm_d_wu[obs1], "llvm_d_obs1");
	}
	std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
	Random rnd(rnd_src);
	// sample some points:
	for(int i=0; i<50; ++i){
		values.set(delta0, rnd.uniform() * 20.0 - 10.0);
		values.set(delta1, rnd.uniform() * 20.0 - 10.0);
		values.set(rv_delta0, rnd.uniform() * 20.0 - 10.0);
		values.set(rv_delta1, rnd.uniform() * 20.0 - 10.0);
		m->get_prediction(d_wu, values);
		llvm_m->get_prediction(llvm_d_wu, values);
		BOOST_CHECK(close_to(d_wu, llvm_d_wu, scale, true));
	}

	// test likelihood:
	Data d_asimov;
	values.set(delta0, 0.0);
	values.set(delta1, 0.0);
	values.set(rv_delta0, 0.0);
	values.set(rv_delta1, 0.0);
	m->get_prediction(d_asimov, values);
	ParValues rvobs_values;
	rvobs_values.set(rvobs0, 0.0);
	rvobs_values.set(rvobs1, 0.0);
	d_asimov.set_rvobs_values(rvobs_values);
	BOOST_CHECKPOINT("likelihood model");
	auto_ptr<NLLikelihood> nll = m->get_nllikelihood(d_asimov);
	BOOST_CHECKPOINT("likelihood llvm model");
	auto_ptr<NLLikelihood> llvm_nll = llvm_m->get_nllikelihood(d_asimov);
	// sample a few points, also extreme ones: +-30 "sigma":
	const double range = 30.0;
	for(int i=0; i<10; ++i){
		values.set(delta0, rnd.uniform() * 2 * range - range);
		values.set(delta1, rnd.uniform() * 2 * range - range);
		values.set(rv_delta0, rnd.uniform() * 2 * range - range);
		values.set(rv_delta1, rnd.uniform() * 2 * range - range);
		double l1 = (*nll)(values);
		double l2 = (*llvm_nll)(values);
		BOOST_CHECK_EQUAL(l1,l2);
	}

}// itest
}

BOOST_AUTO_TEST_SUITE_END()

