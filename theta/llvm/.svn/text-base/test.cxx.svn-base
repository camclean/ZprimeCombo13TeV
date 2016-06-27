#define BOOST_TEST_MODULE llvm_testmod
#include <boost/test/included/unit_test.hpp>
#include "llvm/llvm_interface.hpp"

#include "interface/random.hpp"
#include "interface/model.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/variables-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/distribution.hpp"
#include "test/utils.hpp"

#include "plugins/core.hpp"
#include "llvm/llvm_model.hpp"

#include "llvm/Analysis/Verifier.h"
#include "llvm/Support/IRBuilder.h"

#include <stdexcept>
#include <boost/filesystem.hpp>

using namespace theta;
using namespace llvm;
using namespace std;

namespace fs = boost::filesystem;

namespace{
void verify(llvm::Module * module, bool dump = false){
	if(dump)module->dump();
	if(verifyModule(*module, llvm::PrintMessageAction)){
		std::cout << "module verification failed; module:" << std::endl;
		throw std::invalid_argument("module verify failed");
	}
}

std::string read_file(const string & fname){
    try{
        size_t s = fs::file_size(fname);
        string result;
        result.resize(s);
        ifstream in(fname.c_str());
        if(!in){
            throw string(); // will be catched below and converted to exception with proper message
        }
        in.read(&result[0], s);
        if(!in){
            throw string();
        }
        return result;
    }
    catch(...){
        throw ConfigurationException("error reading file '" + fname + "'");
    }
}
}

BOOST_AUTO_TEST_SUITE(llvm_tests)

BOOST_AUTO_TEST_CASE(basic){
	VarIdManager vm;
	ParId p0 = vm.create_par_id("p0");
	ParId p1 = vm.create_par_id("p1");
	ParId p2 = vm.create_par_id("p2");
	theta::ParIds ids;
	ids.insert(p0);
	ids.insert(p2);
	llvm_module module2(ids);

	BOOST_CHECK(module2.get_parameters()==ids);
	BOOST_CHECK(module2.get_index(p0)==0);
	BOOST_CHECK(module2.get_index(p2)==1);
	BOOST_CHECK_THROW(module2.get_index(p1), std::invalid_argument);
}

// also test global_ddata ...
BOOST_AUTO_TEST_CASE(emit_copy_ddata) {
	const size_t sizes[] = {1, 2, 7, 13, 30, 63};
	for(int i=0; i<6; ++i){
		// the src and destination:
		const size_t n = sizes[i];
		DoubleVector src(n), dest(n);
		for(size_t k=0; k<n; ++k){
			src.set(k, -k + 23.6);
		}
		ParIds ids;
		llvm_module mod(ids);
		LLVMContext & context = mod.module->getContext();
		Type * void_ = Type::getVoidTy(context);
		Type * double_ = Type::getDoubleTy(context);
		std::vector<Type*> arg_types(1, double_->getPointerTo());
		FunctionType * FT = FunctionType::get(void_, arg_types, false);
		Value * src_ = mod.add_global_ddata(src.get_data(), n, "ddata", true);
		llvm::Function * f = llvm::Function::Create(FT, llvm::Function::PrivateLinkage, "test_emit_copy_ddata", mod.module);
		llvm::Function::arg_iterator dest_ = f->arg_begin();
		BOOST_CHECKPOINT("emitting copy_ddata");
		BasicBlock * BB = BasicBlock::Create(context, "entry", f);
		BB = mod.emit_copy_ddata(BB, dest_, src_, n);
		IRBuilder<> Builder(context);
		Builder.SetInsertPoint(BB);
		Builder.CreateRetVoid();
		verify(mod.module);
		BOOST_CHECKPOINT("getting function pointer");
		void (*f_)(double*) = reinterpret_cast<void (*)(double*)>(mod.getFunctionPointer(f));
		BOOST_CHECKPOINT("calling function");
		f_(dest.get_data());
		for(size_t k=0; k<n; ++k){
			BOOST_CHECK(src.get(k) == dest.get(k));
		}
	}
}

BOOST_AUTO_TEST_CASE(add_with_coeff_function){
	const size_t sizes[] = {1, 2, 7, 13, 30, 63};
	for(int i=0; i<6; ++i){
		const size_t n = sizes[i];
		DoubleVector h0(n), hadd(n), hexp(n);
		for(size_t k=0; k<n; ++k){
			h0.set(k, -k*2.6);
			hadd.set(k, k*k + 13.2);
		}
		hexp = h0;
		const double coeff = M_PI;
		hexp.add_with_coeff(coeff, hadd);
		// same within llvm:
		ParIds ids;
		llvm_module mod(ids);
		llvm::Function * awc = mod.module->getFunction("add_with_coeff");
		BOOST_REQUIRE(awc);
		typedef void (*t_add_with_coeff)(double, double *, const double *, int);
		t_add_with_coeff add_with_coeff = reinterpret_cast<t_add_with_coeff>(mod.getFunctionPointer(awc));
		add_with_coeff(coeff, h0.get_data(), hadd.get_data(), n);
		for(size_t k=0; k<n; ++k){
			BOOST_CHECK(h0.get(k) == hexp.get(k));
		}
	}
}

BOOST_AUTO_TEST_CASE(emit_add_with_coeff){
	const size_t sizes[] = {1, 2, 7, 13, 30, 63};
	for(int i=0; i<6; ++i){
		const size_t n = sizes[i];
		DoubleVector hleft(n), hright(n), hexp(n);
		for(size_t k=0; k<n; ++k){
			hleft.set(k, -k*2.6);
			hright.set(k, k*k + 13.2);
		}
		hexp = hleft;
		const double coeff = M_PI;
		hexp.add_with_coeff(coeff, hright);
		ParIds ids;
		llvm_module mod(ids);
		LLVMContext & context = mod.module->getContext();
		// make function void awc(double coeff, double * data_lhs, double * data_rhs) which calcultes data_lhs += coeff * data_rhs
		Type * void_ = Type::getVoidTy(context);
		Type * double_ = Type::getDoubleTy(context);
		std::vector<Type*> arg_types(3);
		arg_types[0] = double_;
		arg_types[2] = arg_types[1] = double_->getPointerTo();
		FunctionType * FT = FunctionType::get(void_, arg_types, false);
		llvm::Function * f = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "awc", mod.module);
		llvm::Function::arg_iterator it = f->arg_begin();
		llvm::Value * lcoeff = it++;
		llvm::Value * ldata_lhs = it++;
		llvm::Value * ldata_rhs = it++;
		BasicBlock * BB = BasicBlock::Create(context, "entry", f);
		BB = mod.emit_add_with_coeff(BB, lcoeff, ldata_lhs, ldata_rhs, n);
		IRBuilder<> Builder(context);
		Builder.SetInsertPoint(BB);
		Builder.CreateRetVoid();
		llvm::Function * llvm_awc = mod.module->getFunction("awc");
		BOOST_REQUIRE(llvm_awc);
		typedef void (*t_awc)(double, double *, double *);
		t_awc awc = reinterpret_cast<t_awc>(mod.getFunctionPointer(llvm_awc));
		awc(coeff, hleft.get_data(), hright.get_data());
		for(size_t k=0; k<n; ++k){
			BOOST_CHECK(hleft.get(k) == hexp.get(k));
		}
	}
}


BOOST_AUTO_TEST_CASE(emit_add_with_coeff_const_rhs){
	const size_t sizes[] = {1, 2, 7, 13, 30, 63};
	for(int i=0; i<6; ++i){
		const size_t n = sizes[i];
		DoubleVector hleft(n), hright(n), hexp(n);
		for(size_t k=0; k<n; ++k){
			hleft.set(k, -k*2.6);
			hright.set(k, k*k + 13.2);
		}
		hexp = hleft;
		const double coeff = M_PI;
		hexp.add_with_coeff(coeff, hright);
		ParIds ids;
		llvm_module mod(ids);
		LLVMContext & context = mod.module->getContext();
		// make function void awc(double coeff, double * data_lhs) which calculates data_lhs += coeff * data_rhs
		Type * void_ = Type::getVoidTy(context);
		Type * double_ = Type::getDoubleTy(context);
		std::vector<Type*> arg_types(2);
		arg_types[0] = double_;
		arg_types[1] = double_->getPointerTo();
		FunctionType * FT = FunctionType::get(void_, arg_types, false);
		llvm::Function * f = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "awc", mod.module);
		llvm::Function::arg_iterator it = f->arg_begin();
		llvm::Value * lcoeff = it++;
		llvm::Value * ldata_lhs = it++;
		BasicBlock * BB = BasicBlock::Create(context, "entry", f);
		BB = mod.emit_add_with_coeff(BB, lcoeff, ldata_lhs, hright);
		IRBuilder<> Builder(context);
		Builder.SetInsertPoint(BB);
		Builder.CreateRetVoid();
		llvm::Function * llvm_awc = mod.module->getFunction("awc");
		BOOST_REQUIRE(llvm_awc);
		typedef void (*t_awc)(double, double *);
		t_awc awc = reinterpret_cast<t_awc>(mod.getFunctionPointer(llvm_awc));
		awc(coeff, hleft.get_data());
		for(size_t k=0; k<n; ++k){
			BOOST_CHECK(hleft.get(k) == hexp.get(k));
		}
	}
}

BOOST_AUTO_TEST_CASE(emit_multiply){
	const size_t sizes[] = {1, 2, 7, 13, 30, 63};
	for(int i=0; i<6; ++i){
		const size_t n = sizes[i];
		DoubleVector h0(n), hexp(n);
		for(size_t k=0; k<n; ++k){
			h0.set(k, -k*2.6);
		}
		hexp = h0;
		const double coeff = M_PI;
		hexp *= coeff;
		ParIds ids;
		llvm_module mod(ids);
		// create a function "void multiply(double * , double)":
		LLVMContext & context = mod.module->getContext();
		Type * void_ = Type::getVoidTy(context);
		Type * double_ = Type::getDoubleTy(context);
		std::vector<Type*> arg_types(2);
		arg_types[0] = double_->getPointerTo();
		arg_types[1] = double_;
		FunctionType * FT = FunctionType::get(void_, arg_types, false);
		llvm::Function * f = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "multiply", mod.module);
		llvm::Function::arg_iterator it = f->arg_begin();
		Value * vector = it++;
		Value * coeffl = it;
		BasicBlock * BB = BasicBlock::Create(context, "entry", f);
		BB = mod.emit_multiply(BB, coeffl, vector, n);
		IRBuilder<> Builder(context);
		Builder.SetInsertPoint(BB);
		Builder.CreateRetVoid();
		//verify(mod.module, true);
		typedef void (*t_multiply)(double, double *);
		t_multiply multiply = reinterpret_cast<t_multiply>(mod.getFunctionPointer(f));
		multiply(coeff, h0.get_data());
		for(size_t k=0; k<n; ++k){
			BOOST_CHECK(h0.get(k) == hexp.get(k));
		}
	}
}

BOOST_AUTO_TEST_CASE(emit_normalize_histo){
	const size_t sizes[] = {1, 2, 7, 13, 30, 63};
	std::auto_ptr<RandomSource> src(new RandomSourceTaus());
	Random rnd(src);
	for(int in=0; in<12; ++in){
		int i = in % 6;
		const size_t n = sizes[i];
		DoubleVector h0(n), hexp(n);
		for(size_t k=0; k<n; ++k){
			h0.set(k, rnd.uniform() - 0.5);
		}
		hexp = h0;
		double sum = 0.0;
		for(size_t k=0; k<n; ++k){
			if(hexp.get(k) < 0.0) hexp.set(k, 0.0);
			else sum += hexp.get(k);
		}
		boost::optional<double> norm_to;
		double exp_coeff = 1.0;
		if(in / 6 > 0){
			norm_to = 2.7;
			if(sum > 0.0){
				hexp *= *norm_to / sum;
				exp_coeff = *norm_to / sum;
			}
			else{
			    //cout << "0.0" << endl;
			    exp_coeff = 0.0;
			}
		}
		ParIds ids;
		llvm_module mod(ids);
		// create a function "void normalize(double *)":
		LLVMContext & context = mod.module->getContext();
		Type * void_ = Type::getVoidTy(context);
		Type * double_ = Type::getDoubleTy(context);
		std::vector<Type*> arg_types(1, double_->getPointerTo());
		FunctionType * FT = FunctionType::get(void_, arg_types, false);
		llvm::Function * f = llvm::Function::Create(FT, llvm::Function::PrivateLinkage, "normalize", mod.module);
		llvm::Function::arg_iterator it = f->arg_begin();
		BasicBlock * BB = BasicBlock::Create(context, "entry", f);
		Value * v = 0;
		BB = mod.emit_normalize_histo(BB, it, n, norm_to, v);
		IRBuilder<> Builder(context);
		Builder.SetInsertPoint(BB);
		if(norm_to){
		    theta_assert(v);
		    Builder.CreateRet(v);
		}
		else{
		    theta_assert(v == 0);
		    Builder.CreateRet(ConstantFP::get(double_, 1.0));
		}
		//verify(mod.module, true);
		typedef double (*t_normalize)(double *);
		t_normalize normalize = reinterpret_cast<t_normalize>(mod.getFunctionPointer(f));
		double coeff = normalize(h0.get_data());
		BOOST_CHECK_EQUAL(coeff, exp_coeff);
		//printf("%d\n", in);
		for(size_t k=0; k<n; ++k){
			//printf("%d: %8.3g %8.3g\n", int(k), h0.get(k), hexp.get(k));
			BOOST_CHECK_EQUAL(h0.get(k), hexp.get(k));
		}
	}
}

// this is to force linking of core-plugins.so
void f(){
 new fixed_gauss();
}


BOOST_AUTO_TEST_CASE(htt_model){
    for(int i=0; i<2; ++i){
    string fname = i==0?"htt-nobbunc.cfg":"htt-bbunc.cfg";
    cout << fname << endl;
    string cfg_string = read_file(fname);
    boost::shared_ptr<SettingUsageRecorder> rec(new SettingUsageRecorder());
    Setting root = LibconfigSetting::parse(cfg_string, rec);
    Configuration config(root);
    config.pm->set("default", boost::shared_ptr<VarIdManager>(new VarIdManager));
    config.pm->set("default", boost::shared_ptr<ProductsSink>(new NullProductsSink));
    config.pm->set("runid", boost::shared_ptr<int>(new int(1)));
    apply_vm_settings(config);

    std::auto_ptr<Model> m(new default_model(Configuration(config, root["main"]["model"])));
    std::auto_ptr<Model> ll_m(new llvm_model(Configuration(config, root["main"]["model"]), true));

    std::auto_ptr<DataSource> ds = PluginManager<DataSource>::build(Configuration(config, root["main"]["data_source"]));
    Data d;
    ds->fill(d);

    std::auto_ptr<NLLikelihood> nll = m->get_nllikelihood(d);
    std::auto_ptr<NLLikelihood> llvm_nll = ll_m->get_nllikelihood(d);

    // sample a few points from the prior distribution and check the value of the likelihood there:
    std::auto_ptr<Distribution> dist = PluginManager<Distribution>::build(Configuration(config, root["main"]["model"]["parameter-distribution"]));

    DataWithUncertainties d_pred, llvm_pred;
    ParValues v;
    dist->mode(v);

    m->get_prediction(d_pred, v);
    ll_m->get_prediction(llvm_pred, v);

    ObsIds obs = d_pred.get_observables();
    for(ObsIds::const_iterator it=obs.begin(); it!=obs.end(); ++it){
        Histogram1DWithUncertainties h_pred = d_pred[*it];
        Histogram1DWithUncertainties hwu_pred = llvm_pred[*it];
        for(size_t i=0; i<h_pred.get_nbins(); ++i){
            double v1 = h_pred.get(i);
            double v2 = hwu_pred.get(i);
            BOOST_CHECK(close_to(v1, v2, 1e-2));
            v1 = h_pred.get_uncertainty2(i);
            v2 = hwu_pred.get_uncertainty2(i);
            BOOST_CHECK(close_to(v1, v2, 1e-2));
        }
    }

    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    for(int i=0; i<10; ++i){
        dist->sample(v, rnd);
        double val1 = (*nll)(v);
        double val2 = (*llvm_nll)(v);
        BOOST_CHECK(close_to_relative(val1, val2));
        //BOOST_CHECK(close_to_relative(val1, val2));
        //cout << val1 << " " << val2 << " " << (val1 - val2) << endl;
    }
    }
}

BOOST_AUTO_TEST_SUITE_END()
