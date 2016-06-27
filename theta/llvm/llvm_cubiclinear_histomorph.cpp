#include "llvm/llvm_cubiclinear_histomorph.hpp"
#include "llvm/llvm_interface.hpp"
#include "interface/plugin.hpp"

#include "llvm/Module.h"
#include "llvm/Function.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Analysis/Verifier.h"

#include "boost/optional.hpp"

#include <iostream>

using namespace theta;
using namespace llvm;
using namespace std;

void llvm_cubiclinear_histomorph::apply_functor(const theta::functor<theta::Histogram1D> & f, const theta::ParValues & values) const{
	h = h0;
	const size_t n_sys = hplus_diff.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]) * parameter_factors[isys];
        if(delta==0.0) continue;
        //linear extrapolation beyond 1 sigma:
        if(fabs(delta) > 1){
            const Histogram1D & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
            h.add_with_coeff(fabs(delta), t_sys);
        }
        else{
            //cubic interpolation:
            const double d2 = delta * delta;
            const double d3 = d2 * fabs(delta);
            h.add_with_coeff(0.5*delta, diff[isys]);
            h.add_with_coeff(d2 - 0.5 * d3, sum[isys]);
        }
    }
    double h_sum = 0.0;
    for(size_t i=0; i < h.get_nbins(); ++i){
        double val = h.get(i);
        if(val < 0.0){
            h.set(i, 0.0);
        }
        else{
            h_sum += val;
        }
    }
    if(normalize_to_nominal && h_sum > 0.0){
       h *= h0_sum / h_sum;
    }
    f(h);
}


llvm::Function * llvm_cubiclinear_histomorph::llvm_codegen(llvm_module & mod, const std::string & prefix) const{
	// 1. create global data; most of it const, but one buffer "h" for writing the
	// intermediate result (required for 0-truncation)
	const size_t nbins = h0.get_nbins();
	const size_t nsys = hplus_diff.size();
	std::vector<GlobalVariable*> llvm_hplus_diff;
	std::vector<GlobalVariable*> llvm_hminus_diff;
	std::vector<GlobalVariable*> llvm_diff;
	std::vector<GlobalVariable*> llvm_sum;
	theta_assert(hminus_diff.size() == nsys);
	theta_assert(diff.size() == nsys);
	theta_assert(sum.size() == nsys);
	GlobalVariable* llvm_h0 = mod.add_global_ddata(h0.get_data(), nbins, prefix + "h0", true);
	GlobalVariable* llvm_h = mod.add_global_ddata(h.get_data(), nbins, prefix + "h", false);
	for(size_t i=0; i<nsys; ++i){
		std::stringstream ss;
		ss << i;
		llvm_hplus_diff.push_back(mod.add_global_ddata(hplus_diff[i].get_data(), nbins, prefix + "hplus_diff" + ss.str(), true));
		llvm_hminus_diff.push_back(mod.add_global_ddata(hminus_diff[i].get_data(), nbins, prefix + "hminus_diff" + ss.str(), true));
		llvm_diff.push_back(mod.add_global_ddata(diff[i].get_data(), nbins, prefix + "diff" + ss.str(), true));
		llvm_sum.push_back(mod.add_global_ddata(sum[i].get_data(), nbins, prefix + "sum" + ss.str(), true));
	}
	// 2. create function and read arguments:
	FunctionType * FT = get_ft_hf_add_with_coeff(mod);
	llvm::Function * F = llvm::Function::Create(FT, llvm::Function::PrivateLinkage, prefix + "_cl_evaluate", mod.module);
	llvm::Function::arg_iterator it = F->arg_begin();
	Argument * coeff = it++;
	Argument * par_values = it++;
	Argument * data = it++;
	data->addAttr(Attribute::constructAlignmentFromInt(16));
	theta_assert(it == F->arg_end());

	// 3. build function body
	LLVMContext & context = mod.module->getContext();
	IRBuilder<> Builder(context);
	BasicBlock * BB = BasicBlock::Create(context, "entry", F);
	Builder.SetInsertPoint(BB);
	Type * double_t = Type::getDoubleTy(context);
	Type * i32_t = Type::getInt32Ty(context);
	llvm::Function * add_with_coeff = mod.module->getFunction("add_with_coeff");
	Value * llvm_h_p = Builder.CreateBitCast(llvm_h, double_t->getPointerTo(), "h_pointer");

	// h = h0:
	BB = mod.emit_copy_ddata(BB, llvm_h, llvm_h0, nbins);
	Builder.SetInsertPoint(BB);
	for (size_t isys = 0; isys < nsys; isys++) {
		Constant * index = ConstantInt::get(i32_t, mod.get_index(vid[isys]));
		Value * p_parameter_value = Builder.CreateGEP(par_values, index, "par_values");
		Value * delta = Builder.CreateLoad(p_parameter_value, "delta");
		if(parameter_factors[isys]!=1.0){
			delta = Builder.CreateFMul(delta, ConstantFP::get(double_t, parameter_factors[isys]), "delta_pf");
		}
		// abs_delta = delta > 0.0 ? delta : -delta   :
		Value * abs_delta;
		Value * delta_gtr0;
		{
			BasicBlock * BB_join = BasicBlock::Create(context, "deltajoin", F);
			BasicBlock * BB_lt0 = BasicBlock::Create(context, "deltalt0", F, BB_join);
			delta_gtr0 = Builder.CreateFCmpOGE(delta, ConstantFP::get(double_t, 0.0), "delta_gtr0");
			Builder.CreateCondBr(delta_gtr0, BB_join, BB_lt0);

			Builder.SetInsertPoint(BB_lt0);
			Value * minus_delta = Builder.CreateFNeg(delta, "minus_delta");
			Builder.CreateBr(BB_join);

			Builder.SetInsertPoint(BB_join);
			PHINode * phi = Builder.CreatePHI(double_t, 2, "abs_delta");
			phi->addIncoming(delta, BB);
			phi->addIncoming(minus_delta, BB_lt0);
			abs_delta = phi;

			Builder.SetInsertPoint(BB = BB_join);
		}
		{
			// if (abs_delta > 1) .. else blocks:
			BasicBlock * BB_join = BasicBlock::Create(context, "ad1j", F);
			BasicBlock * BB_gt1 = BasicBlock::Create(context, "adgt1", F, BB_join);
			BasicBlock * BB_lt1 = BasicBlock::Create(context, "adlt1", F, BB_join);
			Value * dgt1 = Builder.CreateFCmpOGE(abs_delta, ConstantFP::get(double_t, 1.0), "adgtr1");
			Builder.CreateCondBr(dgt1, BB_gt1, BB_lt1);

			// abs_delta >= 1 branch:
			Builder.SetInsertPoint(BB_gt1);
			{
				// t_sys = delta > 0 ? hplus_diff : hminus_diff
				BasicBlock * BB_join0 = BasicBlock::Create(context, "d0j", F);
				BasicBlock * BB_gt0 = BasicBlock::Create(context, "dgt0", F, BB_join0);
				BasicBlock * BB_lt0 = BasicBlock::Create(context, "dlt0", F, BB_join0);
				Builder.CreateCondBr(delta_gtr0, BB_gt0, BB_lt0);

				// delta > 0:
				Builder.SetInsertPoint(BB_gt0);
				Value * t_sys_gt0 = Builder.CreateBitCast(llvm_hplus_diff[isys], double_t->getPointerTo(), "t_sys_gt0");
				Builder.CreateBr(BB_join0);

				// delta < 0:
				Builder.SetInsertPoint(BB_lt0);
				Value * t_sys_lt0 = Builder.CreateBitCast(llvm_hminus_diff[isys], double_t->getPointerTo(), "t_sys_lt0");
				Builder.CreateBr(BB_join0);

				Builder.SetInsertPoint(BB_join0);
				PHINode *t_sys = Builder.CreatePHI(double_t->getPointerTo(), 2, "t_sys");
				t_sys->addIncoming(t_sys_gt0, BB_gt0);
				t_sys->addIncoming(t_sys_lt0, BB_lt0);

				// h += abs_delta * t_sys:
				Builder.CreateCall4(add_with_coeff, abs_delta, llvm_h_p, t_sys, ConstantInt::get(i32_t, nbins));
				Builder.CreateBr(BB_join); // closes BB_gt1
			}
			Builder.SetInsertPoint(BB_lt1);
			{
				//abs_delta < 1 branch:
				//  h.add_with_coeff(0.5*delta, diff[isys]);
				//  const double d2 = delta * delta;
				//  const double d3 = d2 * fabs(delta);
				//  h.add_with_coeff(d2 - 0.5 * d3, sum[isys]);
				Value * delta_o_2 = Builder.CreateFMul(delta, ConstantFP::get(double_t, 0.5));
				Value * diff = Builder.CreateBitCast(llvm_diff[isys], double_t->getPointerTo(), "diff");
				Builder.CreateCall4(add_with_coeff, delta_o_2, llvm_h_p, diff, ConstantInt::get(i32_t, nbins));
				Value * delta2 = Builder.CreateFMul(delta, delta);
				Value * delta3 = Builder.CreateFMul(delta2, abs_delta);
				Value * delta3_o_2 = Builder.CreateFMul(delta3, ConstantFP::get(double_t, 0.5));
				Value * coeff = Builder.CreateFSub(delta2, delta3_o_2);
				Value * sum = Builder.CreateBitCast(llvm_sum[isys], double_t->getPointerTo(), "sum");
				Builder.CreateCall4(add_with_coeff, coeff, llvm_h_p, sum, ConstantInt::get(i32_t, nbins));
				Builder.CreateBr(BB_join); //closes BB_lt1
				Builder.SetInsertPoint(BB = BB_join);
			}
			Builder.SetInsertPoint(BB = BB_join);
		}
	}

	boost::optional<double> normalize_to;
	if(normalize_to_nominal){
		normalize_to = h0_sum;
	}
	// truncate llvm_h_p at zero:
	llvm::Value * coeff_to_norm = 0;
	BB = mod.emit_normalize_histo(BB, llvm_h_p, nbins, normalize_to, coeff_to_norm);
	if(normalize_to_nominal){
	    theta_assert(coeff_to_norm!=0);
	}
	Builder.SetInsertPoint(BB);
	BB = mod.emit_add_with_coeff(BB, coeff, data, llvm_h_p, nbins);
	Builder.SetInsertPoint(BB);
	if(normalize_to_nominal){
	    Builder.CreateRet(coeff_to_norm);
	}
	else{
	    Builder.CreateRet(ConstantFP::get(double_t, 1.0));
	}
	return F;
}

void llvm_cubiclinear_histomorph::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h0.get_nbins();
    xmin = h0.get_xmin();
    xmax = h0.get_xmax();
}

theta::Histogram1D llvm_cubiclinear_histomorph::get_uncertainty2_histogram() const{
	return h0_uncertainty2;
}


llvm_cubiclinear_histomorph::llvm_cubiclinear_histomorph(const Configuration & ctx): normalize_to_nominal(false) {
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    //build nominal histogram:
    Histogram1DWithUncertainties h0_wu = get_constant_histogram(Configuration(ctx, ctx.setting["nominal-histogram"]));
    h0 = h0_wu.get_values_histogram();
    h0_uncertainty2 = h0; // to set dimensions
    for(size_t i=0; i<h0.get_nbins(); ++i){
    	h0_uncertainty2.set(i, h0_wu.get_uncertainty2(i));
    }
    if(ctx.setting.exists("normalize_to_nominal")){
        normalize_to_nominal = ctx.setting["normalize_to_nominal"];
    }
    if(normalize_to_nominal){
    	//throw ConfigurationException("normalize_to_nominal not yet implemented!");
    }
    Setting psetting = ctx.setting["parameters"];
    size_t n = psetting.size();
    parameter_factors.resize(n, 1.0);
    bool have_parameter_factors = ctx.setting.exists("parameter_factors");
    for(size_t i=0; i<n; i++){
        std::string par_name = psetting[i];
        ParId pid = vm->get_par_id(par_name);
        par_ids.insert(pid);
        vid.push_back(pid);
        std::string setting_name;
        //plus:
        setting_name = par_name + "-plus-histogram";
        hplus_diff.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[setting_name])).get_values_histogram());
        hplus_diff.back().check_compatibility(h0);
        hplus_diff.back().add_with_coeff(-1.0, h0);
        //minus:
        setting_name = par_name + "-minus-histogram";
        hminus_diff.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[setting_name])).get_values_histogram());
        hminus_diff.back().check_compatibility(h0);
        hminus_diff.back().add_with_coeff(-1.0, h0);

        sum.push_back(hplus_diff.back());
        sum.back() += hminus_diff.back();
        diff.push_back(hplus_diff.back());
        diff.back().add_with_coeff(-1, hminus_diff.back());

        if(have_parameter_factors){
            parameter_factors[i] = ctx.setting["parameter_factors"][i];
        }
    }
    h0_sum = h0.get_sum();
    h = h0;
}

REGISTER_PLUGIN(llvm_cubiclinear_histomorph)
