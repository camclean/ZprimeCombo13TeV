#include "llvm/llvm_exp_function.hpp"
#include "interface/plugin.hpp"

#include "llvm/DerivedTypes.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"

#include "llvm/Support/IRBuilder.h"

using namespace llvm;
using namespace std;
using namespace theta;

llvm_exp_function::llvm_exp_function(const theta::Configuration & cfg){
    ParValues val_lambdas_plus, val_lambdas_minus;
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    if(cfg.setting.exists("parameters")){
        for(size_t i=0; i<cfg.setting["parameters"].size(); ++i){
            ParId pid = vm->get_par_id(cfg.setting["parameters"][i]);
            par_ids.insert(pid);
            val_lambdas_plus.set(pid, cfg.setting["lambdas_plus"][i]);
            val_lambdas_minus.set(pid, cfg.setting["lambdas_minus"][i]);
        }
    }
    else{
        ParId pid = vm->get_par_id(cfg.setting["parameter"]);
        par_ids.insert(pid);
        if(cfg.setting.exists("lambda_minus")){
            val_lambdas_minus.set(pid, cfg.setting["lambda_minus"]);
            val_lambdas_plus.set(pid, cfg.setting["lambda_plus"]);
        }
        else{
            double lambda = cfg.setting["lambda"];
            val_lambdas_minus.set(pid, lambda);
            val_lambdas_plus.set(pid, lambda);
        }
    }
    //convert to the "vector" version:
    lambdas_plus.reserve(par_ids.size());
    lambdas_minus.reserve(par_ids.size());
    v_pids.reserve(par_ids.size());
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
        v_pids.push_back(*it);
        lambdas_plus.push_back(val_lambdas_plus.get(*it));
        lambdas_minus.push_back(val_lambdas_minus.get(*it));
    }
}

double llvm_exp_function::operator()(const theta::ParValues & values) const{
    double exponent_total = 0.0;
    const size_t n = v_pids.size();
    for(size_t i=0; i<n; ++i){
        double val = values.get_unchecked(v_pids[i]);
        double lambda = val < 0 ? lambdas_minus[i] : lambdas_plus[i];
        exponent_total += lambda * val;
    }
    return exp(exponent_total);
}


llvm::Function * llvm_exp_function::llvm_codegen(llvm_module & mod, const std::string & prefix) const{
    FunctionType * FT = get_ft_function_evaluate(mod);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_evaluate", mod.module);
    LLVMContext & context = mod.module->getContext();
    IRBuilder<> Builder(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    Type * double_t = Type::getDoubleTy(context);
    Type * i32_t = Type::getInt32Ty(context);
    Value * par_values = F->arg_begin();
    Value * exponent_total = Builder.CreateFAdd(ConstantFP::get(double_t, 0.0), ConstantFP::get(double_t, 0.0));
    const size_t n = v_pids.size();
    for(size_t i=0; i<n; ++i){
        size_t pindex = mod.get_index(v_pids[i]);
        Constant * index = ConstantInt::get(i32_t, pindex);
        GetElementPtrInst* p_parameter_value = GetElementPtrInst::Create(par_values, index, "", BB);
        Value * parameter_value = Builder.CreateLoad(p_parameter_value);
        if(lambdas_plus[i] == lambdas_minus[i]){
            Value * multmp = Builder.CreateFMul(parameter_value, ConstantFP::get(double_t, lambdas_plus[i]));
            exponent_total = Builder.CreateFAdd(exponent_total, multmp);
        }
        else{
            BasicBlock * BB_join = BasicBlock::Create(context, "join", F);
            
            BasicBlock * BB_gtr0 = BasicBlock::Create(context, "gtr0", F, BB_join);
            Builder.SetInsertPoint(BB_gtr0);
            Value * multmp0 = Builder.CreateFMul(parameter_value, ConstantFP::get(double_t, lambdas_plus[i]));
            Builder.CreateBr(BB_join);
            
            BasicBlock * BB_lt0 = BasicBlock::Create(context, "lt0", F, BB_join);
            Builder.SetInsertPoint(BB_lt0);
            Value * multmp1 = Builder.CreateFMul(parameter_value, ConstantFP::get(double_t, lambdas_minus[i]));
            Builder.CreateBr(BB_join);
            
            Builder.SetInsertPoint(BB);
            Value * is_gtr_eq = Builder.CreateFCmpOGE(parameter_value, ConstantFP::get(double_t, 0.0));
            Builder.CreateCondBr(is_gtr_eq, BB_gtr0, BB_lt0);
            
            Builder.SetInsertPoint(BB_join);
            // we could also move the FADD into the branch. However, here I have the opportunity
            // for a completely valid use of PHI ...
            PHINode * phi = Builder.CreatePHI(double_t, 2);
            phi->addIncoming(multmp0, BB_gtr0);
            phi->addIncoming(multmp1, BB_lt0);
            exponent_total = Builder.CreateFAdd(exponent_total, phi);
            // continue at BB = BB_join
            Builder.SetInsertPoint(BB = BB_join);
        }
    }
    llvm::Function * exp_function = mod.module->getFunction("exp_function");
    theta_assert(exp_function);
    Value * retval = Builder.CreateCall(exp_function, exponent_total);
    Builder.CreateRet(retval);
    return F;
}

REGISTER_PLUGIN(llvm_exp_function)

