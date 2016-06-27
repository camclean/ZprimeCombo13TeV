#include "llvm/llvm_multiply.hpp"
#include "interface/pm.hpp"
#include "interface/plugin.hpp"

#include "llvm/DerivedTypes.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"

#include "llvm/Support/IRBuilder.h"

using namespace llvm;
using namespace theta;
using namespace std;


llvm_multiply::llvm_multiply(const Configuration & cfg): literal_factor(1.0){
    size_t n = cfg.setting["factors"].size();
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    for(size_t i=0; i<n; ++i){
        Setting::Type t = cfg.setting["factors"][i].get_type();
        if(t==Setting::TypeFloat){
            literal_factor *= static_cast<double>(cfg.setting["factors"][i]);
        }
        else if(t==Setting::TypeString){
           ParId pid = vm->get_par_id(cfg.setting["factors"][i]);
           v_pids.push_back(pid);
           par_ids.insert(pid);
        }
        else if(t==Setting::TypeGroup){
            std::auto_ptr<Function> f = PluginManager<Function>::build(Configuration(cfg, cfg.setting["factors"][i]));
            const ParIds & f_p = f->get_parameters();
            par_ids.insert_all(f_p);
            functions.push_back(f);
        }
        else{
           std::stringstream ss;
           ss << "Invalid config setting type in setting 'factors' at index " << i;
           throw ConfigurationException(ss.str());
        }
    }
}

double llvm_multiply::operator()(const ParValues & v) const{
    double result = literal_factor;
    for(size_t i=0; i<v_pids.size(); ++i){
        result *= v.get_unchecked(v_pids[i]);
    }
    for(size_t i=0; i<functions.size(); ++i){
        result *= functions[i](v);
    }
    return result;
}

llvm::Function * llvm_multiply::llvm_codegen(llvm_module & mod, const std::string & prefix) const{
    FunctionType * FT = get_ft_function_evaluate(mod);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_evaluate", mod.module);
    LLVMContext & context = mod.module->getContext();
    IRBuilder<> Builder(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    Type * double_t = Type::getDoubleTy(context);
    Type * i32_t = Type::getInt32Ty(context);
    Value * val0 = Builder.CreateFAdd(ConstantFP::get(double_t, literal_factor), ConstantFP::get(double_t, 0.0));
    Value * last_value = val0;
    Value * p_par_values = F->arg_begin();
    for(size_t i=0; i<v_pids.size(); ++i){
        unsigned int pindex = mod.get_index(v_pids[i]);
        Value * p_parameter_value = Builder.CreateGEP(p_par_values, ConstantInt::get(i32_t, pindex));
        Value * parameter_value = Builder.CreateLoad(p_parameter_value);
        last_value = Builder.CreateFMul(last_value, parameter_value);
    }
    for(size_t i=0; i<functions.size(); ++i){
        stringstream ss_prefix;
        ss_prefix << prefix << "_f" << i;
        llvm::Function * f = create_llvm_function(&functions[i], mod, ss_prefix.str());
        Value * f_factor = Builder.CreateCall(f, p_par_values);
        last_value = Builder.CreateFMul(last_value, f_factor);
    }
    Builder.CreateRet(last_value);
    return F;
}

REGISTER_PLUGIN(llvm_multiply)

