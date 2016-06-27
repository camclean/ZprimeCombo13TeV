#include "llvm/llvm_interface.hpp"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Target/TargetData.h"

#include "llvm/Analysis/Verifier.h"
#include "llvm/Attributes.h"

#include "llvm/PassManager.h"
#include "llvm/DefaultPasses.h"

#include "llvm/Transforms/Scalar.h"
#include "llvm/Transforms/IPO.h"
#include "llvm/Transforms/Vectorize.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"

#include "interface/plugin.hpp"
#include <dlfcn.h>

#include <iostream>

using namespace llvm;
using namespace theta;
using namespace std;


namespace {


struct llvm_initializer{
	llvm_initializer(){
		InitializeNativeTarget();
	}
} init;

std::map<theta::ParId, size_t> create_pid_to_index(const ParIds & pids){
	std::map<theta::ParId, size_t> result;
	size_t i=0;
	for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
		result[*it] = i++;
	}
	return result;
}

}

void dump_info(double d, int i){
	std::cout << "dump_info: " << i << " " << d << std::endl;
}

/*
void llvm_module::emit_dump(llvm::Value * d_, llvm::Value * int_, IRBuilder<> & Builder){
    Builder.CreateCall2(f_dump_info, d_, int_);
}*/


// create the module function
// inline void add_with_coeff(double coeff, double * data_out, const double * data_in, int n);
// which calculates data_out[i] += coeff * data_in[i] for i = 0..n, where data and rhs_data are 16-byte aligned double-vectors of
// length n  (or n + (n%2)).
void llvm_module::emit_add_with_coeff_function(){
    LLVMContext & context = module->getContext();
    std::vector<Type*> arg_types(4);
    arg_types[0] = double_;
    arg_types[1] = arg_types[2] = double_->getPointerTo();
    arg_types[3] = i32_;
    FunctionType * FT = FunctionType::get(void_, arg_types, false);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "add_with_coeff", module);
    F->addFnAttr(Attribute::InlineHint);
    llvm::Function::arg_iterator iter = F->arg_begin();
    Argument * coeff = iter++;
    Argument * data_ = iter++;
    Argument * rhs_data_ = iter++;
    Argument * n_ = iter;
    data_->addAttr(Attribute::constructAlignmentFromInt(16));
    rhs_data_->addAttr(Attribute::constructAlignmentFromInt(16));
    
    IRBuilder<> Builder(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    std::vector<Constant *> zeros(2);
    zeros[0] = zeros[1] = ConstantFP::get(double_, 0.0);
    Value * coeff2tmp = Builder.CreateInsertElement(ConstantVector::get(zeros), coeff, ConstantInt::get(i32_, 0));
    Value * coeff2 = Builder.CreateInsertElement(coeff2tmp, coeff, ConstantInt::get(i32_, 1), "coeff2");
    // n_half = (n_ + 1)/2:
    Value * n1 = Builder.CreateAdd(n_, ConstantInt::get(i32_, 1), "n1");
    Value * n_half = Builder.CreateUDiv(n1, ConstantInt::get(i32_, 2), "n_half");
    // create data and rhs_data to point to data_ / rhs_data_ but of type pointer to double2_t:
    Value * data = Builder.CreateBitCast(data_, double2_->getPointerTo(), "data");
    Value * rhs_data = Builder.CreateBitCast(rhs_data_, double2_->getPointerTo(), "rhs_data");
    // make a loop i = 0..n_half
    BasicBlock * BB_loop_header = BasicBlock::Create(context, "loop_header", F); // here compare i with n and jump to body or exit
    BasicBlock * BB_loop_body = BasicBlock::Create(context, "loop_body", F);
    BasicBlock * BB_loop_exit = BasicBlock::Create(context, "loop_exit", F);
    Builder.CreateBr(BB_loop_header);
    Builder.SetInsertPoint(BB_loop_header);
    PHINode * i_phi = Builder.CreatePHI(i32_, 2, "i");
    theta_assert(i_phi != 0);
    i_phi->addIncoming(ConstantInt::get(i32_, 0), BB);
    // add the i_phi incoming value for next loop value later
    Value * i_eq_n = Builder.CreateICmpEQ(i_phi, n_half);
    Builder.CreateCondBr(i_eq_n, BB_loop_exit, BB_loop_body);
    Builder.SetInsertPoint(BB_loop_body);
    // get data[i] and rhs_data[i] as double2_t
    Value * p_data_i = GetElementPtrInst::Create(data, i_phi, "p_data_i", BB_loop_body);
    Value * p_rhs_data_i = GetElementPtrInst::Create(rhs_data, i_phi, "p_rhs_data_i", BB_loop_body);
    Value * rhs_data_i = Builder.CreateLoad(p_rhs_data_i, "rhs_data_i");
    Value * multmp = Builder.CreateFMul(rhs_data_i, coeff2);
    Value * data_i = Builder.CreateLoad(p_data_i, "data_i");
    Value * res = Builder.CreateFAdd(multmp, data_i, "res");
    Builder.CreateStore(res, p_data_i);
    Value * next_i = Builder.CreateAdd(i_phi, ConstantInt::get(i32_, 1), "next_i");
    i_phi->addIncoming(next_i, BB_loop_body);
    Builder.CreateBr(BB_loop_header);
    Builder.SetInsertPoint(BB_loop_exit);
    Builder.CreateRetVoid();
    //llvm_verify(F, "add_with_coeff");
}

llvm::BasicBlock * llvm_module::emit_add_with_coeff(llvm::BasicBlock * BB_in, llvm::Value * coeff, llvm::Value * data_out,
		llvm::Value * data_in, size_t n){
	LLVMContext & context = module->getContext();
	IRBuilder<> Builder(context);
	Builder.SetInsertPoint(BB_in);
	llvm::Function * add_with_coeff = module->getFunction("add_with_coeff");
	Builder.CreateCall4(add_with_coeff, coeff, data_out, data_in, ConstantInt::get(i32_, n));
	return BB_in;
}

BasicBlock * llvm_module::emit_add_with_coeff(BasicBlock * BB_in, Value * coeff, Value * data_out, const DoubleVector & h){
	LLVMContext & context = module->getContext();
	IRBuilder<> Builder(context);
	Builder.SetInsertPoint(BB_in);
	const size_t n = h.size();
	theta_assert(n>0);
	if(n==1){
		Value * data_o = Builder.CreateBitCast(data_out, double_->getPointerTo(), "data_out");
		Value * p_data0 = Builder.CreateGEP(data_o, ConstantInt::get(i32_, 0));
		Value * data0 = Builder.CreateLoad(p_data0);
		Value * coeff_times_h = Builder.CreateFMul(coeff, ConstantFP::get(double_, h.get(0)));
		Value * res = Builder.CreateFAdd(data0, coeff_times_h);
		Builder.CreateStore(res, p_data0);
	}
	else if(n<=8){
		std::vector<Constant *> zeros(2);
		zeros[0] = zeros[1] = ConstantFP::get(double_, 0.0);
		Value * coeff2tmp = Builder.CreateInsertElement(ConstantVector::get(zeros), coeff, ConstantInt::get(i32_, 0));
		Value * coeff2 = Builder.CreateInsertElement(coeff2tmp, coeff, ConstantInt::get(i32_, 1), "coeff2");
		Value * data_o = Builder.CreateBitCast(data_out, double2_->getPointerTo(), "data_out");
		const size_t imax = (n + n%2) / 2;
		for(size_t i=0; i < imax; ++i){
			Value * p_data_i = Builder.CreateGEP(data_o, ConstantInt::get(i32_, i));
			Value * data_i = Builder.CreateLoad(p_data_i);
			Value * h2tmp = Builder.CreateInsertElement(ConstantVector::get(zeros), ConstantFP::get(double_, h.get(2*i)), ConstantInt::get(i32_, 0));
			Value * h2 = h2tmp;
			if(2*i+1 < n){
				h2 = Builder.CreateInsertElement(h2tmp, ConstantFP::get(double_, h.get(2*i + 1)), ConstantInt::get(i32_, 1), "coeff2");
			}
			Value * coeff_times_h = Builder.CreateFMul(coeff2, h2);
			Value * res = Builder.CreateFAdd(data_i, coeff_times_h);
			Builder.CreateStore(res, p_data_i);
		}
	}
	else{
		llvm::Function * add_with_coeff = module->getFunction("add_with_coeff");
		GlobalVariable * gv_data = add_global_ddata(h.get_data(), h.size(), "ewc_data", true);
		Value * gv_data_conv = Builder.CreateBitCast(gv_data, double_->getPointerTo());
		Builder.CreateCall4(add_with_coeff, coeff, data_out, gv_data_conv, ConstantInt::get(i32_, h.size()));
	}
	return BB_in;
}

BasicBlock * llvm_module::emit_multiply(BasicBlock * BB_in, Value * coeff, Value * histo_, size_t n){
	LLVMContext & context = module->getContext();
	IRBuilder<> Builder(context);
	llvm::Function * f = BB_in->getParent();
	Builder.SetInsertPoint(BB_in);
	Value * histo = Builder.CreateBitCast(histo_, double2_->getPointerTo(), "histo");

	std::vector<Constant *> zeros(2);
	zeros[0] = zeros[1] = ConstantFP::get(double_, 0.0);
	Value * coeff2tmp = Builder.CreateInsertElement(ConstantVector::get(zeros), coeff, ConstantInt::get(i32_, 0));
	Value * coeff2 = Builder.CreateInsertElement(coeff2tmp, coeff, ConstantInt::get(i32_, 1), "coeff2");

	BasicBlock * BB_loop_exit = BasicBlock::Create(context, "mul_loop_exit", f);
	BasicBlock * BB_loop_header = BasicBlock::Create(context, "mul_loop_header", f, BB_loop_exit);
	BasicBlock * BB_loop_body = BasicBlock::Create(context, "mul_loop_body", f, BB_loop_exit);

	Builder.CreateBr(BB_loop_header);
	Builder.SetInsertPoint(BB_loop_header);
	PHINode * i_phi = Builder.CreatePHI(i32_, 2, "i");
	i_phi->addIncoming(ConstantInt::get(i32_, 0), BB_in);
	Value * i_eq_n = Builder.CreateICmpEQ(i_phi, ConstantInt::get(i32_, (n + n%2) / 2));
	Builder.CreateCondBr(i_eq_n, BB_loop_exit, BB_loop_body);

	Builder.SetInsertPoint(BB_loop_body);
	Value * p_histo_i = Builder.CreateGEP(histo, i_phi, "p_src_i");
	Value * histo_i = Builder.CreateLoad(p_histo_i, "src_i");
	Value * multiplied = Builder.CreateFMul(histo_i, coeff2);
	Builder.CreateStore(multiplied, p_histo_i);
	Value * next_i = Builder.CreateAdd(i_phi, ConstantInt::get(i32_, 1), "next_i");
	i_phi->addIncoming(next_i, BB_loop_body);
	Builder.CreateBr(BB_loop_header);

	return BB_loop_exit;
}


BasicBlock * llvm_module::emit_normalize_histo(llvm::BasicBlock * BB_in, Value * histo_, size_t n, const boost::optional<double> & normalize_to, llvm::Value *& coeff){
	LLVMContext & context = module->getContext();
	IRBuilder<> Builder(context);
	llvm::Function * f = BB_in->getParent();
	BasicBlock * BB = BB_in;
	Builder.SetInsertPoint(BB);
	Value * histo = Builder.CreateBitCast(histo_, double_->getPointerTo(), "histo");
	Value * sum_p = 0;
	if(normalize_to){
		sum_p = Builder.CreateAlloca(double_, ConstantInt::get(i32_, 1), "sum_p");
		Builder.CreateStore(ConstantFP::get(double_, 0.0), sum_p);
	}
	Value * sum = 0;
	if(n <= 8){
		if(normalize_to){
			sum = Builder.CreateLoad(sum_p);
		}
		for(size_t i=0; i<n; ++i){
			BasicBlock * BB_lt0 = BasicBlock::Create(context, "sz_lt0", f);
			BasicBlock * BB_gt0 = BasicBlock::Create(context, "sz_gt0", f);
			BasicBlock * BB_next = BasicBlock::Create(context, "sz_next", f);
			Value * p_data_i = Builder.CreateGEP(histo, ConstantInt::get(i32_, i), "p_histo_i");
			Value * data_i = Builder.CreateLoad(p_data_i, "histo_i");
			Value * data_gtr0 = Builder.CreateFCmpOGE(data_i, ConstantFP::get(double_, 0.0), "data_i_gtr0");
			Builder.CreateCondBr(data_gtr0, BB_gt0, BB_lt0);

			Builder.SetInsertPoint(BB_lt0);
			Builder.CreateStore(ConstantFP::get(double_, 0.0), p_data_i);
			Builder.CreateBr(BB_next);

			Builder.SetInsertPoint(BB_gt0);
			Value * sum_gt0 = 0;
			if(normalize_to){
				sum_gt0 = Builder.CreateFAdd(sum, data_i);
			}
			Builder.CreateBr(BB_next);

			Builder.SetInsertPoint(BB = BB_next);
			if(normalize_to){
				PHINode * s = Builder.CreatePHI(double_, 2, "s");
				s->addIncoming(sum_gt0, BB_gt0);
				s->addIncoming(sum, BB_lt0);
				sum = s;
			}
		}
		if(normalize_to){
			Builder.CreateStore(sum, sum_p);
		}
	}
	else{
		BasicBlock * BB_loop_header = BasicBlock::Create(context, "sz_loop_header", f);
		BasicBlock * BB_loop_body = BasicBlock::Create(context, "sz_loop_body", f);
		BasicBlock * BB_loop_inc = BasicBlock::Create(context, "sz_loop_inc", f);
		BasicBlock * BB_loop_exit = BasicBlock::Create(context, "sz_loop_exit", f);

		Builder.CreateBr(BB_loop_header);
		Builder.SetInsertPoint(BB_loop_header);
		PHINode * i_phi = Builder.CreatePHI(i32_, 2, "i");
		i_phi->addIncoming(ConstantInt::get(i32_, 0), BB_in);
		// add the i_phi incoming value for next loop value later
		Value * i_eq_n = Builder.CreateICmpEQ(i_phi, ConstantInt::get(i32_, n));
		Builder.CreateCondBr(i_eq_n, BB_loop_exit, BB_loop_body);

		Builder.SetInsertPoint(BB_loop_body);
		BasicBlock * BB_gt0 = BasicBlock::Create(context, "sz_gt0", f, BB_loop_exit);
		BasicBlock * BB_lt0 = BasicBlock::Create(context, "sz_lt0", f, BB_loop_exit);
		Value * p_data_i = Builder.CreateGEP(histo, i_phi, "p_histo_i");
		Value * data_i = Builder.CreateLoad(p_data_i, "histo_i");
		Value * data_gtr0 = Builder.CreateFCmpOGE(data_i, ConstantFP::get(double_, 0.0), "data_i_gtr0");
		Builder.CreateCondBr(data_gtr0, BB_gt0, BB_lt0);

		Builder.SetInsertPoint(BB_lt0);
		Builder.CreateStore(ConstantFP::get(double_, 0.0), p_data_i);
		Builder.CreateBr(BB_loop_inc);

		Builder.SetInsertPoint(BB_gt0);
		if(normalize_to){
			Value * sum = Builder.CreateLoad(sum_p);
			sum = Builder.CreateFAdd(sum, data_i);
			Builder.CreateStore(sum, sum_p);
		}
		Builder.CreateBr(BB_loop_inc);

		Builder.SetInsertPoint(BB_loop_inc);
		Value * next_i = Builder.CreateAdd(i_phi, ConstantInt::get(i32_, 1), "next_i");
		i_phi->addIncoming(next_i, BB_loop_inc);
		Builder.CreateBr(BB_loop_header);

		Builder.SetInsertPoint(BB = BB_loop_exit);
	}
	if(normalize_to){
		Value * sum = Builder.CreateLoad(sum_p);
		Value * sum_gt0 = Builder.CreateFCmpOGT(sum, ConstantFP::get(double_, 0.0), "sum_gt0");
		BasicBlock * bb_gt0 = BasicBlock::Create(context, "bb_gt0", f);
		BasicBlock * bb_lt0 = BasicBlock::Create(context, "bb_lt0", f);
		BasicBlock * bb_exit = BasicBlock::Create(context, "bb_exit", f);
		Builder.CreateCondBr(sum_gt0, bb_gt0, bb_lt0);

		Builder.SetInsertPoint(BB = bb_gt0);
		Value * coeff_p = Builder.CreateAlloca(double_, ConstantInt::get(i32_, 1), "coeff_p");
		Builder.CreateStore(ConstantFP::get(double_, *normalize_to), coeff_p);
		Value * vcoeff = Builder.CreateLoad(coeff_p);
		vcoeff = Builder.CreateFDiv(vcoeff, sum);
		BB = bb_gt0 = emit_multiply(BB, vcoeff, histo, n);
		Builder.SetInsertPoint(BB);
		Builder.CreateBr(bb_exit); // closes BB = bb_gt0

		Builder.SetInsertPoint(bb_lt0);
		Builder.CreateBr(bb_exit); // closes bb_lt0

		Builder.SetInsertPoint(BB = bb_exit);
		llvm::PHINode * phi_coeff = Builder.CreatePHI(double_, 2);
		phi_coeff->addIncoming(vcoeff, bb_gt0);
		phi_coeff->addIncoming(ConstantFP::get(double_, 0.0), bb_lt0);
		coeff = phi_coeff;
	}
	return BB;
}


BasicBlock * llvm_module::emit_copy_ddata(BasicBlock * BB_in, Value * dest_, Value * src_, size_t n){
	LLVMContext & context = module->getContext();
	IRBuilder<> Builder(context);
	BasicBlock * BB = BB_in;
	llvm::Function * f = BB_in->getParent();
	Builder.SetInsertPoint(BB);
	Value * dest = Builder.CreateBitCast(dest_, double2_->getPointerTo(), "dest");
	Value * src = Builder.CreateBitCast(src_, double2_->getPointerTo(), "src");
	// process 2 doubles at a time, so make sure to process enough pairs:
	const size_t n_padded = n + (n % 2);
	if(n<=8){
		for(size_t i=0; i<n_padded/2; ++i){
			Value * src_index = Builder.CreateGEP(src, ConstantInt::get(i32_, i));
			Value * dest_index = Builder.CreateGEP(dest, ConstantInt::get(i32_, i));
			Value * src_value = Builder.CreateLoad(src_index);
			Builder.CreateStore(src_value, dest_index);
		}
	}
	else{
		BasicBlock * BB_loop_header = BasicBlock::Create(context, "memcpy_loop_header", f);
		BasicBlock * BB_loop_body = BasicBlock::Create(context, "memcpy_loop_body", f);
		BasicBlock * BB_loop_exit = BasicBlock::Create(context, "memcpy_loop_exit", f);
		Builder.CreateBr(BB_loop_header); // closes BB

		Builder.SetInsertPoint(BB_loop_header);
		PHINode * i_phi = Builder.CreatePHI(i32_, 2, "i");
		i_phi->addIncoming(ConstantInt::get(i32_, 0), BB);
		Value * i_eq_n = Builder.CreateICmpEQ(i_phi, ConstantInt::get(i32_, n_padded / 2));
		Builder.CreateCondBr(i_eq_n, BB_loop_exit, BB_loop_body); // closes loop_header

		Builder.SetInsertPoint(BB_loop_body);
		Value * p_src_i = GetElementPtrInst::Create(src, i_phi, "p_src_i", BB_loop_body);
		Value * src_i = Builder.CreateLoad(p_src_i, "src_i");
		Value * p_dest_i = GetElementPtrInst::Create(dest, i_phi, "p_dest_i", BB_loop_body);
		Builder.CreateStore(src_i, p_dest_i);
		Value * next_i = Builder.CreateAdd(i_phi, ConstantInt::get(i32_, 1), "next_i");
		i_phi->addIncoming(next_i, BB_loop_body);
		Builder.CreateBr(BB_loop_header); // closes loop_body

		BB = BB_loop_exit;
	}
	return BB;
}


llvm_module::llvm_module(const ParIds & pids_): pids(pids_), pid_to_index(create_pid_to_index(pids)){
   LLVMContext & context = getGlobalContext();
   module =  new Module("llvm_module", context);
   std::string err;
   ee =  ExecutionEngine::createJIT(module, &err, 0, CodeGenOpt::Aggressive);
   if(!ee){
       throw Exception("llvm_module: could not create llvm::ExecutionEngine: " + err);
   }
   //define function prototype for the codegen_*_evaluate functions:
   //double codegen_f_evaluate(llvm_module* mod, void* fptr, const double * par_value);
   double_ = Type::getDoubleTy(context);
   double2_ = VectorType::get(double_, 2);
   void_ = Type::getVoidTy(context);
   i32_ = Type::getInt32Ty(context);

   vector<Type*> arg_types(2);
   arg_types[0] = double_;
   arg_types[1] = i32_;
   llvm::FunctionType * dump_info_ty = FunctionType::get(void_, arg_types, false);
   f_dump_info = llvm::Function::Create(dump_info_ty, llvm::Function::ExternalLinkage, "dump_info", module);

   // TODO: maybe get rid of indirection ...
   Type * p_char_t = Type::getInt8Ty(context)->getPointerTo();
   arg_types.resize(3);
   arg_types[0] = arg_types[1] = p_char_t;
   arg_types[2] = double_->getPointerTo();
   FunctionType * FT = FunctionType::get(double_, arg_types, false);
   llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "codegen_f_evaluate", module);
   arg_types.resize(5);
   arg_types[2] = double_;
   arg_types[3] = arg_types[4] = double_->getPointerTo();
   FunctionType * FT2 = FunctionType::get(void_, arg_types, false);
   llvm::Function::Create(FT2, llvm::Function::ExternalLinkage, "codegen_hf_add_with_coeff", module);
   
   // create an alias 'exp_function' which either redirects to amd_exp (if it exists) or to exp:
   arg_types.resize(1);
   arg_types[0] = double_;
   FunctionType * FT_exp = FunctionType::get(double_, arg_types, false);
   void * handle = dlopen(0, 0);
   bool amd_exp_available = handle && dlsym(handle, "amd_exp");
   llvm::Function * f_exp = llvm::Function::Create(FT_exp, llvm::Function::ExternalLinkage, amd_exp_available?"amd_exp":"exp", module);
   llvm::Function * f_exp_function = llvm::Function::Create(FT_exp, llvm::Function::ExternalLinkage, "exp_function", module);
   {
	   IRBuilder<> Builder(context);
	   BasicBlock * BB = BasicBlock::Create(context, "entry", f_exp_function);
	   Builder.SetInsertPoint(BB);
	   Value * ret = Builder.CreateCall(f_exp, f_exp_function->arg_begin());
	   Builder.CreateRet(ret);
   }
   f_exp_function->addFnAttr(Attribute::AlwaysInline);
   theta_assert(module->getFunction("exp_function") != 0);
   emit_add_with_coeff_function();
}

size_t llvm_module::get_index(const theta::ParId & pid) const{
    map<theta::ParId, size_t>::const_iterator it = pid_to_index.find(pid);
    if(it==pid_to_index.end()) throw invalid_argument("llvm_module::get_index");
    return it->second;
}

void* llvm_module::getFunctionPointer(llvm::Function * function){
	theta_assert(function);
    return ee->getPointerToFunction(function);
}

llvm_module::~llvm_module(){
    delete ee;
}

void llvm_module::optimize(){
    PassManager pm;
    pm.add(createFunctionInliningPass());
    pm.add(createCFGSimplificationPass());
    //TODO: find out good combinations ...
    //pm.add(createPromoteMemoryToRegisterPass());
    //pm.add(createInstructionCombiningPass());
    //pm.add(createIndVarSimplifyPass());
    //pm.add(createLoopUnrollPass());
    //pm.add(createInstructionCombiningPass());
    //pm.add(createLoopSimplifyPass());
    //pm.add(createLoopStrengthReducePass());
    //pm.add(createBBVectorizePass());
    //pm.add(createLoopInstSimplifyPass());
    //pm.add(createGVNPass());
    //pm.add(createInstructionSimplifierPass());
    pm.add(createDeadCodeEliminationPass());
    pm.add(createCodeGenPreparePass());
    pm.run(*module);
}



// get llvm function types for the functions created by HistogramFunction and Function code
// generators:     double coeff, double * par_values, double * data
// the generated function should perform:
//     data += coeff * eval(par_values)
// where eval is the result of the HistogramFunction evaluation ...
llvm::FunctionType * get_ft_hf_add_with_coeff(llvm_module & mod){
    LLVMContext & context = mod.module->getContext();
    Type * double_t = Type::getDoubleTy(context);
    std::vector<Type*> arg_types(3);
    arg_types[0] = double_t;
    arg_types[1] = arg_types[2] = double_t->getPointerTo();
    return FunctionType::get(double_t, arg_types, false);
}

llvm::FunctionType * get_ft_function_evaluate(llvm_module & mod){
    LLVMContext & context = mod.module->getContext();
    Type * double_t = Type::getDoubleTy(context);
    std::vector<Type*> arg_types(1);
    arg_types[0] = double_t->getPointerTo();
    return FunctionType::get(double_t, arg_types, false);
}

void llvm_verify(llvm::Function* f, const std::string & fname){
    if(verifyFunction(*f, llvm::PrintMessageAction)){
        f->dump();
        throw invalid_argument("llvm_verify failed for function " + fname);
    }
}

GlobalVariable * llvm_module::add_global_ddata(const double * data, size_t n, const std::string & name, bool const_){
    size_t n_orig = n;
    if(n%2)++n;
    LLVMContext & context = module->getContext();
    std::vector<Constant*> elements(n);
    for(size_t i=0; i<n_orig; ++i){
        elements[i] = ConstantFP::get(context, APFloat(data[i]));
    }
    if(n_orig % 2){
        elements[n_orig] = ConstantFP::get(context, APFloat(0.0));
    }
    ArrayType* typ = ArrayType::get(Type::getDoubleTy(context), n);
    Constant * ini = ConstantArray::get(typ, elements);
    GlobalVariable * gvar = new GlobalVariable(*module, typ, const_, GlobalValue::PrivateLinkage, ini, name);
    gvar->setAlignment(16);
    return gvar;
}

llvm::Function * llvm_generic_codegen(const theta::Function * f, llvm_module & mod, const std::string & prefix){
    FunctionType * FT = get_ft_function_evaluate(mod);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_evaluate", mod.module);
    LLVMContext & context = mod.module->getContext();
    IRBuilder<> Builder(context);
    llvm::Function * llvm_codegen_f_evaluate = mod.module->getFunction("codegen_f_evaluate");
    Type * p_char_t = Type::getInt8Ty(context)->getPointerTo();
    Type * i64_t = Type::getInt64Ty(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    Value * par_values = &*(F->arg_begin());
    Value * mod_ = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(&mod)), p_char_t);
    Value * fptr = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(f)), p_char_t);
    Value * ret = Builder.CreateCall3(llvm_codegen_f_evaluate, mod_, fptr, par_values);
    Builder.CreateRet(ret);
    return F;
}

llvm::Function * llvm_generic_codegen(const theta::HistogramFunction * hf, llvm_module & mod, const std::string & prefix){
    FunctionType * FT = get_ft_hf_add_with_coeff(mod);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, prefix + "_add_with_coeff", mod.module);
    LLVMContext & context = mod.module->getContext();
    IRBuilder<> Builder(context);
    llvm::Function * llvm_codegen_hf_add_with_coeff = mod.module->getFunction("codegen_hf_add_with_coeff");
    Type * p_char_t = Type::getInt8Ty(context)->getPointerTo();
    Type * i64_t = Type::getInt64Ty(context);
    Type * double_t = Type::getDoubleTy(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    llvm::Function::arg_iterator iter = F->arg_begin();
    Value * coeff = iter++;
    Value * par_values = iter++;
    Value * data = iter;
    Value * mod_ = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(&mod)), p_char_t);
    Value * hfptr = Builder.CreateBitCast(ConstantInt::get(i64_t, reinterpret_cast<unsigned long>(hf)), p_char_t);
    Value *Args[] = { mod_, hfptr, coeff, par_values, data };
    Builder.Insert(CallInst::Create(llvm_codegen_hf_add_with_coeff, ArrayRef<Value*>(Args, Args+5)), "");
    Builder.CreateRet(ConstantFP::get(double_t, 1.0));
    return F;
}

llvm::Function * create_llvm_function(const theta::Function * f, llvm_module & mod, const std::string & prefix){
    llvm::Function * result;
    if(const llvm_enabled_function * en = dynamic_cast<const llvm_enabled_function*>(f)){
        result = en->llvm_codegen(mod, prefix);
    }
    else{
        result = llvm_generic_codegen(f, mod, prefix);
    }
    result->addFnAttr(Attribute::AlwaysInline);
    result->setLinkage(GlobalValue::PrivateLinkage);
    return result;
}



llvm::Function * create_llvm_histogram_function(const theta::HistogramFunction * hf, llvm_module & mod, const std::string & prefix){
    llvm::Function * result;
    if(const llvm_enabled_histogram_function* en = dynamic_cast<const llvm_enabled_histogram_function*>(hf)){
        result = en->llvm_codegen(mod, prefix);
    }
    else if(hf->get_parameters().size()==0){
        Histogram1D h0;
        hf->apply_functor(copy_to<Histogram1D>(h0), ParValues());
        llvm::FunctionType * FT = get_ft_hf_add_with_coeff(mod);
        result = llvm::Function::Create(FT, llvm::Function::PrivateLinkage, prefix + "_add_with_coeff", mod.module);
        llvm::Function::arg_iterator iter = result->arg_begin();
        Value * coeff = iter++;
        /*Value * par_values =*/ iter++; // not used ...
        Value * data = iter;
        LLVMContext & context = mod.module->getContext();
        Type * double_t = Type::getDoubleTy(context);
        BasicBlock * BB = BasicBlock::Create(context, "entry", result);
        BB = mod.emit_add_with_coeff(BB, coeff, data, h0);
        IRBuilder<> Builder(context);
        Builder.SetInsertPoint(BB);
        Builder.CreateRet(ConstantFP::get(double_t, 1.0));
        //mod.module->dump();
    }else{
        result = llvm_generic_codegen(hf, mod, prefix);
    }
    result->addFnAttr(Attribute::AlwaysInline);
    //llvm_verify(result, prefix);
    return result;
}



namespace{
    
    class add_to_vdouble: public functor<Histogram1D>{
    private:
        mutable double * data;
        double coeff;
    public:
        add_to_vdouble(double * data_, double coeff_): data(data_), coeff(coeff_){}
        virtual void operator()(const Histogram1D & h) const{
            utils::add_fast_with_coeff(data, h.get_data(), coeff, h.size());
        }
    };
    
}

// the callback functions llvm will call for Functions and HistogramFunctions which are not derived from llvm_enabled<T>:
extern "C" {
double codegen_f_evaluate(llvm_module* mod, void* fptr, const double * par_value);
void codegen_hf_add_with_coeff(llvm_module* mod, void* hfptr, double coeff, const double * par_value, double * data);
}

double codegen_f_evaluate(llvm_module* mod, void* fptr, const double * par_values){
    ParValues values(par_values, mod->get_parameters());
    return static_cast<theta::Function*>(fptr)->operator()(values);
}

void codegen_hf_add_with_coeff(llvm_module* mod, void* hfptr, double coeff, const double * par_values, double * data){
    ParValues values(par_values, mod->get_parameters());
    static_cast<theta::HistogramFunction*>(hfptr)->apply_functor(add_to_vdouble(data, coeff), values);
}


/*  llvm_enable_function */
double llvm_enable_function::operator()(const theta::ParValues & values) const{
    std::vector<double> v_values(par_ids.size());
    size_t i=0;
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
        v_values[i] = values.get(*it);
    }
    return llvmf(&v_values[0]);
}

llvm_enable_function::llvm_enable_function(const theta::Configuration & cfg): f(PluginManager<theta::Function>::build(Configuration(cfg, cfg.setting["function"]))),
   m(f->get_parameters()){
    par_ids = f->get_parameters();
    llvm::Function * tmp_llvmf = create_llvm_function(f.get(), m, "enablef");
    llvmf = reinterpret_cast<t_function_evaluate>(m.getFunctionPointer(tmp_llvmf));
    if(llvmf==0){
        throw ConfigurationException("could not compile llvm");
    }
}

/*  llvm_enable_histogram_function */
void llvm_enable_histogram_function::apply_functor(const functor<Histogram1D> & f, const theta::ParValues & values) const{
    std::vector<double> v_values(par_ids.size());
    size_t i=0;
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
        v_values[i] = values.get(*it);
    }
    h0.set_all_values(0.0);
    llvmhf(1.0, &v_values[0], h0.get_data());
    f(h0);
}

void llvm_enable_histogram_function::apply_functor(const functor<Histogram1DWithUncertainties> & f, const theta::ParValues & values) const{
    Histogram1D h;
    apply_functor(copy_to<Histogram1D>(h), values);
    f(Histogram1DWithUncertainties(h));
}

llvm_enable_histogram_function::llvm_enable_histogram_function(const theta::Configuration & cfg):
  hf(PluginManager<theta::HistogramFunction>::build(Configuration(cfg, cfg.setting["histogram_function"]))), m(hf->get_parameters()){
  par_ids = hf->get_parameters();
  size_t nbins;
  double xmin, xmax;
  hf->get_histogram_dimensions(nbins, xmin, xmax);
  h0 = Histogram1D(nbins, xmin, xmax);
  llvm::Function * tmp_llvmf = create_llvm_histogram_function(hf.get(), m, "enablehf");
  llvmhf = reinterpret_cast<t_hf_add_with_coeff>(m.getFunctionPointer(tmp_llvmf));
  if(llvmhf==0){
      throw ConfigurationException("could not compile llvm");
  }
}

REGISTER_PLUGIN(llvm_enable_function)
REGISTER_PLUGIN(llvm_enable_histogram_function)

    
