#ifndef LLVM_LLVM_INTERFACE_HPP
#define LLVM_LLVM_INTERFACE_HPP

#include <string>
#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/phys.hpp"
#include "interface/histogram-function.hpp"

#include "boost/optional.hpp"

#include "llvm/DerivedTypes.h"
#include "llvm/Function.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"
#include "llvm/Support/IRBuilder.h"

namespace llvm{
class ExecutionEngine;
}

typedef double (*t_function_evaluate)(const double *);
typedef double (*t_hf_add_with_coeff)(double, const double *, double *);

// hold all relevant information about a LLVM module, including any
// data allocations, etc.
// The lifetime of this object must be guaranteed as long as the llvm module
// is used (functions there called, etc.)
class llvm_module{
private:
    const theta::ParIds pids;
    const std::map<theta::ParId, size_t> pid_to_index;
    llvm::ExecutionEngine * ee;

    llvm::Type * double_,  * double2_, * void_,  * i32_;

    llvm::Function * f_dump_info;

    void emit_add_with_coeff_function();

public:
    // making the module public is probably nor very nice. On the other hand, re-inventing and proxying
    // the whole llvm interface isn't either.
    llvm::Module * module;

    explicit llvm_module(const theta::ParIds & pids);
    
    // get the index to use in llvm functions into par_values for the given ParId.
	size_t get_index(const theta::ParId & pid) const;

	const theta::ParIds & get_parameters() const{
		return pids;
	}

	// the emit_* functions ha sthe convention
	// that they start code emission at the end of the given BasicBlock "in", and they return
	// a BasicBlock where code insertion of the caller should continue (might be the same block).
	//
	// The typical call at the site is:
	// BasicBlock * BB = ...
	//...
	// BB = emit_..(BB);
	// Builder.SetInsertPoint(BB);
	//
	// The returned BasicBlock remains "unfinished", i.e. no branch is yet inserted.

    // emits code which normalizes the given histogram:
    // 1. truncates values below zero
    // 2. optionally multiplies with an overall factor such that the sum of the elements
    //    is the given double. In this case, sets the coefficient coeff to the coefficient used in the normalization.
    llvm::BasicBlock * emit_normalize_histo(llvm::BasicBlock * in, llvm::Value * histo_, size_t n,
            const boost::optional<double> & normalize_to, llvm::Value *& coeff);

    // emit code in the BasicBlock to multiply the given histogram with the given coefficient.
    llvm::BasicBlock * emit_multiply(llvm::BasicBlock * in, llvm::Value * coeff, llvm::Value * histo_, size_t n);

    // emits code for data_out[i] += coeff * h.get(i)  for i in 0...h.size()
    llvm::BasicBlock * emit_add_with_coeff(llvm::BasicBlock * in, llvm::Value * coeff, llvm::Value * data_out, const theta::DoubleVector & h);

    llvm::BasicBlock * emit_add_with_coeff(llvm::BasicBlock * in, llvm::Value * coeff, llvm::Value * data_out, llvm::Value * data_in, size_t n);


    // dest and src are pointer to double and have to be aligned. In case of odd n, can still copy even number of bytes
    llvm::BasicBlock* emit_copy_ddata(llvm::BasicBlock * in, llvm::Value * dest_, llvm::Value * src_, size_t n);

    llvm::GlobalVariable * add_global_ddata(const double * data, size_t n, const std::string & name, bool const_);

    //void emit_dump(llvm::Value * d_, llvm::Value * int_, llvm::IRBuilder<> & Builder);

    // run llvm optimizations on all the code generated so far
    void optimize();
    
    // compile module and get function pointer
    void* getFunctionPointer(llvm::Function * function);
    

    ~llvm_module();
};


/** \brief Abstract base class for plugins supporting llvm code generation
 *
 * Only used for T=Function and T=HistogramFunction
 * 
 * Note: making this derived from T [= Function, HistogramFunction], so we can
 * try to dynamic_cast pointers of T down to llvm_enabled<T> during setup.
 */
class llvm_enabled_function: public theta::Function{
public:
    virtual llvm::Function * llvm_codegen(llvm_module & mod, const std::string & prefix) const = 0;
    virtual ~llvm_enabled_function(){}
};


class llvm_enabled_histogram_function: public theta::HistogramFunction{
public:
    virtual llvm::Function * llvm_codegen(llvm_module & mod, const std::string & prefix) const = 0;
    // the histogram with squared uncertainties. It will be added by the model in case of bb_uncertainties = true
    virtual theta::Histogram1D get_uncertainty2_histogram() const = 0;
    virtual ~llvm_enabled_histogram_function(){}
};

/** \brief Function plugin wrapper to execute Functions within llvm
 *
 * This plugin is usually not needed explicitly, as for llvm code generation all non-llvm Functions will
 * be automatically wrapped if used in a llvm_model.
 * This class is provided for testing purposes only (to avoid the test program explicitly depends on the llvm plugin).
 *
 * It generated llvm code in its own module. This will either generate native llvm-code or a wrapper, depending
 * on whether or not the specified function is "llvm_enabled".
 *
 * configuration via:
 * \code
 * some_function = {
 *     type = "llvm_enable_function";
 *     function = { / *  specification of the Function to wrap  * /};
 * };
 * \endcode
 */
class llvm_enable_function: public theta::Function {
private:
   std::auto_ptr<theta::Function> f;
   llvm_module m;
   t_function_evaluate llvmf; // pointer to the llvm-generated function, called by operator()
public:
   virtual double operator()(const theta::ParValues & values) const;
   llvm_enable_function(const theta::Configuration & cfg);
};


/** \brief HistogramFunction plugin wrapper to execute HistogramFunctions within llvm
 *
 * Same as llvm_enable_function but for HistogramFunctions.
 *
 * Configure via:
 * \code
 * {
 *   type = "llvm_enable_histogram_function";
 *   histogram_function = { / * some histogram function to wrap * / };
 * };
 * \endcode
 */
class llvm_enable_histogram_function: public theta::HistogramFunction {
private:
   std::auto_ptr<theta::HistogramFunction> hf;
   llvm_module m;
   t_hf_add_with_coeff llvmhf; // pointer to the llvm-generated function, called by operator()
   mutable theta::Histogram1D h0;
public:
    /// note: uncertainties are not supported, all uncertainties will be set to zero!
   virtual void apply_functor(const theta::functor<theta::Histogram1DWithUncertainties> & f, const theta::ParValues & values) const;
   virtual void apply_functor(const theta::functor<theta::Histogram1D> & f, const theta::ParValues & values) const;
   
   virtual void get_histogram_dimensions(size_t & n, double & xmin, double & xmax) const{
       hf->get_histogram_dimensions(n, xmin, xmax);
   }
   llvm_enable_histogram_function(const theta::Configuration & cfg);
};

// To generate the code for a Function / HistogramFunction, always use these create_llvm_* methods. They will either
// call the object's llvm_codegen function (if it inherits from llvm_enabled<...>) or generate the generic one
// which issue a callback to the operator()().
//
// Also calls a verifier for the function which throws a FatalException if not successful
//
// Note that for constant HistogramFunctions, will always emit the (non-callback) llvm code; the supplied HistogramFunction
// does not need to inherit from llvm_enabled<HistogramFunction>.
llvm::Function * create_llvm_function(const theta::Function * f, llvm_module & mod, const std::string & prefix);
llvm::Function * create_llvm_histogram_function(const theta::HistogramFunction * f, llvm_module & mod, const std::string & prefix);

// get the llvm function types the llvm_enabled<Function>::llvm_codegen and llvm_enabled<HistogramFunction>::llvm_codegen should
// generate.
llvm::FunctionType * get_ft_function_evaluate(llvm_module & mod);
llvm::FunctionType * get_ft_hf_add_with_coeff(llvm_module & mod);

extern "C" void dump_info(double d, int i);


#endif

