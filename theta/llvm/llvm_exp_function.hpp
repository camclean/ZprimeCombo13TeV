#ifndef LLVM_EXP_FUNCTION
#define LLVM_EXP_FUNCTION

#include "interface/decls.hpp"
#include "interface/phys.hpp"
#include "llvm/llvm_interface.hpp"

/** \brief Function returning exp(lambda * p) where p is a configurable parameter and lambda a literal constant
 *
 * This is the llvm_enabled version of exp_function. See documentation there.
 */
class llvm_exp_function: public llvm_enabled_function{
private:
    std::vector<theta::ParId> v_pids;
    std::vector<double> lambdas_minus, lambdas_plus;

public:
    llvm_exp_function(const theta::Configuration & cfg);
    virtual double operator()(const theta::ParValues & v) const;
    virtual llvm::Function * llvm_codegen(llvm_module & mod, const std::string & prefix) const;
};


#endif

