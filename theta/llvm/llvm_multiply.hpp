#include "interface/phys.hpp"
#include "llvm/llvm_interface.hpp"

/** \brief A Function which multiplies different factors which can be parameters, literal values and other Functions
 *
 * It is a generalization of the \link mult mult \endlink plugin and accepts more types of factors, namely
 * constant, literal values, and other Functions.
 * 
 * This is a llvm-enabled version of the \link multiply multiply \endlink plugin. See documentation there.
 */
class llvm_multiply: public llvm_enabled_function {
public:
    llvm_multiply(const theta::Configuration & cfg);
    virtual double operator()(const theta::ParValues & v) const;
    llvm::Function * llvm_codegen(llvm_module & mod, const std::string & prefix) const;
    
private:
    std::vector<theta::ParId> v_pids;
    double literal_factor;
    boost::ptr_vector<theta::Function> functions;
};

