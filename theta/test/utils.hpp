#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include "interface/plugin.hpp"
#include "interface/cfg-utils.hpp"
#include "libconfig/libconfig.h++"

#include <cmath>

void load_core_plugins();

//returns false if loading root plugins was not successful
bool load_root_plugins();
bool load_llvm_plugins();


/** \brief Equality check for floating point numbers using relative comparison
 *
 * This function checks whether a and b are "close" on the sense
 * that the relative difference is very small. In case both are zero, only exact equality will
 * lead to a return value of true.
 */
inline bool close_to_relative(double a, double b){
    return a==b || fabs(a-b) / std::max(fabs(a),fabs(b)) < 1e-15;
}

/** \brief Equality check for floating point numbers using comparison to an external scale
 *
 * This function checks whether a and b are "close"
 * compared to \c scale. Note that \c scale is not a
 * maximal tolerance, but a typical scale which was used to 
 * compute a and b which might be equal as result of this
 * computation, but round-off errors might tell you that a!=b.
 *
 * scale > 0 must hold.
 */
inline bool close_to(double a, double b, double scale){
   if(scale==0.0) return a==b;
   return fabs(a-b) / scale < 1e-15;
}

bool histos_equal(const theta::Histogram1D & h1, const theta::Histogram1D & h2);

class ConfigCreator{
public:
    
    ConfigCreator(const std::string & cfg_string, const boost::shared_ptr<theta::VarIdManager> & vm);
    
    const theta::Configuration & get(){
        return cfg;
    }
    
private:
    theta::Configuration setup_cfg(const std::string & cfg_string);
    
    boost::shared_ptr<theta::SettingUsageRecorder> rec;
    theta::Configuration cfg;
};

#endif
