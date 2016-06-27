#ifndef VARIABLES_UTILS_HPP
#define VARIABLES_UTILS_HPP

#include "interface/decls.hpp"

namespace theta{
      /** \brief Populate VarIdManager from a Setting.
       *
       * This function uses the "observables" and "parameters" setting groups in cfg.setting to
       * populate the VarIdManager cfg.vm.
       *
       * Assumes \c cfg.settings is a setting group as encountered at the 'top level' of a configuration like this:
       * \code
       *  //...
       *  observables = {
       *     mass = {
       *        nbins = 200;
       *        range = [0.0, 10.0];
       *    };
       *  };
       *
       *  //...
       *
       *  parameters = ("p0", "p1", "p2");
       *
       *  //...
       * \endcode
       * then this function will call VarIdManager::createParId("p0", 20, 0, inf) and VarIdManager::createObsId("mass", 200, 0.0, 10.0).
       */
      void apply_vm_settings(theta::Configuration & cfg);
      
      std::ostream & operator<<(std::ostream & out, const theta::ParIds & pids);
}



#endif

