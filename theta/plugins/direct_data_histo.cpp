#include "plugins/direct_data_histo.hpp"

direct_data_histo::direct_data_histo(const theta::Configuration & cfg){
   int nbins_ = cfg.setting["nbins"];
   if(nbins_<=0){
       throw theta::ConfigurationException("nbins <= 0 not allowed.");
   }
   size_t nbins = nbins_;
   double xmin = cfg.setting["range"][0];
   double xmax = cfg.setting["range"][1];
   if(cfg.setting["data"].size() != nbins){
      throw theta::ConfigurationException("The length of " + cfg.setting["data"].get_path() + " and the nbins setting are inconsistent.");
   }
   theta::Histogram1DWithUncertainties h(nbins, xmin, xmax);
   bool have_uncertainties = false;
   if(cfg.setting.exists("uncertainties")){
       if(cfg.setting["uncertainties"].size() != nbins){
           throw theta::ConfigurationException("The length of " + cfg.setting["uncertainties"].get_path() + " and the nbins setting are inconsistent.");
       }
       have_uncertainties = true;
   }
   for(unsigned int i=0; i<nbins; ++i){
       double unc = have_uncertainties ? cfg.setting["uncertainties"][i] : 0.0;
       h.set(i, cfg.setting["data"][i], unc);
   }
   set_histo(h);
}

REGISTER_PLUGIN(direct_data_histo)
