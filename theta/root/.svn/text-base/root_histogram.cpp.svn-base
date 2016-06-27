#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"
#include "root/root_histogram.hpp"

#include "TH1.h"
#include "TFile.h"

using namespace theta;
using namespace std;

root_histogram::root_histogram(const Configuration & ctx){
    string filename = utils::replace_theta_dir(ctx.setting["filename"]);
    string histoname = ctx.setting["histoname"];
    int rebin = 1;
    double range_low = -999;
    double range_high = -999;
    if(ctx.setting.exists("rebin")){
       rebin = ctx.setting["rebin"];
    }
    if(ctx.setting.exists("range")){
       range_low = ctx.setting["range"][0];
       range_high = ctx.setting["range"][1];
       if(range_low >= range_high) throw ConfigurationException("invalid 'range' given");
    }
    bool use_errors = false;
    if(ctx.setting.exists("use_errors")){
         use_errors = ctx.setting["use_errors"];
    }
    
    TFile file(filename.c_str(), "read");
    if(file.IsZombie()){
        stringstream s;
        s << "Could not open file '" << filename << "'";
       throw ConfigurationException(s.str());
    }
    TH1* histo = dynamic_cast<TH1*>(file.Get(histoname.c_str()));
    if(not histo){
       stringstream s;
       s << "Did not find TH1 '" << histoname << "' in file '" << filename << "'";
       throw ConfigurationException(s.str());
    }
    if(ctx.setting.exists("normalize_to")){
       double norm = ctx.setting["normalize_to"];
       double histo_integral = histo->Integral();
       if(histo_integral==0){
           if(norm!=0){
               throw ConfigurationException("specified non-zero 'normalize_to' setting but original Histogram's integral is zero");
           }
       }
       else{
           histo->Scale(norm / histo_integral);
       }
    }

    Histogram1DWithUncertainties h;
    //take care of 1D histograms:
    if(histo->GetDimension() == 1){
       histo->Rebin(rebin);
       // bin_low and bin_high refer to the ROOT binning convention
       int bin_low = 1;
       int bin_high = histo->GetNbinsX();
       double xmin = histo->GetXaxis()->GetXmin();
       double xmax = histo->GetXaxis()->GetXmax();
       if(range_low!=-999){
          bin_low = histo->GetXaxis()->FindBin(range_low);
          if(bin_low==0) xmin = histo->GetXaxis()->GetXmin() - histo->GetXaxis()->GetBinWidth(1);
          else xmin = histo->GetXaxis()->GetBinLowEdge(bin_low);
          if(bin_low > 0 && range_low!=histo->GetXaxis()->GetBinLowEdge(bin_low)){
              throw ConfigurationException("'range' setting incompatible with bin borders.");
          }
       }
       if(range_high!=-999){
          if(range_high > histo->GetXaxis()->GetXmax()){
              bin_high = histo->GetNbinsX()+1;
              xmax = histo->GetXaxis()->GetXmax() + histo->GetXaxis()->GetBinWidth(histo->GetNbinsX());
          }
          else{
             bin_high = histo->GetXaxis()->FindBin(range_high);
             --bin_high;
             xmax = histo->GetXaxis()->GetBinUpEdge(bin_high);
             if(range_high!=histo->GetXaxis()->GetBinUpEdge(bin_high)){
                 throw ConfigurationException("'range' setting incompatible with bin borders.");
             }
          }
       }
       int nbins = bin_high - bin_low + 1;
       h = Histogram1DWithUncertainties(nbins, xmin, xmax);
       for(int i = bin_low; i <= bin_high; i++){
          double content = histo->GetBinContent(i);
          double unc = use_errors ? histo->GetBinError(i) : 0.0;
          h.set(i - bin_low, content, unc);
       }
    }
    //take care of 2D/3D histograms:
    else{
       if(rebin!=1 || range_low!=-999 || range_high!=-999) throw ConfigurationException("rebin and range setting not allowed for multidimensional Histograms!");
       size_t nbins_x = histo->GetNbinsX();
       size_t nbins_y = histo->GetNbinsY();
       size_t nbins_z = histo->GetNbinsZ();
       const size_t n_total = nbins_x * nbins_y * max<size_t>(nbins_z, 1);
       h = Histogram1DWithUncertainties(n_total, 0, n_total);
       size_t ibin_theta=0;
       for(size_t i=1; i <= nbins_x; ++i){
          for(size_t j=1; j<=nbins_y; ++j){
              // to treat both 2D and 3D histos, use k==0 as z-index for 2D histos,
              // but skip k==0 for 3D histos, as then, this is an underflow bin
              for(size_t k=0; k<=nbins_z; ++k){
                  if(k==0 && nbins_z > 0)continue;
                  double content = histo->GetBinContent(i, j, k);
                  double unc = use_errors ? histo->GetBinError(i, j, k) : 0.0;
                  h.set(ibin_theta, content, unc);
                  ibin_theta++;
              }
          }
       }
        
    }
    
    //apply zerobin_fillfactor:
    double zerobin_fillfactor = 0.0;
    if(ctx.setting.exists("zerobin_fillfactor")){
        zerobin_fillfactor = ctx.setting["zerobin_fillfactor"];
        if(zerobin_fillfactor < 0){
           throw ConfigurationException("zerobin_fillfactor must be >= 0.0!");
        }
        double integral = 0.0;
        for(size_t i=0; i<h.get_nbins(); ++i){
            integral += h.get_value(i);
        }
        // the minimum value to set all histogram bin to:
        const double min_val = integral * zerobin_fillfactor / h.get_nbins();
        for(size_t i=0; i<h.get_nbins(); ++i){
            double unc = h.get_uncertainty(i);
            double val = h.get_value(i);
            // set uncertainty to (at least) the difference of new and old value, if we have fill the histogram
            if(val < min_val){
                unc = max(unc, min_val - val);
                val = min_val;
            }
            h.set(i, val, unc);
        }
    }
    set_histo(h);
}

REGISTER_PLUGIN(root_histogram)
