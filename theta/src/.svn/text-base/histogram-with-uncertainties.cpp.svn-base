#include "interface/histogram-with-uncertainties.hpp"
#include "interface/exception.hpp"

#include <sstream>


using namespace theta;
using namespace std;

Histogram1DWithUncertainties::Histogram1DWithUncertainties(size_t bins, double xmin_, double xmax_): xmin(xmin_), xmax(xmax_),
   values(bins), nontrivial_unc(false){
   if(xmin >= xmax) throw invalid_argument("Histogram1DWithUncertainty: xmin >= xmax not allowed");
}

void Histogram1DWithUncertainties::set_nontrivial_unc(){
    if(!nontrivial_unc){
        nontrivial_unc = true;
        // we should not have the uncertainty vector allocated yet:
        theta_assert(sq_uncertainties.size()==0);
        sq_uncertainties = DoubleVector(values.size());
    }
    theta_assert(sq_uncertainties.size()==values.size());
}

void Histogram1DWithUncertainties::fail_check_compatibility(const Histogram1DWithUncertainties & h) const {
    std::stringstream s;
    s <<  "Histogram1DWithUncertainties::check_compatibility: Histograms are not compatible (nbins, xmin, xmax) are: "
            " (" << get_nbins() << ", " << xmin << ", " << xmax << ") and "
            " (" << h.get_nbins() << ", " << h.xmin << ", " << h.xmax << ")";
    throw invalid_argument(s.str());
}

Histogram1DWithUncertainties::Histogram1DWithUncertainties(const Histogram1D & h): xmin(h.get_xmin()), xmax(h.get_xmax()), values(h), nontrivial_unc(false){
}
