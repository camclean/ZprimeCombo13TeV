#ifndef LLVM_CUBICLINEAR_HISTOMORPH
#define LLVM_CUBICLINEAR_HISTOMORPH

#include "interface/decls.hpp"
#include "llvm/llvm_interface.hpp"
#include "interface/histogram-function.hpp"

class llvm_cubiclinear_histomorph: public llvm_enabled_histogram_function{
public:
    llvm_cubiclinear_histomorph(const theta::Configuration & cfg);
    virtual void apply_functor(const theta::functor<theta::Histogram1DWithUncertainties> & f, const theta::ParValues & values) const{
    	throw std::invalid_argument("not implemented");
    }
    virtual void apply_functor(const theta::functor<theta::Histogram1D> & f, const theta::ParValues & values) const;
    virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const;
    virtual llvm::Function * llvm_codegen(llvm_module & mod, const std::string & prefix) const;
    virtual theta::Histogram1D get_uncertainty2_histogram() const;
private:
    theta::Histogram1D h0;
    theta::Histogram1D h0_uncertainty2;
	double h0_sum;

	//the interpolation parameters used to interpolate between hplus and hminus.
	std::vector<theta::ParId> vid;
	std::vector<theta::Histogram1D> hplus_diff; // hplus_diff[i] + h0 yields hplus
	std::vector<theta::Histogram1D> hminus_diff;

	//diff and sum are the difference and sum of the hplus_diff and hminus_diff histos
	std::vector<theta::Histogram1D> diff;
	std::vector<theta::Histogram1D> sum;

	std::vector<double> parameter_factors;
	bool normalize_to_nominal;

	//intermediate histogram for evaluation:
	mutable theta::Histogram1D h;
};

#endif
