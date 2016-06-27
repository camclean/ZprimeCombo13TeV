#ifndef ROOT_ROOT_HISTOGRAM_HPP
#define ROOT_ROOT_HISTOGRAM_HPP

#include "interface/histogram-function.hpp"

/** \brief Plugin to read Histogram from a root file
 *
 * Configuration: anywhere, where a (constant) Histogram has to be defined,
 * use a setting like:
 * \code
 * {
 * type = "root_histogram";
 * filename = "path/to/file.root";
 * histoname = "histoname-in-file";
 * normalize_to = 1.0; //optional; default is not to scale th histogram
 * use_errors = true; // optional; default is false
 * rebin = 2; // optional; default is 1.
 * range = (0.0, 50.); // optional; default is to use the whole range
 * zerobin_fillfactor = 0.0001; // optional, default is 0.0
 * }
 * \endcode
 *
 * Returns a ConstantHistogramFuntion which contains the specified ROOT histogram from the
 * given root file. Underflow and overflow are usually not copied, unless the range setting (see below)
 * explicitely includes them.
 *
 * If \c use_errors is true, the errors in the root Histogram are used. If it is \c false (the default),
 * all uncertainties are set to 0.
 *
 * If rebin is given, TH1::Rebin will be called with this value.
 *
 * If range is given, only bins within this range are copied, not the whole histogram.
 * This can also be used to explicitely include the overflow or underflow bins by specifying
 * a range which goes beyond the border of the histogram. If the range is within the histogram range,
 * the range borders must coincide with bin borders. Otherwise, an exception is thrown. Note that
 * rebinning is done before the range is applied, so make sure to specify a range which is valid after
 * rebinning.
 *
 * Bin entries of exactly zero can be problematic if dicing pseudo data from templates with a non-zero entry (for example,
 * because dicing pseudo data from templates affected by a systematic uncertainty) while calculating
 * the likelihood function for a model which has always zero prediction in that bin (for example, because the likelihood is calculated
 * for a model without systematic uncertainties). Then, the likelihood function is exactly zero and negative log-likelihood is infinity
 * which is a problem if you want to minimize this function. The best way is to prevent zero-prediction is
 * in the phase of template construction where one should avoid zero-bins altogether.
 * The parameter \c zerobin_fillfactor provides a "quick-and-dirty" workaround for this problem: 
 * each bin entry is set to a value at least \c zerobin_fillfactor * integral where integral is the sum of all bins.
 *
 * Note that normalize_to is applied first to the histogram excluding under / overflow, then the rebinning is done, and then
 * the range setting is applied and then the zero bins are filled.
 *
 * It is also possible to read in multidimensional Histograms (TH2 / TH3) with this plugin. In this case, however,
 * <ul>
 *   <li>it is not allowed to specify rebin or range (otherwise, reading the histogram will fail). All bins except overflow / uncerflow will
 *      be used.</li>
 *   <li>internally, the range will always be set to [0, nbinsx*nbinsy*nbinsz]. Usually, this is not relevant, except if you want to
 *     mix root_histograms and Histograms created by other theta plugins and have to make sure that the ranges are the same. </li>
 * </ul>
 *
 * \sa ConstantHistogramFunctionError ConstantHistogramFunction
 */
class root_histogram: public theta::ConstantHistogramFunction{
public:
    /// Constructor used by the plugin system
    root_histogram(const theta::Configuration & ctx);
};

#endif
