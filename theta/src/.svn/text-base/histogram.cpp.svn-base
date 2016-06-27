#include "interface/histogram.hpp"
#include "interface/random.hpp"
#include "interface/exception.hpp"
#include "interface/utils.hpp"

#include <cmath>
#include <sstream>
#include <limits>
#include <new>

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif


using namespace theta;
using std::invalid_argument;

namespace{
  double * allocate_doubles(size_t n){
      double * result = 0;
      const size_t n_orig = n;
      //always allocate an even number of bins:
      if(n_orig % 2) ++n;
      //for the add_fast routine, which might use SSE optimizations, we need this alignment. And
      // while we at it, we should make sure double is as expected:
      BOOST_STATIC_ASSERT(sizeof(double)==8);
      #ifdef __APPLE__
      // note: Mac OS does not have memalign, but does 16-byte alignment anyway. This is all we need, so just use malloc directly in this case:
      result = reinterpret_cast<double*>(malloc(sizeof(double)*n));
      #else
      result = reinterpret_cast<double*>(memalign(2 * sizeof(double), sizeof(double) * n)); // note: while it is not guaranteed that we can call "free" on the result, that's ok for most systems.
      #endif
      if(result==0){
        throw std::bad_alloc();
      }
      //set the extra allocated double to zero by convention.
      if(n_orig % 2) result[n-1] = 0.0;
      return result;
   }
   
   // note: data can be NULL
   void free_doubles(double * data){
      free(data);
   }
}


DoubleVector::DoubleVector(size_t n): data(0), n_data(n), n_alloc(0){
    realloc(n_data);
    set_all_values(0.0);
}

DoubleVector::~DoubleVector(){
    free_doubles(data);
}

DoubleVector::DoubleVector(const DoubleVector & rhs): data(0), n_data(rhs.n_data), n_alloc(0){
    realloc(n_data);
    std::copy(rhs.data, rhs.data + n_data, data);
}

// note: this is the only method which actually sets the data member and does the allocation.
void DoubleVector::realloc(size_t n){
    if(n == n_data && n_alloc >= n_data && data != 0) return;
    n_data = n;
    if(n_alloc < n){
        free_doubles(data);
        if(n > 0){
            data = allocate_doubles(n);
        }else{
            data = 0;
        }
        n_alloc = n;
    }
    else{
        if(n%2 > 0){
            data[n] = 0.0;
        }
    }
}

Histogram1D::Histogram1D(size_t b, double x_min, double x_max) : DoubleVector(b), xmin(x_min), xmax(x_max) {
    if(xmin >= xmax) throw invalid_argument("Histogram: xmin >= xmax not allowed");
    set_all_values(0.0);
}

void Histogram1D::multiply_with_ratio_exponented(const Histogram1D & nominator, const Histogram1D & denominator, double exponent){
   check_compatibility(nominator);
   check_compatibility(denominator);
   const double * n_data = nominator.get_data();
   const double * d_data = denominator.get_data();
   double * data = get_data();
   const size_t n = size();
   for(size_t i=0; i<n; i++){
      if(d_data[i]>0.0)
         data[i] *= pow(n_data[i] / d_data[i], exponent);
   }
}

void Histogram1D::fill(double xvalue, double weight) {
    int bin = static_cast<int> ((xvalue - xmin) * get_nbins() / (xmax - xmin));
    if (bin < 0 || (bin==0 && xvalue < xmin)) return;
    if (static_cast<size_t> (bin) >= size()) return;
    get_data()[bin] += weight;
}

void Histogram1D::fail_check_compatibility(const Histogram1D & h) const {
    std::stringstream s;
    s <<  "Histogram1D::check_compatibility: Histograms are not compatible (nbins, xmin, xmax) are: "
            " (" << get_nbins() << ", " << xmin << ", " << xmax << ") and "
            " (" << h.get_nbins() << ", " << h.xmin << ", " << h.xmax << ")";
    throw invalid_argument(s.str());
}

void Histogram1D::operator*=(const Histogram1D & h) {
    check_compatibility(h);
    const double * hdata = h.get_data();
    double * data = get_data();
    const size_t n = size();
    for (size_t i = 0; i < n; ++i) {
        data[i] *= hdata[i];
    }
}
