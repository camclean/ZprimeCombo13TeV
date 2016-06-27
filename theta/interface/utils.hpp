#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <string>
#include <algorithm>

#ifdef USE_AMDLIBM
#include "amdlibm/include/amdlibm.h"
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif



namespace theta { namespace utils{

extern std::string theta_dir;
void fill_theta_dir(char** argv);

/// Replaces the string "$THETA_DIR" by the theta directory; to be used by plugins resolving filenames
std::string replace_theta_dir(const std::string & path);

double phi_inverse(double p);

/** \brief Calculate the roots of a quadratic equation
 * 
 * numerically solves
 * \code
 * x**2 + b*x + c = 0
 * \endcode
 * for x in a numerically stable way.
 * 
 * Retuns the number of solutions, which is usually either 0 or 2. The case 1 is extremely rare as numerical
 * comparison is done directly and no care is taken for roundoff effects.
 * 
 * The solutions will be written in x1 and x2. In case of no solutions, both are set to NAN, in case of
 * one solution, both will have the same value. In case of two solutions x1 < x2.
 */
int roots_quad(double & x1, double & x2, double b, double c);

/** \brief The lngamma function
 *
 * Forwards to the boost implementation which is thread save (note that
 * C99 implementations need not be thread save).
 */
double lngamma(double x);

/** \brief copy a double vector
 *
 * Effectively as x = y, reading x and y as double vectors of length n.
 * Assumes that dest and src are aligned.
 */
inline void copy_fast(double * dest, const double * src, size_t n){
#ifndef __SSE2__
   std::copy(src, src + n, dest);
#else
   n += (n % 2);
   // 8 doubles at once in an unrolled loop (64 bytes = L1 cache line size ...)
   size_t imax = n / 8;
   for(size_t i=0; i < imax; ++i){
	   __m128d XMM0 = _mm_load_pd(src);
	   _mm_store_pd(dest, XMM0);
	   src += 2;
	   dest += 2;
	   __m128d XMM1 = _mm_load_pd(src);
	   _mm_store_pd(dest, XMM1);
	   src += 2;
	   dest += 2;
	   __m128d XMM2 = _mm_load_pd(src);
	   _mm_store_pd(dest, XMM2);
	   src += 2;
	   dest += 2;
	   __m128d XMM3 = _mm_load_pd(src);
	   _mm_store_pd(dest, XMM3);
	   src += 2;
	   dest += 2;
   }
   // rest:
   imax = (n % 8) / 2;
   for(size_t i=0; i<imax; ++i){
	   __m128d XMM0 = _mm_load_pd(src);
	   _mm_store_pd(dest, XMM0);
	   src += 2;
	   dest += 2;
   }
#endif
}

/** \brief add 2 vectors, possibly with sse optimization
 *
 * Calculates x+=y. Requires x and y to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 */
inline void add_fast(double * x, const double * y, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] += y[i];
   }
#else
  for(size_t i=0; i<n; i+=2){
     __m128d XMM0 = _mm_load_pd(x + i);
     __m128d XMM1 = _mm_load_pd(y + i);
     XMM0 = _mm_add_pd(XMM0, XMM1);
     _mm_store_pd(x+i  , XMM0);
  }
#endif
}

/** \brief multiply a vector with a constant, possibly with sse optimization
 *
 * Calculates x*=cy. Requires x to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 *
 */
inline void mul_fast(double * x, double c, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] *= c;
   }
#else
  __m128d XMM2 = _mm_set1_pd(c);
  for(size_t i=0; i<n; i+=2){
     __m128d XMM0 = _mm_load_pd(x + i);
     XMM0 = _mm_mul_pd(XMM0, XMM2);
     _mm_store_pd(x+i, XMM0);
  }
#endif
}


/** \brief add 2 vectors, possibly with sse optimization
 *
 * Calculates x+=c * y. Requires x and y to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 */
inline void add_fast_with_coeff(double * x, const double * y, double c, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] += c * y[i];
   }
#else
  __m128d XMM2 = _mm_set1_pd(c);
  for(size_t i=0; i<n; i+=2){
     __m128d XMM0 = _mm_load_pd(x + i);
     __m128d XMM1 = _mm_load_pd(y + i);
     XMM1 = _mm_mul_pd(XMM1, XMM2);
     XMM0 = _mm_add_pd(XMM0, XMM1);
     _mm_store_pd(x+i, XMM0);
  }
#endif
}

/** \brief add 2 vectors, possibly with sse optimization
 *
 * Calculates x+=c1 * y1 + c2 * y2. Requires x and y to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 */
inline void add_fast_with_coeff2(double * x, const double * y1, double c1, const double * y2, double c2, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] += c1 * y1[i] + c2 * y2[i];
   }
#else
  __m128d Xc1 = _mm_set1_pd(c1);
  __m128d Xc2 = _mm_set1_pd(c2);
  for(size_t i=0; i<n; i+=2){
     __m128d Xx = _mm_load_pd(x + i);
     __m128d Xy1 = _mm_load_pd(y1 + i);
     Xy1 = _mm_mul_pd(Xy1, Xc1);
     __m128d Xy2 = _mm_load_pd(y2 + i);
     Xy2 = _mm_mul_pd(Xy2, Xc2);
     Xx = _mm_add_pd(Xx, Xy1);
     Xx = _mm_add_pd(Xx, Xy2);
     _mm_store_pd(x+i, Xx);
  }
#endif
}

/** \brief possible redirections of log
 *
 * Tests have shown that the common log function is slow. In order
 * to allow easy switching for testing, all code in %theta should
 * use this log function.
 */
inline double log(double x){
#ifdef USE_AMDLIBM
    return amd_log(x);
#else
    return ::log(x);
#endif
}

inline double exp(double x){
#ifdef USE_AMDLIBM
    return amd_exp(x);
#else
    return ::exp(x);
#endif
}


}}

#endif
