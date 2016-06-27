#include <boost/test/unit_test.hpp>
#include <iostream>

#include <stdint.h>
#include <math.h>

#include "interface/log2_dot.hpp"
#include "interface/histogram.hpp"
#include "test/utils.hpp"

using namespace std;
using namespace theta;

#ifdef USE_AMDLIBM
#include "amdlibm/include/amdlibm.h"
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

namespace {

double template_nllikelihood_reference(const double * data, const double * pred, unsigned int n){
   double result = 0.0;
   double pred_norm = 0.0;
   for(unsigned int i=0; i<n; ++i){
        pred_norm += pred[i];
        if(pred[i] > 0.0){
             if(data[i] > 0.0){
                 result -= data[i] * log(pred[i]);
             }
         }else if(data[i] > 0.0){
             return std::numeric_limits<double>::infinity();
         }
    }
    return result + pred_norm;
}


}


BOOST_AUTO_TEST_SUITE(log2_dot_tests)


BOOST_AUTO_TEST_CASE(tl){
    unsigned const int N=5;
    double data[N];
    double pred[N];

    for(int i=0; i<5; i++){
       data[i] = i;
       pred[i] = i + 0.2;
    }
    double ref = template_nllikelihood_reference(data, pred, N);
    double res = template_nllikelihood(data, pred, N);
    BOOST_CHECK(close_to(ref, res, 1.0));
    
    //check that inf is returned if pred is 0.0:
    pred[N/2] = 0.0;
    res = template_nllikelihood(data, pred, N);
    BOOST_CHECK(std::isinf(res) && res > 0);
    
    //check that bin is skipped if data is 0.0:
    data[N/2] = 0.0;
    ref = template_nllikelihood_reference(data, pred, N);
    res = template_nllikelihood(data, pred, N);
    BOOST_REQUIRE(!std::isinf(res) && !std::isinf(ref));
    BOOST_CHECK(close_to(ref, res, 1.0));
}



BOOST_AUTO_TEST_SUITE_END()
