#include "interface/log2_dot.hpp"
#include "interface/utils.hpp"

#include <limits>
#include <math.h>

using namespace std;

double theta::template_nllikelihood(const double * data, const double * pred, unsigned int n){
   double result = 0.0;
   for(unsigned int i=0; i<n; ++i){
        result += pred[i];
        if(pred[i] > 0.0){
             if(data[i] > 0.0){
                 result -= data[i] * theta::utils::log(pred[i]);
             }
         }else if(data[i] > 0.0){
             return numeric_limits<double>::infinity();
         }
    }
    return result;
}

double theta::template_nllikelihood_robust(const double * data, const double * pred, unsigned int n){
   double result = 0.0;
   for(unsigned int i=0; i<n; ++i){
        result += pred[i] - data[i];
        if(pred[i] > 0.0){
             if(data[i] > 0.0){
                 result -= data[i] * theta::utils::log(pred[i] / data[i]);
             }
         }else if(data[i] > 0.0){
             return numeric_limits<double>::infinity();
         }
    }
    return result;
}

double theta::template_pchisquare(const double * data, const double * pred, unsigned int n){
    double result = 0.0;
    for(unsigned int  i=0; i<n; ++i){
        const double n = data[i];
        const double mu = pred[i];
        if(mu > 0){
            if(n > 0){
                result += n * utils::log(n / mu) + mu - n;
            }
            else{
                result += mu;
            }
        }
        else if(n > 0){
            result = numeric_limits<double>::infinity();
            break;
        }
    }
    result *= 2;
    return result;
}
