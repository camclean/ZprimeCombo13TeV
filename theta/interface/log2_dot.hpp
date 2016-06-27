#ifndef LOG2_DOT_HPP
#define LOG2_DOT_HPP

#include <cmath>

namespace theta{
    
/** \brief calculate the negative-log-likelihood for the data and prediction for the Poisson means
 *
 * calculates
 * sum_{i=0}^{n-1}  pred[i] - data[i] * log(pred[i])
 *
 * If, for some i, data[i] and pred[i] are both zero, this is i is skipped in
 * the sum (i.e., summand is treated as 0).
 *
 * If, for some i, data[i] > 0.0 and pred[i] <= 0.0, +infinity is returned.
 */
double template_nllikelihood(const double * data, const double * pred, unsigned int n);

/** \brief Same as template_nllikelihood, but add constant term to increase robustness
 * 
 * A constant term is added to the negative log-likelihood that should increase numerical accuracy
 * around the minimum.
 */
double template_nllikelihood_robust(const double * data, const double * pred, unsigned int n);


/** \brief Calculate a generalized chi2 test statistic from the given model prediction for the Poisson means and observation
 * 
 * The generalized (pseudo-)chi2 coincides with
 * 
 * chi2 = 0.5 * sum_{i=0}^{n}  (data[i] - pred[i]) / sqrt(pred[i])
 * 
 * if pred[i] is large. In general, it is twice the negative log-likelihood ratio, i.e.
 * 2 * (template_nllikelihood(data, pred, n) - template_nllikelihood(data, data, n)).
 */
double template_pchisquare(const double * data, const double * pred, unsigned int n);

}

#endif
