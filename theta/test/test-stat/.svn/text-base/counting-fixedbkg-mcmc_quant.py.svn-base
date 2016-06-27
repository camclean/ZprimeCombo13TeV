#!/usr/bin/env python
# -*- coding: utf-8 -*-

execfile("../lib.py")

import scipy.special

#returns the posterior quantile of a likelihood poisson for n
# observed events
def true_q(n, theta_est, mu):
    #we want
    # result := 1 - scipy.special.gammaincc(n+1, theta_est + mu) / scipy.special.gammaincc(n+1, mu)
    # this is numerically unstable if approximately 1 - 1 is calculated. The second term is
    # getting to 1 soon if   theta_est + mu < n+1.
    #
    # If setting
    # A:= gamma(n+1, mu)
    # B:= gamma(n+1, mu + theta_est)
    # where gamma is the *lower* incomplete gamma function, result can be rewritten as
    #   result = B + (1 + B) * (sum_n=1^inf  A^n)
    # A is very small in this case such that this coverges rapidly.
    result = 1 - scipy.special.gammaincc(n+1, theta_est + mu) / scipy.special.gammaincc(n+1, mu)
    if result > 1 - 1e-5:
        print "result close to 1: ", result
        A = gammainc(n+1, mu)
        B = gammainc(n+1, mu)
        result = B
        sum = A
        while A > 0:
            A *= A
            sum += A
        result += (1+B)*sum
        print "new result: ", result
    return result

#The meaning of the _old functions is as follows:
# to a first approximation, one might be tempted to subtract the fixed background mu from
# the number of observed events, n and then calculate within a model without background. However, this
# is not correct. To show that here, the deviation will be displayed.

def true_q_old(n, theta_est):
    return scipy.special.gammainc(n+1, theta_est)

quantiles=(0.5,0.16,0.84,0.95)
for mu in (2000.0, 3.0):
    for Theta in (5.0,10000.0):
        execute_checked("sed \"s/__THETA__/%s/g\" counting-fixedbkg.cfg.tpl > counting-fixedbkg.cfg" % Theta)
        execute_checked("sed -i \"s/__THETA_WIDTH__/%s/g\" counting-fixedbkg.cfg" % (Theta * 0.1))
        execute_checked("sed -i \"s/__MU__/%s/g\" counting-fixedbkg.cfg" % mu)
        exec_theta("counting-fixedbkg-mcmc_quant.cfg")
        rows = sql("counting-fixedbkg-mcmc.db", "select writer__n_events_o, mcmc__quant05000, mcmc__quant01600, mcmc__quant08400, mcmc__quant09500 from products")
        events_low=0
        events_hi=0
        events_low_old=0
        events_hi_old=0
        for row in rows:
            n_events=row[0]
            quantiles_est = map(lambda theta_est: true_q(n_events, theta_est, mu), row[1:])
            quantiles_est_old = map(lambda theta_est: true_q_old(n_events - mu, theta_est), row[1:])
            differences = map(lambda (x,y): x-y, zip(quantiles, quantiles_est))
            differences_old = map(lambda (x,y): x-y, zip(quantiles, quantiles_est_old))
            if min(differences) < -0.015: events_low += 1
            if max(differences) > 0.015: events_hi += 1
            if min(differences_old) < -0.015: events_low_old += 1
            if max(differences_old) > 0.015: events_hi_old += 1
        if events_low < 10 and events_hi < 10:
            passed("Theta = %f; mu = %f; events low / high: %d %d; events low / high old: %d %d" % (Theta, mu, events_low, events_hi, events_low_old, events_hi_old))
        else:
            fail("Theta = %f; ; mu = %f; events low / high: %d %d; events low / high old: %d %d" % (Theta, mu, events_low, events_hi, events_low_old, events_hi_old))
