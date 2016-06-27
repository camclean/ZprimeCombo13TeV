#!/usr/bin/env python
# -*- coding: utf-8 -*-

execfile("../lib.py")

import scipy.special

#returns the posterior quantile of a likelihood poisson for n
# observed events
def true_q(n, theta_est):
    return scipy.special.gammainc(n+1, theta_est)

quantiles=(0.5,0.16,0.84,0.95)
for Theta in (5.0,10000.0):
    execute_checked("sed \"s/__THETA__/%s/g\" counting-nobkg.cfg.tpl > counting-nobkg.cfg" % Theta)
    execute_checked("sed -i \"s/__THETA_WIDTH__/%s/g\" counting-nobkg.cfg" % (Theta * 0.1))
    exec_theta("counting-nobkg-mcmc_quant.cfg")
    rows = sql("counting-nobkg-mcmc.db", "select writer__n_events_o, mcmc__quant05000, mcmc__quant01600, mcmc__quant08400, mcmc__quant09500 from products")
    events_low=0
    events_hi=0
    for row in rows:
        n_events=row[0]
        quantiles_est = map(lambda theta_est: true_q(n_events, theta_est), row[1:])
        differences = map(lambda (x,y): x-y, zip(quantiles, quantiles_est))
        if min(differences) < -0.015: events_low += 1
        if max(differences) > 0.015: events_hi += 1
    if events_low < 15 and events_hi < 15:
        passed("Theta = %s; events low / high: %d %d" % (Theta, events_low, events_hi))
    else:
        fail("Theta = %s; events low / high: %d %d" % (Theta, events_low, events_hi))
