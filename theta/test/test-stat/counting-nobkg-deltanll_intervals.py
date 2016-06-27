#!/usr/bin/env python
# -*- coding: utf-8 -*-

execfile("../lib.py")

import scipy.special, scipy.optimize, math, os

#return the interval of a counting experiment according to the delta-nll method
def interval_delta_nll(n, confidence):
    delta_nll = scipy.special.erfinv(confidence)**2
    if n==0.0: return (0, delta_nll)
    nll = lambda theta: n * math.log(1.0*n/theta) - n + theta - delta_nll
    interval_low = n
    while nll(interval_low) <= delta_nll: interval_low /= 2
    interval_high = n
    while nll(interval_high) <= delta_nll: interval_high *= 2
    interval_low = scipy.optimize.bisect(nll, interval_low, n)
    interval_high = scipy.optimize.bisect(nll, n, interval_high)
    return interval_low, interval_high
    
confidence=(0.6827, 0.9545)
for Theta in (5.0, 10000.0):
    execute_checked("sed \"s/__THETA__/%s/g\" counting-nobkg.cfg.tpl > counting-nobkg.cfg" % Theta)
    execute_checked("sed -i \"s/__THETA_WIDTH__/%s/g\" counting-nobkg.cfg" % (Theta * 0.1))
    exec_theta("counting-nobkg-deltanll_intervals.cfg")
    rows = sql("counting-nobkg-deltanll.db", "select writer__n_events_o, int__lower06827, int__upper06827, int__lower09545, int__upper09545 from products")
    for row in rows:
        n_events=row[0]
        intervals_expected = reduce(lambda x, y: x+y, map(lambda c: interval_delta_nll(n_events, c), confidence), ())
        for i in range(4):
            error=abs(intervals_expected[i] - row[i+1]) / float(Theta)
            if error > 0.001: fail("Theta = %s: n_events= %f; %d; expected %f, got %f" % (Theta, n_events, i, intervals_expected[i], row[i+1]) )
        #do additional checks in the asymptotic case:
        if Theta == "10000.0":
            factors = (1, 1, 4, 4)
            for i in range(1, 5):
                f = factors[i-1]
                dev = ((row[i] - n_events)**2 - f * n_events) / (f * n_events)
                if dev > 0.02: fail("Theta = %s; asymptotic check; n_events = %g" % (Theta, n_events))
    passed("Theta = %s " % Theta)
