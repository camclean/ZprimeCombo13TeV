#!/usr/bin/env python
# -*- coding: utf-8 -*-

execfile("../lib.py")

import scipy.special, math

# returns the p-value for a hypothesis test of the "on/off" problem,
# where one makes a 2-channel measurement in a signal-free (off) and a signal channel (on).
# The expected number of events follows a poisson for both channels with expected
# (=modeled) mean values
#
# mu^model_on = mu_B + mu_S
# mu^model_off = mu_off
#
# where    mu_off / mu_B =: tau    is known precisely. The null hypothesis is mu_S = 0 in which case
# the model has one free parameter left; the latrenative hypothesis is the (compound) hypothesis
# mu_S > 0.
#
# More discussion in arxiV 0702156v4. Notation is also taken from there.
# In particular, see Formula (14) there.
def p_Bi(n_on, n_off, tau):
    res = scipy.special.betainc(n_on, 1+n_off, 1.0/(1 + tau))
    return res

def p_to_Z(p):
    return math.sqrt(2) * scipy.special.erfinv(1 - 2*p)


def test_expect(n_on, n_off, tau, expected_Z):
    observed_Z = p_to_Z(p_Bi(n_on, n_off, tau))
    #numbers in the papar are given to two significant digits, so the deviation
    # of the 100-fold difference should be .5 at most, assuming the figures have been
    # rounded correctly.
    if abs((observed_Z - expected_Z) * 100) > 0.5:
        print "FAIL: Z_Bi: got Z = %f, expected %f for n_on, n_off, tau = %f %f %f"  % (observed_Z, expected_Z, n_on, n_off, tau)
    else:
        print "PASS: Z_Bi: got Z = %f; expected Z = %f" % (observed_Z, expected_Z)

#test the values from the above mentioned paper:
test_expect(4, 5, 5.0, 1.66)
test_expect(6, 18.78, 14.44, 2.63)
test_expect(9, 17.83, 4.69, 1.82)
test_expect(17, 40.11, 10.56, 4.46)
test_expect(50, 55, 2.0, 2.93)
test_expect(67, 15, 0.5, 2.89)
test_expect(200, 10, 0.1, 2.20)
test_expect(523, 2327, 5.99, 5.93)
test_expect(498426, 493434, 1.0, 5.01)
test_expect(2119449, 23650096, 11.21, 6.40)
