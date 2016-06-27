#!/usr/bin/env python
# -*- coding: utf-8 -*-

# on/off problem, profile likelihood, see arXiv 0702156v4, section 5, formula (20).

execfile("../lib.py")

import math

#tuples n_on, n_off, tau, expected_Z to consider. From arXiv paper above (table 1)
tuples = [
    (4, 5, 5.0, 1.95),
    (6, 18.78, 14.44, 2.82),
    (9, 17.83, 4.69, 1.99),
    (17, 40.11, 10.56, 4.57),
    (50, 55, 2.0, 3.02),
    (67, 15, 0.5, 3.04),
    (200, 10, 0.1, 2.38),
    (523, 2327, 5.99, 5.95),
    (498426, 493434, 1.0, 5.01),
    (2119449, 23650096, 11.21, 6.40)
   ]

for t in tuples:
    n_on, n_off, tau, expected_Z = t
    execute_checked("sed \"s/__N_ON__/%e/g\" onoff.cfg.tpl > onoff.cfg" % n_on)
    execute_checked("sed -i \"s/__N_OFF__/%e/g\" onoff.cfg" % n_off)
    execute_checked("sed -i \"s/__TAU__/%e/g\" onoff.cfg" % tau)
    exec_theta("onoff_pl.cfg")
    query = "select hypotest__nll_sb - hypotest__nll_b from products;"
    rows = sql("onoff_pl.db", query)
    log_lr = rows[0][0]
    got_Z = math.sqrt(abs(2 * log_lr))
    if abs(got_Z - expected_Z) * 100 > 0.5:
        fail("Z off too much: got %f, expected %f" % (got_Z, expected_Z))
    else:
        passed("n_on, n_off = %f, %f: got %f, expected %f" % (n_on, n_off, got_Z, expected_Z))
