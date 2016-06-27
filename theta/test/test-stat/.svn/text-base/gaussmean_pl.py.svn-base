#!/usr/bin/env python
# -*- coding: utf-8 -*-

# on/off problem, profile likelihood, see arXiv 0702156v4, section 5, formula (21).

execfile("../lib.py")

import math

#tuples muhat_b, sigma_b, n_on, expected_Z to consider. From arXiv paper above (table 1)
tuples = [
    (1.0, 0.447, 4, 2.00),
    (1.3, 0.3, 6, 2.83),
    (3.8, 0.9, 9, 2.02),
    (3.8, 0.6, 17, 4.62),
    (27.5, 3.71, 50, 3.10),
    (30.0, 7.75, 67, 3.45),
    (100.0, 31.6, 200, 2.90),
    (388.6, 8.1, 523, 5.96),
    (493434, 702.4, 498426, 5.02),
    (2109732, 433.8, 2119449, 6.40)
   ]

for t in tuples:
    muhat_b, sigma_b, n_on, expected_Z = t
    execute_checked("sed \"s/__MUHAT_B__/%e/g\" gaussmean_model.cfg.tpl > gaussmean_model.cfg" % muhat_b)
    execute_checked("sed -i \"s/__SIGMA_B__/%e/g\" gaussmean_model.cfg" % sigma_b)
    execute_checked("sed -i \"s/__N_ON__/%e/g\" gaussmean_model.cfg" % n_on)
    exec_theta("gaussmean_pl.cfg")
    query = "select hypotest__nll_sb - hypotest__nll_b from products;"
    rows = sql("gaussmean_pl.db", query)
    log_lr = rows[0][0]
    got_Z = math.sqrt(abs(2 * log_lr))
    if abs(got_Z - expected_Z) * 100 > 0.5:
        fail("Z off too much: got %f, expected %f" % (got_Z, expected_Z))
    else:
        passed("muhat_b, sigma_b = %f, %f: got %f, expected %f" % (muhat_b, sigma_b, got_Z, expected_Z))
