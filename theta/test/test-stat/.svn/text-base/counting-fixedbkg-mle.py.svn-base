#!/usr/bin/env python
# -*- coding: utf-8 -*-

execfile("../lib.py")

for mu in (2000.0, 3.0):
    for Theta in (5.0, 8000.0):
        execute_checked("sed \"s/__THETA__/%s/g\" counting-fixedbkg.cfg.tpl > counting-fixedbkg.cfg" % Theta)
        execute_checked("sed -i \"s/__MU__/%s/g\" counting-fixedbkg.cfg" % mu)
        exec_theta("counting-fixedbkg-mle.cfg")
        #select the failing ones, i.e., pseudo experiments for which
        # n > mu,  but  \hat\Theta != n - \mu   or
        # n <= mu, but \hat\Theta != 0
        #where != means that the difference is larger than the 1e-3 times the stat. error
        stat_error = math.sqrt(Theta + mu)
        query = "select count(*) from products where (writer__n_events_o > %(mu)g and abs(%(mu)g + mle__Theta - writer__n_events_o) > %(tol)g) or "\
                   "(writer__n_events_o <= %(mu)g and abs(mle__Theta) > %(tol)g)" % {'mu': mu, 'tol': 1e-2 * stat_error} 
        rows = sql("counting-fixedbkg-mle.db", query)
        num_fail = rows[0][0]
        if num_fail > 0: fail("maximum likelihood estimate off too much in %d cases (SQL: %s)" % (num_fail, query))
        else: passed("Theta = %f; mu = %f" % (Theta, mu))
        if Theta > 5000:
            rows = sql("counting-fixedbkg-mle.db", "select count(*) from products where abs(mle__Theta_error * mle__Theta_error - mle__Theta - %(mu)g) > %(tol)g" % {'mu': mu,
                'tol': 1e-2 * stat_error})
            num_fail = rows[0][0]
            if num_fail > 0: fail("asymptotic error estimate in mle off too much in %d cases" % num_fail)
            else: passed("Theta = %f asymptotic error" % Theta)
