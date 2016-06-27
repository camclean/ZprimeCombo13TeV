.. _likelihood:

************************
Likelihood-based methods
************************

This page summarises those statistical methods which exploit asymptotic
properties of the likelihood function, most importantly the maximum likelihood method
for point estimation, the profile likeilhood method for interval estimation, and the approximate Z-value (significance)
determination.

The methods documented here are defined in ``theta_auto/likelihood.py``.


.. module:: theta_auto


Maximum likelihood estimate
============================

.. autofunction:: mle


Profile Likelihood method
==========================

.. autofunction:: pl_interval(model, input, n, cls = [cl_1sigma, cl_2sigma], signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, parameter = 'beta_signal')

.. autofunction:: pl_interval_coveragetest(model, spid, beta_signal_values = [0.1*i for i in range(21)], n = 1000, cl = cl_1sigma, nuisance_prior_toys = None, nuisance_constraint = None, options = None)

.. autofunction:: nll_scan


Approximate/asymptotic Z-value
==============================

.. autofunction:: zvalue_approx
