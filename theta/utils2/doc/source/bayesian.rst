.. _bayesian:

****************
Bayesian Methods
****************

.. module:: theta_auto

Bayesian methods construct the posterior -- given by the likelihood function times the prior -- and make statistical inferences from this posterior. This
involves "integrating out" the nuisance parameters which is done in theta using a Markov Chain Monte-Carlo method.

The provided functions in this module can be used to:

 * Calculate expected and observed Bayesian upper limits on the signal strength parameter ``beta_signal``. For this, use the :func:`bayesian_limits` function which evaluates both
   expected and observed limit calculation and reports the expected 1sigma and 2sigma bands and the observed limit in a form suitable for plotting.
 * Calculate Quantiles of the marginal posterior of ``beta_signal`` directly using :func:`bayesian_quantiles`. This is used internally by :func:`bayesian_limits` and allows more control,
   e.g. specifying which quantile you want, and on which input dataset to run on. This can be used to derive quantile-based Bayesian confidence intervals (where upper limits are just
   a special case). This method is also the right entry point for distributed computing, as it gives access to the individual theta configuration files by setting ``run_theta = False``.
 * Derive the marginal posterior density in a parameter (typically ``beta_signal``). This can be done with the :func:`bayesian_posteriors` method.
 * Calculate the posterior model prediction, i.e. mean and standard deviation of the Poisson mean in each bin. This can be used to visualize how well the model fit the data.
 
As all of these methods make use of Markov-Chain Monte-Carlo, all the options in the "[mcmc]" section apply, discussed in :ref:`options_class`.


Expected and observed limits
============================
.. autofunction:: bayesian_limits


Posterior quantiles
===================
.. autofunction:: bayesian_quantiles


Marginal posteriors
===================
.. autofunction:: bayesian_posteriors


Posterior model prediction
==========================
.. autofunction:: bayesian_posterior_model_prediction