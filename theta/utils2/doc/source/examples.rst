.. _examples:

********
Examples
********

This section contains many common examples. Most examples are shown using a simple counting experiment model constructed with the method
``test_model.simple_counting`` which has one signal process called ``s``, and an optional log-normal uncertainty on the background, modeled
via a nuisance parameter called ``"bunc"``. The model prediction for the Poisson mean is

.. math::

\mu = \beta_{\rm signal} \cdot s + e^{{\rm bunc}\cdot{\rm b_uncertainty}} \cdot b

where ``s``, `b`` and ``b_uncertainty`` are constants, passed as parameters to ``test_model.simple_counting``. The parameters of this model
are ``beta_signal`` and ``bunc``.


.. _examples_mle:

=======================
Maximum Likelihood Fits
=======================


.. _examples_mle_fit:

Perform a maximum likelihood fit
--------------------------------

See ``utils2/examples/mle/fit.py`` ::

 model = test_model.simple_counting(s = 10.0, n_obs = 12, b = 2.0)
 result = mle(model, input = 'data', n = 1)
 print result
 
This will print::

 {'s': {'__nll': [-17.818879797456002], 'beta_signal': [(0.9999999999999996, 0.34463876178582953)]}}
 
So the result is a python dictionary. As explained in :ref:`return_values`, the first-level key specifies
the signal group id, here "s". The dictionary `result['s']` has to keys: "__nll" contains the
negative log-likelihood values at the minimum. These values are not usually not useful on their own, but can be used to construct
likelihood ratios by making several maximum likelihood fit. The "beta_signal" contains a list of ``n = 1`` results, each
result is a two-tuple of the fitted parameter values and its uncertainty, so the result of the maximum likelihood fit
in this example can be summarized as

.. math::

 \hat{\beta}_{\rm signal} = 1.00 \pm 0.34


Making the fit more robust
--------------------------

In some cases, the fit might not converge. In this case, you can use the ``Options class`` (see :ref:`options_class`) for more control
about the minimizer settings. Also, in case you do not need the error, you can disable the error calculation.
So a robust maximum likelihood fit would be (see ``utils2/examples/mle/robust-fit.py``) ::

 model = test_model.bernstein_model(10) # a model of 10 bernstein polynomials, difficult to minimize
 mle(model, 'toys:0.0', 100) # will have many failures
 mle(model, 'toys:0.0', 100, with_error = False)
 options = Options()
 options.set('minimizer', 'strategy', 'robust')
 mle(model, 'toys:0.0', 100, with_error = False, options = options)

If you execute theta-auto on this python file, you'll note that the first theta execution will have many errors.
Specifying ``with_error = False`` for the ``mle`` method removes these errors already: this will not run MINUIT's migrad,
but you will not have an error estimate (if you need one, you can try using profile likelihood intervals via :meth:`theta_auto.pl_intervals` ).

The third call of ``mle`` in this example then uses different options for the minimizer.

.. warning:: If you find yourself in the situation where you the minimization fails, even the "robust" minimizer probably does not return
 the correct minimum. In this case, it might be a good idea to understand why the minimizer fails and
 to check that the result of the minimizer makes sense (e.g., by checking the pull distribution for toys).

.. _examples_mle_evaluate:
 
Get the histograms at the maximum likelihood fit
------------------------------------------------

For displaying the result of the fit, it can be useful to get the histograms of the model scaled and morphed according to the parameterization
in the model. This can be done with the :meth:`theta_auto.evaluate_prediction` method. The resulting histograms can either be inspected directly in the python
code. In many cases, you probably only want to write them to a root file to plot them with your external script. In this case, use the :meth:`theta_auto.write_histograms_to_rootfile` method.
The resulting histogram names have the same naming convention as the nominal histograms for :meth:`theta_auto.build_model_from_rootfile`, so they are called "<channel>__<process>".

A complete example which performs a background-only maximum likelihood fit and writes the morphed/scaled histograms for these parameter values to a root file
is (see ``utils2/examples/mle/evaluate.py``) ::

 model = test_model.simple_counting(s = 2.0, b = 10.0, n_obs = 7, b_uncertainty = 3.0)
 signal_process_groups = {'': []}
 fit = mle(model, 'data', n = 1, signal_process_groups = signal_process_groups)
 parameter_values = {}
 for p in model.get_parameters([]):
     parameter_values[p] = fit[''][p][0][0]
 histos = evaluate_prediction(model, parameter_values, include_signal = False)
 write_histograms_to_rootfile(histos, 'histos-mle.root')
 
The first line builds a (trivial) example model with only one observable called "obs", with one signal process called "s" and one background process called "b".
The second line defines the signal process groups to run on: here, we only want to
make a background-only fit, so an empty string is used as signal process id, and the corresponding list of signal processes is empty as well; see
:ref:`what_is_signal` for details. The third line runs the actual maximum -likelihood fit on data. The next three lines build the dictionary of parameter
values. The expression ``fit[''][p][0][0]`` might look a bit complicated and deserves some explanation: the keys to the return value of ``mle`` are
(i) the signal process id (the keys to ``signal_process_groups``), (ii) the parameter name, (iii) the
toy index (from 0 to ``n-1`` where ``n`` is the value passed to ``mle``). Using these three keys, the value is a two-tuple of ``(fitted value, uncertainty)``, so
the fourth and last index ``0`` selects the fitted value. See the documentation of :meth:`theta_auto.mle` for details.

The second-to-last line does the actual evaluation of the model prediction, and the last line writes the histograms to the root file "histos-mle.root".

The resulting root file will only have one histogram in this example, ``obs__b``, but this example can also be applied to more complicated models.


==========
CLs Limits
==========

Calculate the expected asymptotic CLs limits with injected signal 
-----------------------------------------------------------------

To calculate the expected asymptotic CLs limits with a non-zero signal, use the ``beta_signal_expected`` parameter of :meth:`asymptotic_cls_limits`::

  model = test_model.simple_counting(s = 2.0, b = 10.0, n_obs = 7, b_uncertainty = 3.0)
  exp, obs = asymptotic_cls_limits(model, beta_signal_expected = 1.0)
  print exp, obs

The resulting expected limits have been evaluated with the signal strength parameter ``beta_signal`` set to 1.0. This method, however, only works if
you want to inject the same signal that you want to exclude.

In general, you might want to inject a signal different from the one you fit. To calculate the expected limit in such a case, the toy data is generated according to
the "signal-injected" model in a first stage. In the second stage, this toy data is passed as ``input`` parameter to the asymptotic CLs limit calculation::

  model_toys = test_model.simple_counting(s = 1.2, b = 10.0, n_obs = 7, b_uncertainty = 3.0)
  model_toys = get_bootstrapped_model(model_toys)
  N_toy = 1000
  data = make_data(model_toys, 'toys:1.0', N_toy)
  
  model_limit = test_model.simple_counting(s = 2.0, b = 10.0, n_obs = 7, b_uncertainty = 3.0)
  exp, obs = asymptotic_cls_limits(model_limit, input = data['s'], n = N_toy)
  print obs
  
Some things to note:
 * the models for generating toy data and for deriving the CLs limits are completely decoupled. This is Ok as long as they define the same observables (including the same binning, etc.).
 * in the second line, the toy model is adapted from the default "prior-predictive" mode (in which each nuisance parameter has a Bayesian prior) to the "frequentist" view in which for each nuisance parameter, there is a (pseudo-)measurement. The central values are bootsrapped from data by making a maximum likelihood fit. This is required for the data generation to generate the observations of the nuisance parameter measurements consistently with what is done in ``asymptotic_cls_limits``.
 * the asymptotic CLs limits for the toy data are reported as the "observed" limits; in general, anything specified in the ``input`` parameter will contribute to the ``obs`` result. The ``exp`` variable is unchanged, i.e. it contains the asymptotic CLs limits for an ensemble of (toy-)data distributed according to ``model_limit`` with ``beta_signal = 0.0``.
 
See also: :meth:`theta_auto.make_data`, :meth:`theta_auto.asymptotic_cls_limits`.
