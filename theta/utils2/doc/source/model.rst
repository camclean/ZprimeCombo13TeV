.. _model:

*********************
The Statistical Model
*********************

As discussed in :ref:`model_intro` in more detail, the statistical model specifies the predicted event yields (Poisson means) in each channel
as a function of the model parameters. This section summarizes the :class:`Model` class and all its building blocks, as well as the factory functions
usually used to build new :class:`Model`\ s.


.. _what_is_signal:

Specifying what is 'signal'
===========================

The Model class and statistical methods in theta-auto are designed to run them for situations in which there is more than one signal hypothesis, e.g.,
a Zprime with different masses. In this case, the background model is always the same, and only the histogram for signal changes. In this
case, one often wants to apply a method (such as setting limits) to all signal histograms in turn.

Note that this concept is specific to ``theta-auto`` not present in ``theta``\ . This means that each signal scenario will lead
to another theta configuration file and ``theta`` execution.

For example, one might have the signal processes 'zp1000', 'zp2000', and 'zp3000' for Zprime models
with different masses. In order for theta-auto to be able to correctly construct the background-only model, it must know what is signal. Use a call like ::

  model.set_signal_processes('zp*')
  
to declare which processes exactly are the signal processes. Everything else is considered the background-only model. To
calculate Bayesin upper limits on all of them in turn, call::

  bayesian_limits(model)
  
In this example, it will execute ``theta`` 6 times: twice for each signal (once for the observed limit and once for the expected limit).

Sometimes it makes sense to split the signal into several components which are conceptually the same process of interest
(e.g., if the process was generated according to production process or decay or if the uncertainties for different signal components are different). For example,
one might have split the production and processing of the Zprime signal sample according to production process into 'ggzp1000' and 'qqzp1000', etc., and in the example above
with three masses one would have six signal processes in total.
In such a case, one can specify more than one process name, all of which are scaled with the cross section scaling parameter "beta_signal". Such a group of processes to be
used together at the same time is called "signal process group" and is specified via a list of processes which make up the group, i.e., a list of strings.
For different signal hypotheses, such as different masses,
one has several such signal process groups. Therefore, a complete specification for which processes to run the statistical method on is a list of
signal process groups. The python type to do so is a dictionary which uses a user-defined string as a key (called the "signal process group id"), and the list of processes
which make up the signal process group as value. In the example above, a possible choice would be to use::

    signal_process_groups = {'zp1000': ['qqzp1000', 'ggzp1000'], 'zp2000': ['qqzp2000', 'ggzp2000'], 'zp3000': ['qqzp3000', 'ggzp3000']}

and the signal process group ids implicitly defined in this example are 'zp1000', 'zp2000' and 'zp3000'. These ids are used internally to construct file names and also
in the return value of the method. They should not contain any spaces other other special characters.


.. module:: theta_auto

  
Building Models
===============

While it is possible to build :class:`Model`\ s directly using constructors and the methods in :class:`Model`\ , you will usually
want to use one of the two methods documented here.

.. autofunction:: build_model_from_rootfile
.. autofunction:: theta_auto.higgs_datacard.build_model



Inspecting the Model
====================

After building the model, it is usually a good idea to check that the model actually contains the data you expect. For
this purpose, the ``model_summary`` method can be used.

.. autofunction:: model_summary




Evaluating the model prediction
===============================

For displaying the result of a fit, it is useful to evaluate the model prediction at the fitted value, or -- more generally -- for
a given set of model parameters. This can be done with the function :meth:`evaluate_prediction`. See :ref:`examples_mle_evaluate` for a
complete example of how to use this function.

.. autofunction:: evaluate_prediction





The Model class
===============

.. autoclass:: Model
  :members:
  

The Histogram class
===================

.. autofunction:: write_histograms_to_rootfile

.. autoclass:: Histogram
   :members:

The HistogramFunction class
===========================


.. autoclass:: HistogramFunction
  :members:
