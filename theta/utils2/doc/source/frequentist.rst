*****************************
Toy-based Frequentist Methods
*****************************

This page summarises those statistical methods which use frequentist concepts and toy datasets
as part of the evaluation. The most important one is p-value evaluation using the tail distribution of a likelihood ratio test statistic
for the null hypothesis. The methods documented here assume that you want to test the null hypothesis ``beta_signal = 0.0`` versus the alternative
``beta_signal > 0.0``.

.. note:: Toy-based methods often have large running times for more complicated models. In some cases, you can use asymptotic methods instead, see :ref:`likelihood`.

.. module:: theta_auto

The methods documented here are defined in ``theta_auto/frequentist.py``.

To calculate p-values from toys, there are three different "levels" of the interface where you should stick to one. From "higher" to lower level:

1. :meth:`theta_auto.discovery` - this method performs toys adaptively until a certain (user-specified) accuracy on the Z value is reached. It will run :program:`theta` locally,
   so this method is suitable if working interactively. However, it will take very long if the significance to be determined is large or if the calculation
   of the test statistic takes long.
2. If you need to toss a large number of background toys, you might want to run :program:`theta` in parallel, e.g. on a cluster.
   You can use :meth:`theta_auto.pvalue_bkgtoys_runs` to prepare :program:`theta` configuration files, on which you can
   run :program:`theta` distributed in a cluster. You can then call :meth:`theta_auto.pvalue` to calculate the actual p-value using
   the background-only toys. For details, refer to the :meth:`theta_auto.pvalue_bkgtoys_runs`, :meth:`theta_auto.Run`, and
   :ref:`distributed_running`.
3. The lowest level methods are :meth:`deltanll` and :meth:`dernll` which calculates the likelihood ratio and log-likelihood derivative, respectively for an arbitrary input dataset.
   Alterantively, they return the corresponding :class:`theta_auto.Run` instance. You can use this if you need full control about how the toys are performed.


High-level p-values calculation
===============================

.. autofunction:: discovery


Toy generation
==============


.. autofunction:: pvalue_bkgtoys_runs

.. autofunction:: pvalue


Test statistic Routines
=======================

.. autofunction:: deltanll

.. autofunction:: derll

