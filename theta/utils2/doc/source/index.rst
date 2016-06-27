.. theta-auto documentation master file, created by
   sphinx-quickstart on Thu Jun 21 14:25:10 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

theta
=====

``theta`` is a framework for template-based statistical modeling and inference, focussing on problems
in high-energy physics. It provides the possibility for the user to express a "model", i.e.,
the expected data distribution, as function of physical parameters. This model can be used
to make statistical inference about the physical parameter of interest. Modeling is "template-based"
in the sense that the expected data distribution is always expressed as a sum of histograms.

.. _doc_overview:

Documentation overview
======================

There are several sources of documentation for ``theta``, listed with decreasing relevance:
 * The documentation on these pages cover *theta-auto*, the python interface to theta, which should be the right starting place for most users.
 * A general, physics-focused introduction is available as pdf `here <http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/theta.pdf>`_. It includes some documentation of the cfg-file layer, but does not mention *theta-auto* yet.
 * For a documentation of the C++ code (including the configuration file format), refer to the `doxygen documentation <http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/testing/html/index.html>`_.
 * The source code is also available online via `trac <https://ekptrac.physik.uni-karlsruhe.de/trac/theta/>`_.


Contents
========

.. toctree::
   :maxdepth: 2

   installation
   intro
   model
   likelihood
   frequentist
   cls
   bayesian
   utility
   common-parameters
   examples
   low-level

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

