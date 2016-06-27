*************************************************
Low-level theta interface and distributed running
*************************************************

The previous sections focussed on using higher-level methods of theta-auto which usually reqwuire little or no knowledge about the
architecture of :program:`theta` itself. In some special cases, you might want to use a lower-level interface to have more control about
the .cfg files created. This low-level interface is partly described here, although you might need to read the source code to fully understand
the details.

The first section discusses how to run theta distributed on a cluster.

.. _distributed_running:

Distributed Running
====================

Before explaining how to actually run :program:`theta` distributed, here is a outline of what happens in theta-auto "behind the scenes" if you call a statistical method:
 * One or more (text-based, human readable) self-contained configuration files with suffix ".cfg" are created in the workdir.
   The file name is based on the statistical method run, the value of the "input" parameter, the signal process
   group id, and a hash based on the content of the configuration file. Each configuration file is set up such that :program:`theta`
   creates a .db file of the same name where the result is saved.
 * The "cache" subdirectory of the analysis workdir is checked whether it contains the config file just created (by comparing the name and the content)
   and the according output .db file of the same name. If this is the case,
   the .db file is used directly and :program:`theta` is not executed.
   Otherwise, :program:`theta` is executed locally on this .cfg file. In case :program:`theta` was successful, the .cfg file and the resulting .db file
   are copied to the cache directory, for later re-use.
   
So if running theta-auto twice on the same analysis.py file, theta-auto will notice in the first step that the .db file already exists and not execute :program:`theta` again.

Running :program:`theta` distributed makes use of the fact that the .cfg files created in the first step are self-contained. Typically you would:

1. call the python methods to create the .cfg files, but prevent the execution of theta. How this is done depends on the method, but some of the more time-consuming
   methods have the option ``run_theta`` which you can set to ``False``, e.g. :meth:`~theta_auto.bayesian_quantiles` and the test statistics methods :meth:`~theta_auto.deltanll` and
   :meth:`~theta_auto.derll`.
   If setting ``run_theta`` to ``False``, the return value of those methods are instances of :class:`~theta_auto.theta_interface.Run` which
   you can use to generate theta .cfg files.
2. use the "gridpack" directory to build a self-contained theta binary .tar.gz file (see also the README file in this directory)
3. submit jobs to a cluster which use the gridpack and execute "theta <cfg file>" on the worker node. Make sure to copy the resulting .db file back.
4. copy the .db output file *and* the .cfg file file to the "cache" subdirectory in the analysis workdir
5. execute theta-auto again, this time setting "run_theta" to True. theta-auto should notice that the file in the cache is already present in the cache directory and will use the
   output .db file instead of running :program:`theta` locally
   
While this procedure involves quite some steps which have to be done manually (or set up manually), it allows full control
over the process and should therefore enable you to make it work with many different distributed systems.


.. _run_object:

Run object
==========

A ``Run`` instance corresponds to one :program:`theta` configuration file, to one execution of :program:`theta` and to one resulting .db file. Interacting
with ``Run`` instances directly is usually not necessary, unless you want to have detailed control about how/where :program:`theta` is executed, e.g. for running
:program:`theta` on the cluster, as outlined in :ref:`distributed_running`.


.. module:: theta_auto.theta_interface

.. autoclass:: Run
  :members: __init__,get_configfile,run_theta,wait_for_result_available,get_db_fname,get_products

