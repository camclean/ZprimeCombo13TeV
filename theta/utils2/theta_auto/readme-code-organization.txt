The code has mainly two parts: 1. model building, and 2. statistical inference. There is a third, much smaller part
of utility functions.


1. Model Building
------------------
The "model building" part is concerned
with manipulating (creating, modifying, etc.) instances of Model and related classes (HistogramFunction, Distribution,
Function).

The model building interface is defined in
 * Model.py
 * Function.py


2. Statistical Inference
------------------------
The second part is the "statistical inference" part which consists of methods which the user calls to make statistical inference. This is
the code which actually calls theta. This code is structured in two layers: (i) the low-level interface which creates the config files, calls theta,
reads in the theta result, etc. and (ii) the "high-level" part which provides the user a simple interface for many common tasks.

2. (i) The low-level interface
------------------------------
is the one which
 * creates theta config files
 * calls theta (if requested)
 * returns the output data produced by theta, in a low-level format
In general, the low-level interface reflects the structure of theta modules and plugins.

All low-level methods are in the file theta_interface2.py. The main ingredients of this interface are
classes for the "Main" theta run as well as producers which can be used.

The methods of the low-level interface are meant to be as independent as possible from the exact use. For example, they should
avoid unneeded dependencies on plotting packages, should not write anything to the global "Report" object, etc.


2. (ii) The high-level interface
--------------------------------
contains the methods the user usually uses. It does not execute theta directly
but always goes through the low-level interfac. Otherwise, high-level interfaces can do anything they want in particular call
low-level methods as many times as required by the method, write to the "Report" object, make plots, write to standard out, etc.

The high-level methods are grouped by statistical method and are in:
 * bayesian.py for bayesian methods
 * cls_limits.py for the CLs method, asymptotic_cls.py for the asymptotic CLs method
 * neyman.py for the Neyman construction, including Feldman-Cousins intervals
 * likelihood.py for likelihood methods such as mle, profile likelihood intervals, approximate likelihood based, etc.
 * frequentist.py for frequentist inference, such as p-values from toys, re-interpreting a (Bayesian) model frequentistically, etc.
 
In general, the methods 
 * take a parameter "signal_process groups" which is a dictionary spid -> list of signal process names; the default
   is to use model.signal_process_groups. In general, it is only valid to specify a subset of model.signal_process_groups,
   although the methods themselves do not necessarily check that
 * take the parameters "input", "nuisance_prior_toys", "nuisance_constraint", "signal_prior", "n", whenever appropriate
The default choice for all parameter (None) is to use the value from the model. There should be no default value for parameters
"input" and "n" where this is not well-defined.
   
 * write_report: a boolean indicating whether or not to write the results to the global report object, default True
 * options: an theta_interface2.Options object. The default is None and mean that the default is used (by calling the default-constructor)

The return value of these methods should be a (nested) dictionary with spid as (first) key. The value in this dictionary highly depends
on the method, so only very few more rules exist. If the result is reported on a "per-toy" basis, the next key should be
the name of the "result column", the value is a list of results, and the length of this list is the number of toys "n" (although can be less in case
there was an error for some toys). This way around -- a dictionary with lists instead of a list with dictionaries --
saves some space in case of large results.

 
There are also some methods which can be considered "high-level" not associated to a statistical method:
 * model_summary defined in model_summary.py
 
 
The high-level methods make use of
 * plotutil.py for plotting
 * Report.py for generating (html) reports
 

3. Utility functions
--------------------
A number of utility functions required in several places are defined in
 * utils.py

