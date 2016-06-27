************
Installation
************


Prerequisites
=============

``theta`` should run on all modern linux-based operating systems, although due to limitations of developer
resources, it is only regularly tested on recent versions of Ubuntu, as well as scientific linux 5, using g++ (and sometimes
clang). There have also been reports of successful builds on Mac OS X, but no regular testing of newer versions is done.

For compiling ``theta``, you need:
  #. sqlite3
  #. boost (tested with 1.49)
  #. ROOT (tested with 5.28 - 5.32)
   
Be sure to install the header/"-dev" packages. For experimental llvm support, you need llvm-3.1 (although chances are you can compile it
with more recent llvm versions with small changes).

For running *theta-auto*, you need python 2.7 and in addition numpy, scipy, matplotlib, sqlite3, and ROOT bindings.

If you are CMS user, note that CMSSW comes with all necessary dependencies for compiling and running theta, so you don't
need to install any of them; see :ref:`compile_cmssw`.
   
Obtaining theta
================

``theta`` is available as source-code distribution via subversion only. The latest, recommended version can be obtained by running::

 svn co https://ekptrac.physik.uni-karlsruhe.de/public/theta/tags/testing theta

The testing tag is regularly updated and is the main tag suitable for most users.

The "theta" directory just created by the subversion checkout is referred to as the "main theta directory" below.

.. _compile_make:

Compiling theta using Makefile
==============================

You can compile theta by running "make" in the main theta directory. You can configure some details in the file Makefile.options,
see the comments there. The standard settings should be appropriate for most users.

.. _compile_cmssw:

Compiling theta using CMSSW
============================

It is possible to compile *theta*, using the packages CMSSW provides. This is built into the standard theta Makefiles.
You just have to set up CMSSW before you type "make" in the main theta directory.

.. _compile_cmake:

Compiling theta using cmake
===========================

For cmake, use::

 mkdir build
 cmake ..
 make
   
For available build options, have a look at the top of the top-level CMakeLists.txt



