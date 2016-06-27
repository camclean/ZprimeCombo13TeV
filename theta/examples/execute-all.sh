#!/bin/bash

# Execute all examples. Run this script from the main theta directory.
# This script might not be very useful to learn from the examples very much.
#
# Its main purpose is to ensure that all examples run without error prior to releasing theta.

fail()
{
  echo "FAIL: $*";
  exit 1;
}

execute_checked()
{
   ($*)
   if [ $? -gt 0 ]; then
      fail "executing $* returned error"
   fi
}
               
exec_theta()
{
   echo "Executing theta $*..."
   execute_checked bin/theta $*
}

[ -f bin/theta ] || fail did not find theta in bin/theta. Make sure to compile theta and execute this script from the main theta directory.

exec_theta examples/CountingExp.cfg
exec_theta examples/gaussoverflat.cfg
exec_theta examples/gaussoverflat-intervals2.cfg
exec_theta examples/gaussoverflat-intervals.cfg
exec_theta examples/gaussoverflat-mcmc.cfg
exec_theta examples/gaussoverflat-mcmc_histo.cfg
exec_theta examples/gaussoverflat-mcmc_quantiles.cfg
exec_theta examples/gaussoverflat-mle.cfg
exec_theta examples/gaussoverflat-nll_scan.cfg
exec_theta examples/gaussoverflat-rootfile.cfg
exec_theta examples/gaussoverflat-writer.cfg
exec_theta examples/gaussoverflat-cls.cfg
