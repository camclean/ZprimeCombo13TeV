#!/bin/bash

#
# WARNING: running this script takes quite long (about 60 minutes) as it produces *all* data and plots of the paper.
#
# This script executes all steps necessary to produce the plots shown in the paper.
# It must be invoked from the directory it resides in and needs:
# * that theta was compiled, including the root plugins and the sqlite plugin
# * a working setup of ROOT, in order to compile scripts/make_histos.cc
# * python with matplotlib, numpy, scipy and sqlite3 for the plotting / analysis routines in scripts/plot_results.py
#
# Note that the numbers might not match 100% with those given in the papaer as all results use
# pseudo random numbers with a time-based seed.
#
# The configuration files for theta reside in the "cfg" directory. Running these will write the results to a file
# with ending of .db of the same name into the "results" directory.
# After theta is run on all config files, the analysis script scripts/plot_results.py is executed which
# requires the .db-files as input and writes the resulting plots as .pdf files to "results" as well.
# The numbers produced by the plot_results script are put in the text file "results/log.txt".
#
# This script will check whether the results files are already present and only run the examples
# again if they are not. To force re-generation or to clean up, simply delete all files in the "results" directory.

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
  execute_checked time ../../bin/theta $*
}

exec_theta_if_not_exists()
{
filename=$1
shift
[ -f $filename ] || exec_theta $*
}

# generate the input templates of the example model:
[ -f results/templates.root ] || ( execute_checked make -C scripts; execute_checked scripts/make_histos )

# Z_est distribution for signal + background for Figure 2.
# 100,000 pseudo experiments, runtime about 1 minute.
exec_theta_if_not_exists results/deltanll_hypo.db cfg/deltanll_hypo.cfg

# high statistics Z_est distribution for background only for Figure 3.
# 1,000,000 pseudo experiments, runtime about 10 minutes.
exec_theta_if_not_exists results/deltanll_hypo_bonly.db cfg/deltanll_hypo_bonly.cfg

# Bayes factors for signal + background for Figure 4.
# 10,000 pseudo experiments with Markov Chains with 10,000 iterations each; runtime about 20 minutes
exec_theta_if_not_exists results/bayesfactor.db cfg/bayesfactor.cfg

# Profile liklihood intervals for the text between Figure 4 and 5.
# 100,000 pseudo experiments, runtime about 18 minutes
exec_theta_if_not_exists results/profile_likelihood_intervals.db cfg/profile_likelihood_intervals.cfg

# Neyman construction using Z_est for Figure 5.
# 200,000 pseudo experiments, runtime about 2 minutes.
exec_theta_if_not_exists results/neyman_z_est.db cfg/neyman_z_est.cfg

# Neyman construction using the maximum likelihood estimate for mu_s for Figure 6:
# 200,000 pseudo experiments, runtime about 90 seconds.
exec_theta_if_not_exists results/neyman_mle.db cfg/neyman_mle.cfg

# Posterior in mu_s for Figure 7.
# 1 pseudo experiment (average data) with a Markov chain of 5,000,000 iterations; runtime about 30 seconds.
exec_theta_if_not_exists results/posterior.db cfg/posterior.cfg

# Produce maximum likelihood estimate for Neyman construction with systematic (rate only) uncertainties for Figure 9.
# 200,000 pseudo experiments, runtime about 90 seconds.
exec_theta_if_not_exists results/neyman_mle_syst.db cfg/neyman_mle_syst.cfg

# Produce maximum likelihood estimate for Neyman construction with systematic (template interpolation) uncertainties for Figure 9:
# 200,000 pseudo experiments, runtime about 6 minutes.
exec_theta_if_not_exists results/neyman_mle_syst_interp.db cfg/neyman_mle_syst_interp.cfg

# Produce all the pdf plots from the .db files generated previously.
# Runtime about 30 seconds.
./scripts/plot_results.py > ./results/log.txt

