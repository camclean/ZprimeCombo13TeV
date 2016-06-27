#!/bin/bash

curr_dir=$PWD
tmpdir=$(mktemp -d)
logfile=$tmpdir/log.txt

fail()
{
 echo "FAIL: " $* | tee -a $logfile;
 cp $logfile $curr_dir;
 echo "note: logfile was copied to ${curr_dir}/log.txt"
 exit 1;
}

if [ $? -gt 0 ] || [ ! -d "$tmpdir" ]; then
   fail creating tempdir
fi
echo [0] Using directory $tmpdir
cd $tmpdir
echo [1] Subversion checkout
svn --non-interactive co https://ekptrac.physik.uni-karlsruhe.de/public/theta/tags/testing theta &> $logfile
[ $? -eq 0 ] || fail subversion checkout
cd theta
echo [2] Compiling
make &> $logfile || fail Compiling
echo [3] Running tests
test/testall.sh &> $logfile || fail tests
echo SUCCESS

#TODO: copy back logfile
#valgrind is not run because ROOT spoils it all ...
#valgrind --leak-check=full --show-reachable=yes test/test >> $logfile

#lcov -c --directory build-coverage --output-file theta.info
#genhtml theta.info --no-function-coverage -o doc/coverage

