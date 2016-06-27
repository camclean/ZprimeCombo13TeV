#!/usr/bin/env bash

#run this from the theta main directory!

[ -x bin/theta ] || { echo "theta main executable in 'bin/theta' not found (are we in the theta root directory?)"; exit 1; }

THETA_DIR=$PWD

. test/lib.sh

[ -x root/create_testhistos ] && root/create_testhistos
execute_checked bin/test
rm -f testhistos.root

[ "$*" = "unit" ] && exit 0;

fail=0

for i in test/test-stat/*.py; do
   i=`basename $i`
   echo "START `date +%s`;     `date -R`     $i"
   export output=$( cd test/test-stat; ./$i 2>&1 )
   echo "output: ${output}"
   output_nonpass=$( echo "$output" | grep -v "PASS:" )
   if [ -n "$output_nonpass" ]; then
       echo -e "\n\e[0;31m $i FAILED: ${output_nonpass}\e[m \n";
       fail=$(($fail+1))
   fi
   echo "$output" | awk "{print \"[$i]\", \$_}"
   echo "END `date +%s`;     `date -R`     $i"
done

echo "executing theta-auto tests"
cd $THETA_DIR/utils/test || { echo "did not find utils test dir!"; exit 1; }
rm -rf ./test/cache
a=$( ../theta-auto.py test.py | grep Failures )
[ $? -gt 0 ] && { echo "error executing theta-auto tests!"; exit 1; }
eval $a
fail=$(($fail+$Failures))

echo "executing theta-auto2 tests"
cd $THETA_DIR/utils2/test || { echo "did not find utils2 test dir!"; exit 1; }
rm -rf ./test/cache
a=$( ../theta-auto.py test.py | grep Failures )
[ $? -gt 0 ] && { echo "error executing theta-auto2 tests!"; exit 1; }
eval $a
fail=$(($fail+$Failures))



echo "Failures: ${fail}"
[ $fail -eq 0 ] || exit 1
