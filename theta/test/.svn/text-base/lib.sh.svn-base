
fail()
{
   echo "FAIL: $*";
   exit 1;
}

warn()
{
   echo "WARNING: $*"
}

pass()
{
   echo "PASS: $*"
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
    execute_checked ../../bin/theta -q $*
}
