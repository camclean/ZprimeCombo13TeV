# -*- coding: utf-8 -*-

import sys, os, sqlite3, math

def fail(s):
    print "FAIL: %s" % s
    sys.exit(1)

def warn(s):
    print "WARNING: %s" % s
    sys.exit(1)

def passed(s):
    print "PASS: %s" % s


def execute_checked(command):
    i=os.system(command)
    if i>0: fail("executing %s returned error" % command)

def exec_theta(s, quiet = True):
    prefix = ""
    if os.getenv("THETA_MEMCHECK") is not None: prefix="valgrind --leak-check=full"
    if os.getenv("THETA_PROFILE") is not None: prefix="valgrind --tool=callgrind"
    execute_checked("%s ../../bin/theta %s %s" % (prefix, s, '-q' if quiet else ''))


def sql(filename, query):
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    c.execute(query)
    result = c.fetchall()
    c.close()
    conn.close()
    return result

