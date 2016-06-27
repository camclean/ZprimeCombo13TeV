#!/usr/bin/env python

# usage: convert-databse.py <input file 1> ... <input file n>  <output file>
#
# converts and merges input files 1 through n to the output database file specified as
# last command line argument. The file type is guessed from the extension: it is assumed to
# be a ROOT file for the "root" extension and a sqlite file for "db" extension.
# As input file type, currently only sqlite is supported. As output file type, root files and
# sqlite files are supported.

from theta_auto.theta_interface import settingvalue_to_cfg
import sys, os, os.path, tempfile

if len(sys.argv) < 3:
    print "Usage: %s <infile 1> .. <infile n>  <outfile>" % sys.argv[0]
    print "Specify at least one input file; the last argument will be used as output file."
    sys.exit(1)

infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

for ifile in infiles:
    if not ifile.endswith('.db'):
        raise RuntimeError, "only .db files supported as input"

outtype = ''
if outfile.endswith('.db'): outtype = 'sqlite_database'
elif outfile.endswith('.root'): outtype = 'rootfile_database'

if outtype == '': raise RuntimeError, "unknown output file extension"

main = dict(type = 'convert_database', input_database = dict(type = 'sqlite_database_in', filenames = infiles),
   output_database = dict(type = outtype, filename = outfile))

options = {'plugin_files': ['$THETA_DIR/lib/core-plugins.so']}
if outtype == 'rootfile_database': options['plugin_files'].append('$THETA_DIR/lib/root.so')

f = tempfile.NamedTemporaryFile(delete = False)
print 'Creating theta config file ', f.name

print >> f, 'main = ', settingvalue_to_cfg(main), ';'
print >> f, 'options = ', settingvalue_to_cfg(options), ';'
f.flush()

theta_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
retval = os.system(os.path.join(theta_dir, 'bin', 'theta') + " " + f.name)

if retval != 0: raise RuntimeError, "executing theta failed with exit code %d" % retval
else: print "done"

os.unlink(f.name)


