#!/usr/bin/env python

# Usage:
# copydeps.py <binary> <target path> <ld path1> <ld path 2> <ld path 3> ...
# copies the shared object the binary depends on to the target path, ignoring any
# libraries which are resolved to the given ld paths, if using these as the LD_LIBRARY_PATH.

import sys, os, os.path, subprocess, shutil, stat

if len(sys.argv) < 3: raise RuntimeError, "too few arguments"
binary = os.path.realpath(sys.argv[1])
target_path = os.path.realpath(sys.argv[2])
ld_paths = map(os.path.realpath, sys.argv[3:])

old_ld_path = os.environ.get('LD_LIBRARY_PATH', None)
os.environ['LD_LIBRARY_PATH'] = ':'.join(ld_paths)

if old_ld_path:
    os.environ['LD_LIBRARY_PATH'] += ':' + old_ld_path

out, err = '',''
p = subprocess.Popen(['ldd', binary], stderr = subprocess.PIPE, stdout=subprocess.PIPE)
while p.poll() is None:
   (tmpout, tmperr) = p.communicate()
   out += tmpout
   err += tmperr
code = p.wait()
if code!=0:
   print "error from 'ldd %s': " % binary, out, err
   sys.exit(1)

#print out
lines = out.split('\n')
for line in lines:
   if line.find('ld-linux') > -1:
      p = line.find('ld-linux')
      p_end = line.find(' ', p)
      p_start = p
      while p_start > 0 and line[p_start].strip()!="": p_start -= 1
      ld_path = line[p_start+1:p_end]
      shutil.copy2(ld_path, os.path.join(target_path, 'ld-linux.so'))
      continue
   if line.find('=>') == -1: continue
   so_name, so_path = line.split('=>', 2)
   so_path = so_path.strip()
   so_name = so_name.strip()
   if 'not found' in so_path: raise RuntimeError, "dependency %s not found." % so_name
   if so_path.find('(') != -1:
      so_path, dummy = so_path.split('(', 2)
      so_path = so_path.strip()
   if so_path=='': continue
   so_path = os.path.realpath(so_path)
   in_ld_paths = False
   for ld_path in ld_paths:
      if so_path.startswith(ld_path): in_ld_paths = True
   if in_ld_paths:
       #print "Skipping %s as it resolved already to one ld path" % so_path
       continue
   #print "Copying '%s' to %s" % (so_path, target_path)
   target = os.path.join(target_path, so_name)
   shutil.copy(so_path, target)
   sres = os.stat(target)
   os.chmod(target, sres.st_mode | stat.S_IRUSR | stat.S_IWUSR)

