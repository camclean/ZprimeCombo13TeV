#!/bin/bash

# keep same directory structure: lib for plugins or theta-built shared objects and bin for
# binaries.
# add extlib for external libraries

[ -e ../bin/theta ] || { echo "../bin/theta not found. Compile theta and execute this script from within its directory!"; exit 1; }


echo "[0] creating gridpack-tmp"
rm -rf gridpack-tmp
mkdir -p gridpack-tmp/lib
mkdir -p gridpack-tmp/bin
mkdir -p gridpack-tmp/extlib
mkdir -p gridpack-tmp/etc/plugins

if [ ! -z "$ROOTSYS" ]; then
   echo "[0.5] Copying essential ROOT dependencies"
   cp -r $ROOTSYS/etc/plugins/TVirtualStreamerInfo gridpack-tmp/etc/plugins/
   cp $ROOTSYS/etc/system.rootrc gridpack-tmp/etc/
fi

# add the dependencies:
echo "[1] copying theta to gridpack-tmp"
cp scripts/theta gridpack-tmp/bin/
cp ../bin/theta gridpack-tmp/bin/theta.exe
cp ../lib/*.so gridpack-tmp/lib/

echo "[2] copying dependencies to gridpack-tmp"
scripts/copydeps.py gridpack-tmp/bin/theta.exe gridpack-tmp/extlib gridpack-tmp/lib || { echo "error executing copydeps"; exit 1; }
for plugin in gridpack-tmp/lib/*.so; do
   scripts/copydeps.py $plugin gridpack-tmp/extlib gridpack-tmp/lib gridpack-tmp/extlib || { echo "error executing copydeps"; exit 1; }
done
for so_file in gridpack-tmp/extlib/*.so; do
   if [ $so_file == "gridpack-tmp/extlib/ld-linux.so" ]; then continue; fi
   scripts/copydeps.py $so_file gridpack-tmp/extlib gridpack-tmp/lib gridpack-tmp/extlib  || { echo "error executing copydeps"; exit 1; }
done

if [ ! -f gridpack-tmp/extlib/ld-linux.so ]; then
   echo "Did not find ld-linux!";
   exit 1;
fi

echo "[3] stripping all symbols to reduce size"
strip -s gridpack-tmp/bin/theta.exe gridpack-tmp/lib/* gridpack-tmp/extlib/*

echo "[4] creating gridpack.tgz from gridpack-tmp"
rm -f gridpack.tgz
cd gridpack-tmp
tar zcf ../gridpack.tgz *

FILEOUT=`file ./bin/theta.exe`
VER=`echo ${FILEOUT/, for GNU\//%} | cut -d% -f2 | cut -d, -f1`
echo "NOTE: minimal version this will run on:" $VER

