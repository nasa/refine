#!/bin/csh
# Run this to generate all the initial makefiles, etc.

set arch = `uname`

echo "Configure $arch..."

echo "Linking required AutoTool files..."
( cd $arch && ln -s ../AutoConfInput configure.ac )
( cd $arch && ln -s ../AUTHORS . )
( cd $arch && ln -s ../ChangeLog . )
( cd $arch && ln -s ../NEWS . )
( cd $arch && ln -s ../COPYING . )
( cd $arch && ln -s ../INSTALL . )
( cd $arch && ln -s ../README . )

echo "Running aclocal ..."
( cd $arch && aclocal )

echo "Running automake ..."
( cd $arch && automake --add-missing )

echo "Running autoconf ..."
( cd $arch && autoconf )

( cd $arch && echo "Running ./configure --prefix=`pwd` " $* )
( cd $arch && ./configure --prefix=`pwd` $* )
