#!/bin/sh
# Run this to generate all the initial makefiles, etc.

echo "Running aclocal ..."
aclocal

echo "Running automake ..."
automake --add-missing

echo "Running autoconf ..."
autoconf

echo "Running ./configure --prefix=`pwd` ..."
./configure --prefix=`pwd` $1
