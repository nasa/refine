#!/bin/sh
# Run this to generate all the initial makefiles, etc.

echo "linking configure.ac to configure.in..."
ln -s configure.ac configure.in

echo "Running aclocal..."
aclocal

echo "Running automake..."
automake --add-missing

echo "Running autoconf..."
autoconf

echo "Running ./configure..."
./configure --prefix=`pwd`
