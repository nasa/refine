#!/bin/sh
#Run to generate bootstrp AutoTools and run configure in a archetecture dir

arch=`uname`

echo "Setting up $arch..."

if [ ! -d $arch ]; then
  echo
  echo "Creating \"$arch\" directory..."
  mkdir $arch
  echo
fi

echo "Running bootstrap ..."
./bootstrap

echo "Running ../configure --prefix=`pwd` $* ..."
( cd $arch && ../configure --prefix=`pwd` $* )
