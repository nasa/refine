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
if [ ! -e $arch/acinclude.m4 ] ; then
  echo "***********************************************************************"
  echo "**"
  echo "** WARNING: No \"acinclude.m4\" file found in \"$arch/\"."
  echo "**          You may wish to copy \"libtool.m4\" into"
  echo "**          \"$arch/acinclude.m4\" from your libtool/aclocal"
  echo "**          installation for the most current libtool configuration."
  echo "**"
  echo "**          Check your libtool installation for location."
  echo "**"
  echo "**          Example: \"/usr/local/libtool/share/aclocal/libtool.m4\""
  echo "**                   or \"/usr/share/aclocal/libtool.m4\""
  echo "**"
  echo "***********************************************************************"
fi

./bootstrap

echo "Running ../configure --prefix=`pwd` $* ..."
( cd $arch && ../configure --prefix=`pwd` $* )
