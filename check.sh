#!/usr/bin/env bash

set -x

script_name_with_path=`which "${0}"`
script_path=`dirname "${script_name_with_path}"`

# TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full'

if test -d ${script_path}/mpi
then
  ( cd ${script_path}/mpi/two && \
     make -j check ) || exit
fi

if test -d ${script_path}/strict
then
  ( cd ${script_path}/strict/two && \
    make -j check TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full' ) || exit

fi

