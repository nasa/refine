#!/usr/bin/env bash

set -x

script_name_with_path=`which "${0}"`
script_path=`dirname "${script_name_with_path}"`

# TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full'

( cd ${script_path}/mpi/two && \
    make check ) && \
( cd ${script_path}/strict/two && \
    make check TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full' ) && \
true

