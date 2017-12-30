#!/usr/bin/env bash

set -x
set -u

jobs=${1:-}

script_name_with_path=`which "${0}"`
script_path=`dirname "${script_name_with_path}"`

if test -d ${script_path}/egads
then
  ( cd ${script_path}/egads && \
     make -j ${jobs} check ) || exit
  ( cd ${script_path}/egads && \
     make -j ${jobs} install ) || exit
fi

if test -d ${script_path}/zoltan
then
  ( cd ${script_path}/zoltan && \
     make -j ${jobs} check ) || exit
  ( cd ${script_path}/zoltan && \
     make -j ${jobs} install ) || exit
fi

if test -d ${script_path}/parmetis
then
  ( cd ${script_path}/parmetis && \
     make -j ${jobs} check ) || exit
  ( cd ${script_path}/parmetis && \
     make -j ${jobs} install ) || exit
fi

if test -d ${script_path}/strict
then
  ( cd ${script_path}/strict && \
    make -j ${jobs} check TESTS_ENVIRONMENT='valgrind --quiet  --error-exitcode=1 --leak-check=full' ) || exit
  ( cd ${script_path}/strict && \
    make -j ${jobs} install ) || exit
fi

