#!/usr/bin/env bash

CC=gcc
CFLAGS="-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused"

if [ -z "$1" ]
then
  tests=`ls -1 *_test.c`
else
  tests=$*
fi

for test in $tests;
do
  root=${test%_test.c}
  echo \$ $0 $root

  dependencies=''
  if [ -a ${root}.h ]; then
    for dep_header in `grep '#include "ref_' ${root}.h`
    do
      source_with_leading_quote=${dep_header%.h\"}.c
      source=${source_with_leading_quote#\"}
      if [ -a ${source} ]; then
        if ! echo ${dependencies} | grep -q "${source}" ; then
          dependencies=${dependencies}' '${source}
        fi
      fi
    done
  fi
  for dep_header in `grep '#include "ref_' ${root}_test.c`
  do
    source_with_leading_quote=${dep_header%.h\"}.c
    source=${source_with_leading_quote#\"}
    if [ -a ${source} ]; then
      if ! echo ${dependencies} | grep -q "${source}" ; then
        dependencies=${dependencies}' '${source}
      fi
    fi
  done

  compile="${CC} ${CFLAGS} -o ${root}_test ${dependencies} ${root}_test.c -lm "
  (eval ${compile} && eval valgrind --quiet --leak-check=full ./${root}_test) || \
      ( echo ${compile} && echo FAIL: to re-run, \$ $0 ${root} ; exit 1 ) ||  exit 1
done

