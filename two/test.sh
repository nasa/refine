#!/usr/bin/env bash

CC=gcc
CFLAGS="-g -O0 -pedantic-errors -Wall -Wextra -Werror -Wunused"

tests=`ls -1 *_test.c`

for test in $tests;
do
  root=${test%_test.c}
  echo $root
  compile="${CC} ${CFLAGS} -o ${root}_test ${root}.c ${root}_test.c"
  eval ${compile} && eval ./${root}_test
done