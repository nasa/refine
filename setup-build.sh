#!/usr/bin/env bash

set -x

mkdir -p strict
( cd strict && ../configure CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused' TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full' )