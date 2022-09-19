#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

field="-u sin50xy"

${src}/ref_fixture_test --prism prism.lb8.ugrid \
		 0 1 0 1 0 1 \
		 11 11 11 0

${src}/ref_shard_test prism.lb8.ugrid 1 5

${src}/ref translate ref_shard_test.b8.ugrid mixed01.lb8.ugrid

${src}/ref_acceptance ${field} mixed01.lb8.ugrid \
      mixed01.solb

${src}/ref multiscale mixed01.lb8.ugrid mixed01.solb \
	  500 mixed01-metric.solb

${src}/ref adapt mixed01.lb8.ugrid \
      --metric mixed01-metric.solb \
      -x mixed02.lb8.ugrid \
      -x mixed02.plt -s 10

