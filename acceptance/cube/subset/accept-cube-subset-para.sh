#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/parmetis/src
fi

cp ../initial-cell/cube01.meshb cube.meshb

${src}/ref_acceptance -u u5 cube.meshb cube.solb

${src}/ref_acceptance 1 whole.meshb
${src}/ref translate whole.meshb target.meshb \
      --scale 0.5 --shift 0.25 0.25 0.25
${src}/ref_acceptance -u u5 target.meshb target.solb

mpiexec -np 4 ${src}/ref_gather_test --subset cube.meshb cube.solb \
	0.24 0.2 0.2 0.02 0.6 0.6 \
	subset.meshb subset.solb
${src}/ref_gather_test cube.meshb cube.solb whole.tec
${src}/ref_gather_test subset.meshb subset.solb subset.tec

${src}/ref translate subset.meshb subset.lb8.ugrid
${src}/ref translate subset.meshb subset-trans.tec

mpiexec -np 4 ${src}/ref interpolate subset.meshb subset.solb \
      target.meshb target_merged.solb --face 1 target.solb

${src}/ref_gather_test target.meshb \
      target_merged.solb target_merged.tec

exit
