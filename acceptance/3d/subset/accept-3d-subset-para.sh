#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/parmetis/src
fi

${src}/ref_acceptance 1 whole.meshb
${src}/ref_acceptance -u u5 whole.meshb whole.solb

${src}/ref translate whole.meshb target.meshb \
      --scale 0.5 --shift 0.25 0.25 0.25
${src}/ref_acceptance -u u5 target.meshb target.solb

mpiexec -np 4 ${src}/ref_gather_test --subset whole.meshb whole.solb \
	0.24 0.2 0.2 0.02 0.6 0.6 \
	subset.meshb subset.solb
${src}/ref_gather_test whole.meshb whole.solb whole.tec
${src}/ref_gather_test subset.meshb subset.solb subset.tec

${src}/ref translate subset.meshb subset.lb8.ugrid
${src}/ref translate subset.meshb subset-trans.tec

mpiexec -np 4 ${src}/ref_interp_test --face subset.meshb subset.solb \
      target.meshb target.solb target_merged.solb 1 

${src}/ref_gather_test target.meshb \
      target_merged.solb target_merged.tec

exit
