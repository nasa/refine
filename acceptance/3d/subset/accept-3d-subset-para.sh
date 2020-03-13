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

mpiexec -np 4 ${src}/ref_gather_test --subset whole.meshb whole.solb -1 -1 -1 1.4 1.6 1.8 subset.meshb subset.solb
${src}/ref_gather_test whole.meshb whole.solb whole.tec
${src}/ref_gather_test subset.meshb subset.solb subset.tec

exit
