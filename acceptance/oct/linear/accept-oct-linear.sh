#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

function adapt_cycle {
    input=$1
    output=$2
    hmin=0.01

    ${src}/ref_acceptance ${input}.meshb ${input}-metric.solb ${hmin}

    ${src}/ref_oct_test --adapt ${input}.meshb ${input}-metric.solb \
	  ${output}.meshb
    
    ${src}/ref_acceptance ${output}.meshb ${output}-metric.solb ${hmin}
}

${src}/ref_oct_test --box 0.05 cycle00.meshb

adapt_cycle cycle00 cycle01
adapt_cycle cycle01 cycle02

