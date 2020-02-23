#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/parmetis/one
    two=${HOME}/refine/parmetis/src
fi

tecplot=-t
metric="-twod linear-0001"

function adapt_para {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_acceptance ${metric} ${inproj}.meshb \
	  ${inproj}.solb

    mpiexec -np 2 ${two}/ref_driver -i ${inproj}.meshb -m ${inproj}.solb \
          -x ${outproj}.meshb \
	  -s ${sweeps} ${tecplot}

    mv ref_gather_histo.tec ${outproj}_histo.tec
    mv ref_gather_movie.tec ${outproj}_movie.tec
    ${two}/ref_acceptance ${metric} ${outproj}.meshb \
	  ${outproj}.solb
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.solb \
	  > ${outproj}.status
}

inproj=para00
${two}/ref_acceptance 2 para00.meshb

adapt_para para00 para01 10
adapt_para para01 para02 10
adapt_para para02 para03 10
adapt_para para03 para04 10

cat para04.status
../../check.rb para04.status 0.40 1.6

