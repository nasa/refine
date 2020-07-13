#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

tecplot=-t
metric="-twod polar-2"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${src}/ref_acceptance ${metric} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref_driver \
	  -i ${inproj}.meshb \
	  -m ${inproj}.solb \
	  ${egads} \
          -x ${outproj}.meshb \
	  -s ${sweeps} ${tecplot}

    mv ref_gather_histo.tec ${outproj}_histo.tec
    mv ref_gather_movie.tec ${outproj}_movie.tec
    ${src}/ref_acceptance ${metric} ${outproj}.meshb \
	  ${outproj}.solb
    ${src}/ref_metric_test ${outproj}.meshb ${outproj}.solb \
	  > ${outproj}.status
}

cp square.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 10
adapt_cycle cycle01 cycle02 10
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 10

cat cycle04.status
../../check.rb cycle04.status 0.30 2.5

