#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/parmetis/src
fi

tecplot=-t
metric="-twod side"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3
    cores=$4

    ${src}/ref_acceptance ${metric} ${inproj}.meshb \
	  ${inproj}.solb

    mpiexec -np ${cores} \
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

serveCSM -batch square.csm
${src}/ref bootstrap square.egads
mv square-vol.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 10 2
adapt_cycle cycle01 cycle02 10 2

cat cycle02.status
../../check.rb cycle02.status 0.5 1.75

