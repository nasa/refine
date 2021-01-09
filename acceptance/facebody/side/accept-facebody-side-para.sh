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
metric="-twod side"
egads="-g square.egads"

function adapt_para {
    inproj=$1
    outproj=$2
    sweeps=$3
    cores=$4

    ${src}/ref_acceptance ${metric} ${inproj}.meshb \
	  ${inproj}.solb

    mpiexec -np ${cores} \
	    ${src}/refmpi adapt \
	    ${inproj}.meshb \
	    -m ${inproj}.solb \
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
cp square.egads para.egads
${src}/ref bootstrap para.egads
mv para-vol.meshb para00.meshb

adapt_para para00 para01 10 2
adapt_para para01 para02 10 2

cat para02.status
../../check.rb para02.status 0.5 1.75

