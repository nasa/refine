#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

field="-ugawg side"

tecplot=-t

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${src}/ref_acceptance ${field} ${inproj}.meshb ${inproj}.solb
    ${src}/ref adapt ${inproj}.meshb \
	  -g sphcyl.egads \
	  --blend sphcyl-surrogate.meshb \
	  -m ${inproj}.solb \
	  -x ${outproj}.meshb \
	  -f ${outproj}-final.tec \
	  -s ${sweeps} \
	  ${tecplot}
    cp ref_gather_movie.tec ${inproj}_movie.tec
    ${src}/ref_acceptance ${field} ${outproj}.meshb ${outproj}.solb
    ${src}/ref_metric_test ${outproj}.meshb ${outproj}.solb > ${outproj}.status
}

serveCSM -batch sphcyl.csm

${src}/ref boostrap sphcyl.egads \
      --blend sphcyl-blend.meshb

${src}/ref adapt sphcyl-adapt-surf.meshb \
      -g sphcyl.egads \
      --blend-metric 100 \
      -s 5 \
      -x sphcyl-surrogate.meshb \
      -f sphcyl-surrogate-final.tec

${src}/ref adapt sphcyl-vol.meshb \
      -g sphcyl.egads \
      --blend sphcyl-surrogate.meshb \
      -x sphcyl.meshb \
      -f sphcyl-final.tec

adapt_cycle sphcyl cycle01 2
adapt_cycle cycle01 cycle02 15

cat cycle02.status
../../check.rb cycle02.status 0.3 3.0

