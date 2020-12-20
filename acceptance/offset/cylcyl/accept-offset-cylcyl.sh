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
	  -g cylcyl.egads \
	  --facelift cylcyl-surrogate.meshb \
	  -m ${inproj}.solb \
	  -x ${outproj}.meshb \
	  -f ${outproj}-final.tec \
	  -s ${sweeps} \
	  ${tecplot}
    cp ref_gather_movie.tec ${inproj}_movie.tec
    ${src}/ref_acceptance ${field} ${outproj}.meshb ${outproj}.solb
    ${src}/ref_metric_test ${outproj}.meshb ${outproj}.solb > ${outproj}.status
}

# serveCSM -batch cylcyl.csm

${src}/ref boostrap cylcyl.egads \
      --facelift cylcyl-facelift.meshb

${src}/ref adapt cylcyl-adapt-surf.meshb \
      -g cylcyl.egads \
      --facelift-metric 100 \
      -s 5 \
      -x cylcyl-surrogate.meshb \
      -f cylcyl-surrogate-final.tec

${src}/ref_facelift_test --viz cylcyl-surrogate.meshb cylcyl.egads

${src}/ref adapt cylcyl-vol.meshb \
      -g cylcyl.egads \
      --facelift cylcyl-surrogate.meshb \
      -x cylcyl.meshb \
      -f cylcyl-final.tec

adapt_cycle cylcyl cycle01 2
adapt_cycle cycle01 cycle02 15

cat cycle02.status
../../check.rb cycle02.status 0.3 10.0

