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
field="-u parabola"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref multiscale ${inproj}.meshb ${inproj}.solb \
	  1000 ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}.solb
}

serveCSM -batch square.csm
${src}/ref bootstrap square.egads
mv square-vol.meshb func00.meshb

adapt_cycle func00 func01
adapt_cycle func01 func02
adapt_cycle func02 func03
adapt_cycle func03 func04


