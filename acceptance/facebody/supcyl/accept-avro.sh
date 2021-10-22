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
field="supcyl"
project="supcyl"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3

    ${src}/ref_acceptance -u ${field} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref multiscale ${inproj}.meshb ${inproj}.solb \
	  ${complexity} ${inproj}-metric.solb

    avro -adapt ${inproj}.meshb ${project}.egads ${inproj}-metric.solb \
	 ${outproj}.meshb
    transmesh ${outproj}_0.mesh ${outproj}.meshb

#    ${src}/ref adapt ${inproj}.meshb \
#	  -g ${project}.egads \
#	  -m ${inproj}-metric.solb \
#	  -x ${outproj}.meshb

    ${src}/ref_acceptance -u ${field} ${outproj}.meshb \
	  ${outproj}.solb
    ${src}/ref multiscale ${outproj}.meshb ${outproj}.solb \
          ${complexity} ${outproj}-metric.solb

    ${src}/ref visualize ${outproj}.meshb ${outproj}.solb ${outproj}.plt
}

serveCSM -batch ${project}.csm
${src}/ref bootstrap ${project}.egads
cp ${project}-vol.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 1000
adapt_cycle cycle01 cycle02 1000
adapt_cycle cycle02 cycle03 1000
