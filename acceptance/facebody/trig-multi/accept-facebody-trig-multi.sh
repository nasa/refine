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
field="-q trig"
egads="--egads square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}_volume.solb
    ${src}/ref_gather_test ${inproj}.meshb \
	  ${inproj}_volume.solb ${inproj}_volume.tec

    ${src}/ref loop ${inproj} ${outproj} 1000 \
	  ${egads} -s 5 > ${inproj}-loop.txt

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}_volume.solb
    ${src}/ref_gather_test ${outproj}.meshb \
	  ${outproj}_volume.solb ${outproj}_volume.tec
}

serveCSM -batch square.csm
${src}/ref bootstrap square.egads
mv square-vol.meshb cycle00.meshb

adapt_cycle cycle00 cycle01
adapt_cycle cycle01 cycle02
adapt_cycle cycle02 cycle03
adapt_cycle cycle03 cycle04
adapt_cycle cycle04 cycle05
adapt_cycle cycle05 cycle06

../../check.rb cycle05-loop.txt 0.5 2.0

