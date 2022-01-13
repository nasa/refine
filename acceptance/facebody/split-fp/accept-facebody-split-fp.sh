\0;95;0c#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

tecplot=-t
field="-q fp-sa"
scalar="-u fp-sa"
egads="--egads split-fp.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3
    hrles="$4"

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}_volume.solb
    ${src}/ref_gather_test ${inproj}.meshb \
	  ${inproj}_volume.solb ${inproj}_volume.plt
    ${src}/ref_acceptance ${scalar} ${inproj}.meshb \
          ${inproj}_mach_100.solb
   ${src}/ref_acceptance ${scalar} ${inproj}.meshb \
          ${inproj}_mach_200.solb

    ${src}/ref loop ${inproj} ${outproj} ${complexity} \
	  ${egads} -s 5 \
	  --fixed-point _mach_ 100 100 200 \
	  ${hrles} | tee ${inproj}-loop.txt

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}_volume.solb
    ${src}/ref_gather_test ${outproj}.meshb \
	  ${outproj}_volume.solb ${outproj}_volume.plt
}

serveCSM -batch split-fp.csm > split-fp-servecsm.txt
${src}/ref bootstrap split-fp.egads > spit-fp-bootstrap.txt
mv split-fp-vol.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 1000 ""
adapt_cycle cycle01 cycle02 1000 ""
adapt_cycle cycle02 cycle03 2000 ""
adapt_cycle cycle03 cycle04 8000 \
	    "--hrles 0.2 5000000 --fun3d-mapbc split-fp-vol.mapbc"

${src}/ref distance cycle04.meshb cycle04-distance.solb \
      --fun3d-mapbc split-fp-vol.mapbc
${src}/ref_metric_test --hrles \
      cycle04.meshb cycle04-distance.solb cycle04_volume.solb

