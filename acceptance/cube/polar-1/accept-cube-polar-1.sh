#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/egads/one
    two=${HOME}/refine/egads/src
fi

field=polar-1

tecplot=-t

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_acceptance -ugawg ${field} ${inproj}.meshb ${inproj}.solb
    ${two}/ref_driver -i ${inproj}.meshb -g cube.egads -m ${inproj}.solb -o ${outproj} -s ${sweeps} ${tecplot}
    cp ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.solb
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.solb > ${outproj}.status
}

# ./cube.sh

adapt_cycle cube cycle01 2
adapt_cycle cycle01 cycle02 15
adapt_cycle cycle02 cycle03 15

cat cycle03.status
../../check.rb cycle03.status 0.3 3.0

