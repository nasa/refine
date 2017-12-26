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

field=polar-2

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status
}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid

${two}/ref_acceptance -ugawg ${field} ega.meshb ega.metric

adapt_cycle ega cycle01 2
adapt_cycle cycle01 cycle02 5
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 15

cat cycle04.status
../../../check.rb cycle10.status 0.05 2.2

exit

adapt_cycle cycle11 cycle12 2
adapt_cycle cycle12 cycle13 2
adapt_cycle cycle13 cycle14 2
adapt_cycle cycle14 cycle15 2
adapt_cycle cycle15 cycle16 2
adapt_cycle cycle16 cycle17 2
adapt_cycle cycle17 cycle18 2
adapt_cycle cycle18 cycle19 2
adapt_cycle cycle19 cycle20 2
adapt_cycle cycle20 cycle21 2
adapt_cycle cycle21 cycle22 2
adapt_cycle cycle22 cycle23 2
adapt_cycle cycle23 cycle24 2
adapt_cycle cycle24 cycle25 2
adapt_cycle cycle25 cycle26 2
adapt_cycle cycle26 cycle27 2
adapt_cycle cycle27 cycle28 2
adapt_cycle cycle28 cycle29 2
adapt_cycle cycle29 cycle30 2

