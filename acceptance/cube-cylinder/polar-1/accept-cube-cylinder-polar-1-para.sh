#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/zoltan/one
    two=${HOME}/refine/zoltan/src
fi

field=polar-1

function adapt_cycle_sant {
    inproj=$1
    outproj=$2
    sweeps=$3

    mpiexec -np 4 ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -l -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}
function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    mpiexec -np 4 ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}

# ${two}/ref_geom_test ega.egads ega.meshb

${two}/ref_acceptance -ugawg ${field} ega.meshb ega.metric

adapt_cycle_sant ega cycle01 2
adapt_cycle_sant cycle01 cycle02 2
adapt_cycle_sant cycle02 cycle03 2
adapt_cycle_sant cycle03 cycle04 2
adapt_cycle_sant cycle04 cycle05 2
adapt_cycle_sant cycle05 cycle06 2
adapt_cycle cycle06 cycle07 2
adapt_cycle cycle07 cycle08 2
adapt_cycle cycle08 cycle09 2
adapt_cycle cycle09 cycle10 2
adapt_cycle cycle10 cycle11 2
adapt_cycle cycle11 cycle12 2
adapt_cycle cycle12 cycle13 2
adapt_cycle cycle13 cycle14 2
adapt_cycle cycle14 cycle15 2
adapt_cycle cycle15 cycle16 2
adapt_cycle cycle16 cycle17 2
adapt_cycle cycle17 cycle18 2
adapt_cycle cycle18 cycle19 2
adapt_cycle cycle19 cycle20 2
adapt_cycle cycle20 cycle21 10
adapt_cycle cycle21 cycle22 10
adapt_cycle cycle22 cycle23 10
adapt_cycle cycle23 cycle24 10
adapt_cycle cycle24 cycle25 10
adapt_cycle cycle25 cycle26 10
adapt_cycle cycle26 cycle27 10
adapt_cycle cycle27 cycle28 10
adapt_cycle cycle28 cycle29 10
adapt_cycle cycle29 cycle30 10

cat cycle30.status
../../../check.rb cycle30.status 0.060 10.0

