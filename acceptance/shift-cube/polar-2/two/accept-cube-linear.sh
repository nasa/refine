#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/src
    two=$1/two
else
    one=${HOME}/refine/egads/src
    two=${HOME}/refine/egads/two
fi

field=polar-2

tecplot=-t

function adapt_cycle_sant {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.b8.ugrid -g cube.egads -p ${inproj}.gas -m ${inproj}.metric -o ${outproj} -s ${sweeps} -l ${tecplot}
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.b8.ugrid ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric > ${outproj}.status

}
function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.b8.ugrid -g cube.egads -p ${inproj}.gas -m ${inproj}.metric -o ${outproj} -s ${sweeps} ${tecplot}
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.b8.ugrid ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric > ${outproj}.status

}

# ./cube.sh

${two}/ref_translate cube.ugrid cube.b8.ugrid
${two}/ref_acceptance -ugawg ${field} cube.b8.ugrid cube.metric

adapt_cycle_sant cube cycle01 2
adapt_cycle_sant cycle01 cycle02 2
adapt_cycle_sant cycle02 cycle03 2
adapt_cycle_sant cycle03 cycle04 2
adapt_cycle_sant cycle04 cycle05 2
adapt_cycle_sant cycle05 cycle06 2
adapt_cycle cycle06 cycle07 10
adapt_cycle cycle07 cycle08 10
adapt_cycle cycle08 cycle09 10
adapt_cycle cycle09 cycle10 10
adapt_cycle cycle10 cycle11 10
adapt_cycle cycle11 cycle12 10
adapt_cycle cycle12 cycle13 10
adapt_cycle cycle13 cycle14 10
adapt_cycle cycle14 cycle15 10
adapt_cycle cycle15 cycle16 10
adapt_cycle cycle16 cycle17 10
adapt_cycle cycle17 cycle18 10
adapt_cycle cycle18 cycle19 10
adapt_cycle cycle19 cycle20 10
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
../../../check.rb cycle30.status 0.08 1.8

