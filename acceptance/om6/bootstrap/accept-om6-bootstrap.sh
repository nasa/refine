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

cp ../recon/onera-m6-sharp-te.egads .

${two}/ref bootstrap onera-m6-sharp-te.egads --blend om6-blend.meshb

${two}/ref adapt \
    onera-m6-sharp-te-vol.meshb \
    -g onera-m6-sharp-te.egads \
    -x onera-m6-sharp-te01.meshb \
    --export-metric-as onera-m6-sharp-te01-final-metric.solb

${two}/ref_metric_test \
    onera-m6-sharp-te01.meshb \
    onera-m6-sharp-te01-final-metric.solb \
    > vol-adapt.status

../../check.rb vol-adapt.status 0.075 3.3


