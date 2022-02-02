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

geomfile=cone-cone.egads

# ${two}/ref_geom_test ${geomfile} cone-cone.meshb

${two}/ref_acceptance cone-cone.meshb cone-cone.metric 0.1
${two}/ref adapt cone-cone.meshb -g ${geomfile} -m cone-cone.metric \
      -x cone-cone1.meshb \
      --export-metric-as cone-cone1-final-metric.solb \
      -t
mv ref_gather_movie.tec cone-cone1_movie.tec
${two}/ref_acceptance cone-cone1.meshb cone-cone1.metric 0.1
${two}/ref_metric_test cone-cone1.meshb cone-cone1-final-metric.solb > accept-cone-cone-uniform-01.status

${two}/ref adapt cone-cone1.meshb -g ${geomfile} -m cone-cone1.metric \
      -x cone-cone2.meshb \
      --export-metric-as cone-cone2-final-metric.solb \
      -t
mv ref_gather_movie.tec cone-cone2_movie.tec
${two}/ref_acceptance cone-cone2.meshb cone-cone2.metric 0.1
${two}/ref_metric_test cone-cone2.meshb cone-cone2-final-metric.solb > accept-cone-cone-uniform-02.status

cat accept-cone-cone-uniform-02.status
../../check.rb accept-cone-cone-uniform-02.status 0.20 6.0

