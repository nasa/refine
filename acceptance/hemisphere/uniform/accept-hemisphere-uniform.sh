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

geomfile=hemisphere.egads

# ~/esp/EngSketchPad/bin/serveCSM -batch hemisphere.csm
# ref_geom_test hemisphere.egads hemisphere.meshb
# ref adapt hemisphere.meshb -g hemisphere.egads -r 4 -x hemicurve.meshb

${two}/ref_acceptance hemicurve.meshb hemicurve-metric.solb 0.1
${two}/ref adapt hemicurve.meshb -g ${geomfile} -m hemicurve-metric.solb -x hemicurve1.meshb
${two}/ref_acceptance hemicurve1.meshb hemicurve1-metric.solb 0.1
${two}/ref_metric_test hemicurve1.meshb hemicurve1-metric.solb > accept-hemisphere-uniform-01.status

cat accept-hemisphere-uniform-01.status
../../check.rb accept-hemisphere-uniform-01.status 0.3 3.0

${two}/ref_interp_test hemicurve.meshb hemicurve1.meshb
${two}/ref_interp_test hemicurve1.meshb hemicurve.meshb

