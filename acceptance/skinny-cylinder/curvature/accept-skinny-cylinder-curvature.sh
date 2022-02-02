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

geomfile=skinny-cylinder.egads

# ${two}/ref_geom_test ${geomfile} skinny-cylinder.meshb 1 1 90

${two}/ref adapt skinny-cylinder.meshb -g ${geomfile} -x skinny-cylinder1.meshb -t
mv ref_gather_movie.tec skinny-cylinder1_movie.tec

${two}/ref adapt skinny-cylinder1.meshb -g ${geomfile} \
      -x skinny-cylinder2.meshb \
      --export-metric-as skinny-cylinder2-final-metric.solb \
      -t
mv ref_gather_movie.tec skinny-cylinder2_movie.tec

${two}/ref_histogram_test skinny-cylinder2.meshb skinny-cylinder2-final-metric.solb \
 > accept-skinny-cylinder-curvature-02.status

cat accept-skinny-cylinder-curvature-02.status
../../check.rb accept-skinny-cylinder-curvature-02.status 0.1 2.2



