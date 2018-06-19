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

${two}/ref_driver -i skinny-cylinder.meshb -g ${geomfile} -o ref_driver1 -r 1 -t
mv ref_gather_movie.tec ref_driver1_movie.tec

${two}/ref_driver -i ref_driver1.meshb -g ${geomfile} -o ref_driver2 -r 10 -t
mv ref_gather_movie.tec ref_driver2_movie.tec

${two}/ref_histogram_test ref_driver2.meshb ref_driver2-final-metric.solb \
 > accept-skinny-cylinder-curvature-02.status

cat accept-skinny-cylinder-curvature-02.status
../../check.rb accept-skinny-cylinder-curvature-02.status 0.1 2.2



