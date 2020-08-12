#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

cp ../recon/onera-m6-sharp-te.egads .
cp ../bootstrap/onera-m6-sharp-te-adapt-surf.meshb .

${src}/ref adapt \
      onera-m6-sharp-te-adapt-surf.meshb \
      -g onera-m6-sharp-te.egads \
      --blend-metric 10000 \
      -x blend-adapt.meshb \
      -t -s 5

${src}/ref_blend_test \
      --viz blend-adapt.meshb onera-m6-sharp-te.egads

