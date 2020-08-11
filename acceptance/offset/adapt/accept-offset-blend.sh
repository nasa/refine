#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

cp ../side/offset.egads .
cp ../side/offset-adapt-surf.meshb .

${src}/ref_blend_test \
      --metric offset-adapt-surf.meshb offset.egads 100

${src}/ref adapt \
      offset-adapt-surf.meshb \
      -g offset.egads \
      --blend-metric 100 \
      -x blend-adapt.meshb \
      -t -s 10

${src}/ref_blend_test \
      --viz blend-adapt.meshb offset.egads
