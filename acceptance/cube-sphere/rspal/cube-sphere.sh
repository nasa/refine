#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

rm -f cube-sphere.egads
serveCSM -batch cube-sphere.csm

${src}/ref boostrap cube-sphere.egads
${src}/ref adapt cube-sphere-vol.meshb \
      -g cube-sphere.egads \
      -x cube-sphere00.meshb




