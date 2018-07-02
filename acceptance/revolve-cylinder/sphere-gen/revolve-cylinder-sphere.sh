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

serveCSM -batch revolve-cylinder-sphere.csm

${two}/ref_geom_test --tess \
    revolve-cylinder-sphere.egads \
    revolve-cylinder-sphere.meshb \
    5 5 90

