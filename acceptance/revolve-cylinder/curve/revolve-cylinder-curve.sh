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

${two}/ref adapt \
    revolve-cylinder-sphere.meshb \
    -g ../sphere-gen/revolve-cylinder-sphere.egads \
    -x revolve-cylinder-sphere-curve.meshb \
    -t

