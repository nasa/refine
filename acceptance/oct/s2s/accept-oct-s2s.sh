#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

project=s2s

serveCSM -batch -skipTess ${project}.csm
${src}/ref bootstrap ${project}.egads -s 5

${src}/ref_oct_test --surf s2s-adapt-surf.meshb oct-s2s.tec

