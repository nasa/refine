#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

project=diamond

serveCSM -batch -skipTess ${project}.csm
${src}/ref bootstrap ${project}.egads -s 5

${src}/ref translate \
      ${project}-vol.meshb \
      ${project}-extrude-hex.ugrid \
      --blockhead 

${src}/ref translate \
      ${project}-extrude-hex.ugrid \
      ${project}-extrude-hex.plt

