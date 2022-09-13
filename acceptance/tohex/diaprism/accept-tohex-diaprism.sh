#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

project=diaprism

serveCSM -batch -skipTess ${project}.csm
${src}/ref bootstrap ${project}.egads -s 5

${src}/ref translate \
      ${project}-vol.meshb \
      ${project}-hex.ugrid \
      --blockhead 

${src}/ref translate \
      ${project}-hex.ugrid \
      ${project}-hex.plt

