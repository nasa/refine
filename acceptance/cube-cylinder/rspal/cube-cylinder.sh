#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

rm -f cube-cylinder.egads
serveCSM -batch cube-cylinder.csm

${src}/ref boostrap cube-cylinder.egads
${src}/ref adapt cube-cylinder-vol.meshb \
      -g cube-cylinder.egads \
      -x cube-cylinder00.meshb




