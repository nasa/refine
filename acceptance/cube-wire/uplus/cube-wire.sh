#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

rm -f cube-wire.egads
serveCSM -batch cube-wire.csm

${src}/ref boostrap cube-wire.egads
${src}/ref adapt cube-wire-vol.meshb \
      -g cube-wire.egads \
      -x cube-wire00.meshb




