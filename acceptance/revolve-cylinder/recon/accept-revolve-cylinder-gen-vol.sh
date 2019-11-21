#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

geomfile=cylinder.egads

${src}/ref boostrap ${geomfile} cylinder-vol.meshb \
    | tee accept-revolve-cylinder-gen-vol.status

