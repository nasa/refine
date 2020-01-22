#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

serveCSM -batch cube.csm

ref boostrap cube.egads

ref adapt cube-vol.meshb -g cube.egads -x cube01.meshb -f cube01-prop.tec

