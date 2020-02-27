#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

cp ../geom/cube-wire.csm .

serveCSM -batch cube-wire.csm \
	 > accept-cube-wire-coarse-csm.txt
ref bootstrap cube-wire.egads --global 10000 10000 90 \
    > accept-cube-wire-coarse-bootstrap.txt
ref adapt cube-wire-vol.meshb -g cube-wire.egads -f cube-wire-final.tec \
    > accept-cube-wire-coarse-status.txt

cat accept-cube-wire-coarse-status.txt
../../check.rb accept-cube-wire-coarse-status.txt 0.3 3.0



