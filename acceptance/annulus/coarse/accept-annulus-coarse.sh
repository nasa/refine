#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

serveCSM -batch annulus.csm \
	 > accept-annulus-coarse-csm.txt
${src}/ref bootstrap annulus.egads --global 10000 10000 90 \
    > accept-annulus-coarse-bootstrap.txt
${src}/ref adapt annulus-vol.meshb -g annulus.egads -f annulus-final.tec \
    > accept-annulus-coarse-status.txt

cat accept-annulus-coarse-status.txt
../../check.rb accept-annulus-coarse-status.txt 0.3 3.0



