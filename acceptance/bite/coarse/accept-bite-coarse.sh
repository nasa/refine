#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

serveCSM -batch bite.csm \
	 > accept-bite-coarse-csm.txt
ref bootstrap bite.egads --global 10000 10000 90 -t \
    > accept-bite-coarse-bootstrap.txt

exit 0

ref adapt bite-vol.meshb -g bite.egads -f bite-final.tec \
    > accept-bite-coarse-status.txt

cat accept-bite-coarse-status.txt
../../check.rb accept-bite-coarse-status.txt 0.3 3.0



