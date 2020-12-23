#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

# serveCSM -batch cylcyl.csm

${src}/ref boostrap cylcyl.egads

${src}/ref_geom_test --enrich2 cylcyl-adapt-surf.meshb cylcyl.egads

${src}/ref adapt cylcyl-vol.meshb \
      -g cylcyl.egads \
      --surrogate ref_geom_enrich2.meshb \
      -x cylcyl-curve.meshb \
      -f cylcyl-curve.tec \
      -s 10

${src}/ref adapt cylcyl-curve.meshb \
      -g cylcyl.egads \
      --implied-complexity 4000 \
      --surrogate ref_geom_enrich2.meshb \
      -x cylcyl-fine.meshb \
      -f cylcyl-fine.tec \
      -s 10

