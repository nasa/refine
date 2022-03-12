
#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

tecplot=-t
egads="-g narrow-cyl.egads"
mapbc="--fun3d-mapbc narrow-cyl-vol.mapbc"
sweeps="-s 10"

function adapt_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref adapt ${inproj}.meshb ${egads} \
	  --spalding 0.01 1000 \
	  ${mapbc} \
	  -x ${outproj}.meshb -x ${outproj}.plt \
	  ${sweeps} ${tecplot}

    cp ref_gather_movie.tec ${inproj}-movie.tec
}

serveCSM -skipTess narrow-cyl.csm

${src}/ref bootstrap narrow-cyl.egads

adapt_cycle narrow-cyl-vol cycle01

