#! /bin/bash -xue

PACKAGE='refine-MPT'
VERSION="$(git describe --always --tag | tr -d '\n')"
. ./test-module-string.sh
VERSION="$(module-string ${VERSION})"

if [ $# -gt 0 ] ; then
   . common.sh  $1
else
   . common.sh
fi

TOPDIR='../../..'
PARMETIS="ParMETIS/${PARMETIS_VERSION}_mpt-${MPT_VERSION}_ifort-${INTEL_VERSION}"
ZOLTAN="Zoltan/${ZOLTAN_VERSION}_mpt-${MPT_VERSION}_ifort-${INTEL_VERSION}"
ESP="ESP/${ESP_VERSION}"

echo Build ${PACKAGE} ${VERSION}

module purge
module load ${INTEL_MODULE}
module load ${MPT_MODULE}
module list

( cd ${TOPDIR} && ./bootstrap )

rm -rf   _build_$VERSION
mkdir -p _build_$VERSION
cd       _build_$VERSION

../${TOPDIR}/configure \
  --prefix=${MODULE_DEST} \
  --with-zoltan=${MODULE_ROOT}/${ZOLTAN} \
  --with-metis=${MODULE_ROOT}/${PARMETIS} \
  --with-parmetis=${MODULE_ROOT}/${PARMETIS} \
  --with-EGADS=${MODULE_ROOT}/${ESP}/EngSketchPad \
  --enable-lite \
  CC=icc \
  CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,49,1419 -ftrapuv' \
  LIBS=-lmpi

make -j 12 
make install

echo Copy build log
cp  config.log ${MODULE_DEST}/.

echo Copy and make exec test util
cp  src/*_test ${MODULE_DEST}/bin
chmod g+x ${MODULE_DEST}/bin/*_test

echo Make modulefile
mkdir -p ${MODFILE_BASE}
cat > ${MODFILE_DEST} << EOF
#%Module#
proc ModulesHelp { } { puts stderr "$PACKAGE $VERSION" }

set sys      [uname sysname]
set modname  [module-info name]
set modmode  [module-info mode]

set base    $MODULE_BASE
set version $VERSION

prereq ${MPT_MODULE}

set logr "/bin"

if { \$modmode == "switch1" } {
  set modmode "switchfrom"
}
if { \$modmode == "switch2" } {
  set modmode "switchto"
}
if { \$modmode != "switch3" } {
  system  "\$logr/logger -p local2.info envmodule \$modmode \$modname"
}

prepend-path PATH \$base/\$version/bin

setenv MPI_ADJUST_ALLREDUCE 1

EOF

echo Set group ownership and permssions
chgrp -R ${GROUP}  ${MODULE_DEST}
chmod -R g+rX      ${MODULE_DEST}
chmod -R g-w,o-rwx ${MODULE_DEST}

chgrp -R ${GROUP}  ${MODFILE_DEST}
chmod -R g+rX      ${MODFILE_DEST}
chmod -R g-w,o-rwx ${MODFILE_DEST}
