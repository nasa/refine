#! /bin/bash -xue

PACKAGE='ParMETIS'
VERSION='4.0.3'

if [ $# -gt 0 ] ; then
   . common.sh  $1
else
   . common.sh
fi

if test ${PARMETIS_VERSION} != ${VERSION} ; then 
  echo WARNING: common.sh PARMETIS_VERSION should be updated from ${PARMETIS_VERSION} to ${VERSION}
fi

echo Override defaults
LONGVERSION="${VERSION}-${MPT_MODULE}-${GCC_MODULE}"
MODULE_DEST="${MODULE_BASE}/${LONGVERSION}"
MODFILE_DEST="${MODFILE_BASE}/${LONGVERSION}"

skipbuild $PACKAGE $MODFILE_DEST

echo Grabbing ${PACKAGE} ${VERSION}
rm -rf parmetis-${VERSION}*
module load git_2.10.1
git lfs clone git@gitlab.larc.nasa.gov:fun3d-developers/fun3d-dependencies.git
mv fun3d-dependencies/parmetis/parmetis-${VERSION}.tar.gz .
rm -rf fun3d-dependencies

tar zxf parmetis-${VERSION}.tar.gz
cd parmetis-$VERSION

echo Building and Installing ${PACKAGE} ${LONGVERSION}

module purge
module load ${GCC_MODULE}
module load ${MPT_MODULE}
module list

make config prefix=${MODULE_DEST}
make install
cd metis
  make config prefix=${MODULE_DEST}
  make install
cd ..

echo Creating module file
mkdir -p ${MODFILE_BASE}
cat > ${MODFILE_DEST} << EOF
#%Module#
proc ModulesHelp { } { puts stderr "${PACKAGE} ${LONGVERSION}." }

set sys      [uname sysname]
set modname  [module-info name]
set modmode  [module-info mode]

set base    ${MODULE_BASE}
set version ${LONGVERSION}

prereq ${INTEL_MODULE}
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
prepend-path C_INCLUDE_PATH \$base/\$version/include
prepend-path MANPATH \$base/\$version/man
EOF

chgrp -R ${GROUP}  ${MODFILE_BASE} ${MODULE_BASE}
chmod -R g+rX      ${MODFILE_BASE} ${MODULE_BASE}
chmod -R g-w,o-rwx ${MODFILE_BASE} ${MODULE_BASE}
