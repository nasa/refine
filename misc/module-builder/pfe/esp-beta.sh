#! /bin/bash -xue

PACKAGE='ESP'

rm -f ESPbeta.tgz
wget -N https://acdl.mit.edu/ESP/archive/ESPbeta.tgz
raw=$(stat -c %Y ESPbeta.tgz)
timestamp=$(date -d @${raw} +"%Y.%m.%d.%H%M")
VERSION="119-beta.${timestamp}"

if [ $# -gt 0 ] ; then
   . common.sh  $1
else
   . common.sh
fi

OCC_COPY_SOURCE="${MODULE_BASE}/117/OpenCASCADE-7.3.1"
OCC_COPY_DEST=${MODULE_DEST}/OpenCASCADE-7.3.1

echo Build ${PACKAGE} ${VERSION}

module purge

mkdir ${MODULE_DEST}
cp -r ${OCC_COPY_SOURCE} ${MODULE_DEST}

rm -rf EngSketchPad
tar xzf ESPbeta.tgz
( cd EngSketchPad/config && ./makeEnv ${OCC_COPY_DEST} )
( cd EngSketchPad/src && . ../ESPenv.sh && make CC=gcc CXX=g++ FCOMP=gfortran)

mkdir ${MODULE_DEST}/EngSketchPad
cp -r EngSketchPad/include EngSketchPad/lib EngSketchPad/bin ${MODULE_DEST}/EngSketchPad

mkdir -p ${MODFILE_BASE}
cat > ${MODFILE_DEST} << EOF
#%Module#
proc ModulesHelp { } { puts stderr "$PACKAGE $VERSION" }

set sys      [uname sysname]
set modname  [module-info name]
set modmode  [module-info mode]

set base    $MODULE_BASE
set version $VERSION

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

setenv ESP_ROOT \$base/\$version/EngSketchPad
setenv CASROOT \$base/\$version/OpenCASCADE-7.3.1

prepend-path PATH \$base/\$version/EngSketchPad/bin

prepend-path LD_LIBRARY_PATH \$base/\$version/EngSketchPad/lib
prepend-path LD_LIBRARY_PATH \$base/\$version/OpenCASCADE-7.3.1/lib

EOF

echo Set group ownership and permssions
chgrp -R ${GROUP}  ${MODULE_DEST}
chmod -R g+rX      ${MODULE_DEST}
chmod -R g-w,o-rwx ${MODULE_DEST}

chgrp -R ${GROUP}  ${MODFILE_DEST}
chmod -R g+rX      ${MODFILE_DEST}
chmod -R g-w,o-rwx ${MODFILE_DEST}
