#! /bin/bash -xue

PACKAGE='ESP'
VERSION='118'

if [ $# -gt 0 ] ; then
   . common.sh  $1
else
   . common.sh
fi

echo Build ${PACKAGE} ${VERSION}

# https://acdl.mit.edu/ESP/PreBuilts/
# https://acdl.mit.edu/ESP/archive/

rm -f ESP${VERSION}Lin.tgz
wget https://acdl.mit.edu/ESP/archive/ESP${VERSION}Lin.tgz
mkdir ${MODULE_DEST}
tar xzf ESP${VERSION}Lin.tgz -C ${MODULE_DEST} --strip-components=1

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
