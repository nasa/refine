#! /bin/bash -xue

PACKAGE='ESP'
VERSION='116'

if [ $# -gt 0 ] ; then
   . common.sh  $1
else
   . common.sh
fi

echo Build ${PACKAGE} ${VERSION}

# https://acdl.mit.edu/ESP/PreBuilts/
# https://acdl.mit.edu/ESP/archive/

echo wget https://acdl.mit.edu/ESP/PreBuilts/ESP${VERSION}Lin.tgz

echo ${MODULE_DEST}

exit

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

prepend-path PATH \$base/\$version/EngSketchPad/bin
prepend-path LD_LIBRARY_PATH \$base/\$version/EngSketchPad/lib
prepend-path LD_LIBRARY_PATH \$base/\$version/OpenCASCADE-6.8.1/Linux/lib

EOF

echo Set group ownership and permssions
chgrp -R ${GROUP}  ${MODULE_DEST}
chmod -R g+rX      ${MODULE_DEST}
chmod -R g-w,o-rwx ${MODULE_DEST}

chgrp -R ${GROUP}  ${MODFILE_DEST}
chmod -R g+rX      ${MODFILE_DEST}
chmod -R g-w,o-rwx ${MODFILE_DEST}
