# set variables to be used by all module builders
# expects PACKAGE and VERSION to be set

GROUP='fun3d_users'

GCC_VERSION='6.2.0'
OPENMPI_VERSION='2.1.1'
INTEL_VERSION='2017.2.174'
PARMETIS_VERSION='4.0.3'
ESP_VERSION='119-beta.2021.05.05.1543'

GCC_MODULE="gcc_${GCC_VERSION}"
INTEL_MODULE="intel.${INTEL_VERSION}"
OPENMPI_MODULE="openmpi_${OPENMPI_VERSION}_intel_2017"
ESP_MODULE="ESP/${ESP_VERSION}"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# -gt 0 ] && [ $1 == "jenkinsbuild" ] ; then
  echo "install modules in local directory"
  PREFIX="$script_dir/_modules_jenkinsbuild/${GROUP}"
  mkdir -p $PREFIX/modules
  mkdir -p $PREFIX/modulefiles
else
  PREFIX="/ump/fldmd/home/casb-shared/${USER}/${GROUP}"     # where everything is anchored
fi


MODULE_ROOT="${PREFIX}/modules"         # where the built artifacts reside
MODFILE_ROOT="${PREFIX}/modulefiles"    # where the modulefiles reside

# artifacts
MODULE_BASE="${MODULE_ROOT}/${PACKAGE}"
MODULE_DEST="${MODULE_BASE}/${VERSION}"

# module system file
MODFILE_BASE="${MODFILE_ROOT}/${PACKAGE}"
MODFILE_DEST="${MODFILE_BASE}/${VERSION}"

. /usr/local/pkgs/modules/init/bash
module use ${MODFILE_ROOT}

skipbuild(){
   if [ -e $2 ] ; then
       echo "$2 exists."
       echo "===Exit $1 build.=="
       exit 
   fi
}
