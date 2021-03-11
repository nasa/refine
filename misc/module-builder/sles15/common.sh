# set variables to be used by all module builders
# expects PACKAGE and VERSION to be set

GROUP='n1337'

MPT_VERSION='mpt.2.23'
INTEL_VERSION='2018.3.222'
ESP_VERSION='119-beta.2021.03.05.1211'

INTEL_MODULE="comp-intel/${INTEL_VERSION}"

PARMETIS="/swbuild/fun3d/shared/fun3d_users/sles15/test_modules/ParMETIS-64/4.0.3_mpt-2.23_ifort-2018.3.222"
ZOLTAN="/swbuild/fun3d/shared/fun3d_users/sles15/test_modules/Zoltan/3.82_mpt-2.23_ifort-2018.3.222"
ESP="/swbuild/fun3d/shared/fun3d_users/modules/ESP/119-beta.2021.03.05.1211"
MPI="/nasa/hpe/mpt/2.23_sles15_patch11654"

PREFIX="/swbuild/fun3d/shared/fun3d_users/sles15" # where everything is anchored

MODULE_ROOT="${PREFIX}/test_modules"         # where the built artifacts reside
MODFILE_ROOT="${PREFIX}/test_modulefiles"    # where the modulefiles reside

# artifacts
MODULE_BASE="${MODULE_ROOT}/${PACKAGE}"
MODULE_DEST="${MODULE_BASE}/${VERSION}"

# module system file
MODFILE_BASE="${MODFILE_ROOT}/${PACKAGE}"
MODFILE_DEST="${MODFILE_BASE}/${VERSION}"

. /usr/share/modules/init/bash
module use ${MODFILE_ROOT}

