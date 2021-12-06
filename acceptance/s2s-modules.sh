
module use --append /u/shared/fun3d/fun3d_users/modulefiles
module use --append /u/shared/fun3d/fun3d_users/test_modulefiles

module load --auto FUN3D_INTG
module load ESP/120
module load tetgen

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/opt/hpe/hpc/mpt/mpt-2.23"

export parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.23-intel_2018.3.222"

export egads_path="${module_path}/ESP/120/EngSketchPad"
export opencascade_path="${module_path}/ESP/120/OpenCASCADE"

