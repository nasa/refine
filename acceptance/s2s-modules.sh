
module use --append /u/shared/fun3d/fun3d_users/modulefiles

module load gcc_6.2.0
module load ESP/119-beta.2021.03.05.1211
module load tetgen

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/opt/hpe/hpc/mpt/mpt-2.19"

export parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.17r14-gcc_6.2.0"

export egads_path="${module_path}/ESP/119-beta.2021.03.05.1211/EngSketchPad"
export opencascade_path="${module_path}/ESP/119-beta.2021.03.05.1211/OpenCASCADE-7.3.1"

