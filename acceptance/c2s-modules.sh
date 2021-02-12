
module use --append /u/shared/fun3d/fun3d_users/modulefiles

module load FUN3D/13.6

module load gcc_6.2.0
module load ESP/119-beta.2020.09.29.1106
module load tetgen/1.5.0

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/opt/hpe/hpc/mpt/mpt-2.19"

export parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.17r14-gcc_6.2.0"

export egads_path="${module_path}/ESP/119-beta.2020.09.29.1106/EngSketchPad"
export opencascade_path="${module_path}/ESP/119-beta.2020.09.29.1106/OpenCASCADE-7.3.1"

