
module use --append /u/shared/fun3d/fun3d_users/modulefiles

esp_module=ESP/121-beta.2022.08.18.1703

module load mpt-2.23

module load gcc_6.2.0
module load ${esp_module}
module load tetgen

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/opt/hpe/hpc/mpt/mpt-2.23"

export parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.23-gcc_6.2.0"

export egads_path="${module_path}/${esp_module}/EngSketchPad"
export opencascade_path="${module_path}/${esp_module}/OpenCASCADE"

