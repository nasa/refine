
esp_module=ESP/120-beta.2022.07.13.0553

module use --append /u/shared/fun3d/fun3d_users/modulefiles
module use --append /u/shared/fun3d/fun3d_users/test_modulefiles

module load --auto FUN3D_INTG
module load ${esp_module}
module load tetgen

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/opt/hpe/hpc/mpt/mpt-2.23"

export parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.23-intel_2018.3.222"

export egads_path="${module_path}/${esp_module}/EngSketchPad"
export opencascade_path="${module_path}/${esp_module}/OpenCASCADE"

