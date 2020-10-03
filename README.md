# Description

`refine` is a 3D mesh adaptation framework implemented in
the C language.

# Quick Start Compile from Git Repo and Basic Usage

`refine` can function without depencies, but the typical use cases of
parallel execution and geometry evaluation require an MPI implementation
and [Engineering Sketch Pad](https://acdl.mit.edu/ESP/ESPreadme.txt) (ESP).
A native implementaion of a recursive coordinate bisection partition
algorithm is included, but better results are expected with
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview).
Initial mesh generation assumes
[TetGen](http://wias-berlin.de/software/tetgen/) or
[AFLR](http://www.simcenter.msstate.edu/research/cavs_cfd/aflr.php) is in
the shell path.
Configuration and compliation is supported with Autoconf and CMake.

## Automake 1.7 (or later) and Autoconf 2.53 (or later):
```
 ./bootstrap
 mkdir -p build
 cd build
 ../configure --prefix=`pwd` \
   --with-mpi=/mpi/path \
   --with-parmetis=/parmetis/path \
   --with-EGADS=/egads/path \
   --with-OpenCASCADE=/opencascade/path
 make
 make install
```
See the INSTALL file for further autoconf build instructions.

## CMake 3.0 (or later):
```
 mkdir -p build
 cd build
 cmake .. -DCMAKE_INSTALL_PREFIX=`pwd` \
   -DCMAKE_PREFIX_PATH="/mpi/path;/parmetis/path;/egads/path;/opencascade/path"
 make
 make install
```

## Usage

The installed `bin`-directory will include the `ref` executable.
Invoking `ref` with no arguments will list available subcommands.
Help on a particular subcommand is available via a `-h`, i.e.,
`ref adapt -h`. If MPI is provided, `refmpi` will allow for parallel
execution. If ESP is provided, `ref` and `refmpifull` includes
EGADS built with OpenCASCADE and `refmpi` includes EGADSlite.

# Examples

The following examples assume that `ref` is in your shell path.
`mpiexec ... refmpi` or `mpiexec ... rempifull` can be substituted for
`ref` in each of these examples if MPI and/or ESP is configured. The
[.meshb and .solb file extensions](https://github.com/LoicMarechal/libMeshb)
are used generically. Other formats are supported, e.g.,
AFLR `*.ugrid`.

## Bootstrapping Mesh Adaptation on an EGADS Model

An `.egads` file can be dumped from OpenCSM in the ESP package.
`ref bootstrap project.egads`
or
`mpiexec ... refmpifull bootstrap project.egads`
which assume that `tetgen` is in your shell path or
`aflr3` is in your shell path with `--mesher aflr` option.
A `project-vol.meshb` is output that includes the surface mesh,
volume mesh, mesh-to-geometry associtivity, and EGADSlite data.

## Mesh Adaptation

The mesh is adapted with
`ref adapt input.meshb -x output.meshb [-m metric.solb] [-g geometry.egads]`
or
`mpiexec ... refmpi adapt input.meshb -x output.meshb [-m metric.solb]`
where a surface curvature metric is used if the `-m` argument is not present.

## Multiscale Metric for Control of Interpolation Error in Lp-norm

In conjunction with the
[Unstructured Grid Adaptation Working Group](https://ugawg.github.io/),
an implementation of the multiscale metric is provided.
`ref multiscale input.meshb scalar.solb complexity output-metric.solb`
or
`mpiexec ... refmpi multiscale input.meshb scalar.solb complexity output-metric.solb`

## Field Interpolation
The fields in a .solb file paired with a donor mesh can be interpolated to
a receptor mesh. This utility can be executed in serial or parallel.

`ref interp donor-mesh.ext donor-field.solb receptor-mesh.ext receptor-field.solb`
or 
`mpiexec ... refmpi interp donor-mesh.ext donor-field.solb receptor-mesh.ext receptor-field.solb`
where the output is `receptor-field.solb`.

