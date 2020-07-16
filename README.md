# Description

Refine is a 3D mixed-element grid adaptation framework implemented in
the C language.

# Quick Start

If you checked out the repository with git, 
[requires automake >= 1.7 and autoconf >= 2.53]
```
 ./bootstrap
 mkdir -p build
 cd build
 ../configure --prefix=`pwd`
 make
 make install
```
See the [INSTALL file](https://github.com/nasa/refine/blob/master/INSTALL) for further build instructions. Here are some examples:
with MPI and Zoltan,
```
../configure --with-zoltan=/zoltan/path CC=mpicc FC=mpif90 CFLAGS='-DHAVE_MPI'
```
with EGADS and OpenCASCADE,
```
../configure --with-EGADS=/egads/path --with-OpenCASCADE=/opencascade/path
```
with EGADSlite,
```
../configure --with-EGADS=/egads/path --enable-lite
```
These can be combined, e.g., MPI+Zoltan+EGADSlite.

# Multiscale Metric for Control of Interpolation Error in Lp-norm
In conjunction with the
[Unstructured Grid Adaptation Working Group](https://ugawg.github.io/),
an implementation of the multiscale metric is provided.
```
./build/src/ref_metric_test --lp project.meshb project-mach.solb 2 1.5 3.0e4 project-metric.solb [--kexact] [--hmax max_edge_length]
```
Where,
 - `project.meshb` is the grid in libMeshb format
 - `project-mach.solb` is a scalar field (Mach number) in libMeshb format
 - `2` is the norm order
 - `1.5` is the gradation limit
 - `3.e4` is the complexity (C), where the new mesh will have approximately 2C vertices and 12C tetrahedra
 - `project-metric.solb` is the output metric in libMeshb format
 - use (sequential only) k-exact Hessian reconstruction with `--kexact`,
   otherwise default to L2-projection Hessian reconstruction
 - `max_edge_length` is an optional limit on maximum edge length with an attempt to hold complexity

See [LoicMarechal/libMeshb](https://github.com/LoicMarechal/libMeshb)
for details on the libMeshb format.

# Metric Conformity Evaluation
In conjunction with the
[Unstructured Grid Adaptation Working Group](https://ugawg.github.io/),
an implementation of descriptive statistics for
edge lengths and mean ratio is provided.
```
./build/src/ref_histogram_test project.meshb project-metric.solb
```
Where,
 - `project.meshb` is the grid in libMeshb format
 - `project-metric.solb` is the metric in libMeshb format

See [LoicMarechal/libMeshb](https://github.com/LoicMarechal/libMeshb)
for details on the libMeshb format.
Histograms are written as `ref_histogram_quality.tec` and
`ref_histogram_ratio.tec`, where each of the two the columns are
mean ratio/edge length and normalized count.  

# Scalar Field for Interpolation Error Verification
To create an initial unit cube domain
```
./build/src/ref_acceptance 1 initial.meshb
```
To compute the scalar on a domain,
```
./build/src/ref_acceptance -u tanh3 project.meshb project-scalar.solb
```
Where,
 - `tanh3` can be `sinfun3`, `sinatan3`, or `tanh3`
 - `project.meshb` is the grid in libMeshb format
 - `project-scalar.solb` is a scalar field in libMeshb format

# Interpolation Error Evaluation
The norm of interpolation error of a scalar function on a "candidate" grid can
be computed based on a scalar function on "truth" grid.
The candidate solution is interpolated to the truth grid,
assuming the solution is linear in each element, and
a norm of the difference of the interpolated candidate solution and the
truth solution is integrated on the truth grid.
```
./build/src/ref_interp_test --error truth.meshb truth-scalar.solb candidate.meshb candidate-scalar.solb 2
```
Where,
 - `truth.meshb` is the truth grid in libMeshb format
 - `truth-scalar.solb` is truth scalar field in libMeshb format
 - `candidate.meshb` is the candidate grid in libMeshb format
 - `candidate-scalar.solb` is candidate scalar field in libMeshb format
 - `2` is the power of the norm

The output is two numbers: the characteristic edge length
(number of vertices raised to the -1/3 power) and the interpolation error norm.

# Field Interpolation
The fields in a .solb file paired with a donor mesh can be interpolated to
a receptor mesh. This utility can be executed in serial or parallel.
```
./build/src/ref_interp_test --field donor-mesh.ext donor-field.solb receptor-mesh.ext receptor-field.solb
```
Where,
 - `donor-mesh.meshb` is the donor mesh in binary AFLR (.lb8.ugrid,.b8.ugrid) or libMeshb (.meshb) format
 - `donor-field.solb` is the donor field(s) in libMeshb format
 - `receptor-mesh.ext` is the receptor mesh in binary AFLR (.lb8.ugrid,.b8.ugrid) or libMeshb (.meshb) format
 - `receptor-field.solb` is the receptor field(s) in libMeshb format

The output is `receptor-field.solb`.

# Grid Adaptation
```
./build/src/ref_driver -i project.meshb -m project-metric.solb [-g project.egads] -x output.meshb
```
Where,
 - `-i project.meshb` is the grid in libMeshb format
 - `-m project-metric.solb` is the metric in libMeshb format
 - `-g project.egads` is optional geometry (when compiled with EGADS)
 - `-x output.meshb` is the output grid name
