# Description

Refine is a 3D mixed-element grid adaptation framework implemented in
the C language.

# Quick Start

If you checked out the repository with git, 
[requires automake >= 1.7 and autoconf >= 2.53]

 % ./bootstrap
 % mkdir -p build
 % cd build
 % ../configure --prefix=`pwd`
 % make
 % make install

See the INSTALL file for further build instructions.

# Metric for Control of Interpolation Error in Lp-norm
In conjunction with the
[Unstructured Grid Adaptation Working Group](https://ugawg.github.io/),
an implementation of the Lp-norm metric is provided.
```
./build/src/ref_metric_test --lp project.meshb project-mach.solb 2 1.5 3.0e4 project-metric.solb
```
Where:
 - `project.meshb` is the grid in libMeshb format
 - `project-mach.solb` is a scalar field (Mach number) in libMeshb format
 - `2` is the norm order
 - `1.5` is the gradation limit
 - `3.e4` is the complexity
 - `project-metric.solb` is the output metric in libMeshb format
 See [LoicMarechal/libMeshb](https://github.com/LoicMarechal/libMeshb)
 for details on the libMeshb format.

# Metric Conformity Evaluation
In conjunction with the
[Unstructured Grid Adaptation Working Group](https://ugawg.github.io/),
an implementation of descriptive statistics for edge lengths and mean ratio is
provided.
```
./build/src/ref_histogram_test project.meshb project-metric.solb
```
 - `project.meshb` is the grid in libMeshb format
 - `project-metric.solb` is the metric in libMeshb format
 See [LoicMarechal/libMeshb](https://github.com/LoicMarechal/libMeshb) for details on the libMeshb format.
Histograms are written as `ref_histogram_quality.tec` and
`ref_histogram_ratio.tec`, where each of the two the columns are
mean ratio/edge length and normalized count.  
