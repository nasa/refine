#!/usr/bin/env ruby

starting_path = Dir.pwd

refine_path = "#{File.dirname $0}/../src"

$:.push refine_path

Dir.chdir refine_path

require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Grid').build
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'

Dir.chdir starting_path

mesh3d_filename = ARGV[0]
fast_filename = ARGV[1] || mesh3d_filename.sub(/\.mesh3D$/,'.fgrid')

grid = Grid.from_mesh3d  File.expand_path(mesh3d_filename)
grid.exportFAST File.expand_path(grid_filename)

__END__

#Number of Verts
nVerts
nVert xyzs
#Number of Elements
nTets nPyramids nPrisms nHexas
nTet indices (4 per element -- 1 bias)
nPyramid indices (5 per element -- 1 bias)
nPrism indices (6 per element -- 1 bias)
nHexa indices (8 per element -- 1 bias)
#Number of Boundaries, volume element flag
nBounds eflag
nBound times:
 Boundary Name
 Boundary type
 Geom index
 nTris nQuads
 nTri indices (3 per facet -- 1 bias)
 nQuad indices (4 per facet -- 1 bias)

The "Geom index" will contain some indication of the originating
geometry (for adaptation and higher-order vertex fill-in). The
"Boundary type" should be something like "Wall", "InFlow", "OutFlow",
"FarField", "Symmetry" or an index pointing to that...

If the eflag is set to 1 then the lines containing the surface element
indices will contain the number of the attached volume element.
Otherwise the triangle connectivity will contain 3 numbers and the
quadrilaterals will contain 4 numbers.

