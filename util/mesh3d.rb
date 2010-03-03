#!/usr/bin/env ruby

starting_path = Dir.pwd

refine_path = "#{File.dirname $0}/../src"

$:.push refine_path

Dir.chdir refine_path
require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Grid').build
require 'Adj/Adj'
require 'Line/Line'
require 'Sort/Sort'
require 'Grid/Grid'
require 'GridMath/GridMath'

RubyExtensionBuilder.new('GridMetric').build
require 'GridMetric/GridMetric'
class Grid
 include GridMetric
end

Dir.chdir starting_path

root = ARGV[0] || 'full_rounded' 
grid = Grid.from_mesh3d(root+'.mesh3D')

grid.writeVTK(root+'.vtk')
grid.writeTecplotVolumeGeom

