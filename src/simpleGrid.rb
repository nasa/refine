#!/usr/bin/env ruby
#
# $Id$
#
# make 2 cell grids for Eric

require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('GridMetric').build

require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMath/GridMath'
require 'GridMetric/GridMetric'

class Grid
 include GridMetric
end

dz = 0.4
grid = Grid.new(10,10,10,10)
grid.addNode(0,0,0)
grid.addNode(1,0,0)
grid.addNode(0,1,0)
grid.addNode(0.3,0.3,dz)
grid.addNode(0.3,0.3,-dz)
grid.addCell(0,1,2,3)
grid.addCell(0,2,1,4)
faceId = 1
grid.addFace(1,3,2,faceId)
grid.addFace(0,2,3,faceId)
grid.addFace(0,3,1,faceId)
faceId = 1
grid.addFace(0,4,2,faceId)
grid.addFace(2,4,1,faceId)
grid.addFace(1,4,0,faceId)

puts 'volume'+grid.minVolume.to_s
puts "original boundary is not right handed" if !grid.rightHandedBoundary

grid.exportFAST
