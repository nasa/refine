#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridMetric adj.c grid.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridCAD FAKEGeom adj.c grid.c gridStruct.h master_header.h'

require 'test/unit'
require 'Grid/Grid'
require 'GridMetric/GridMetric'
require 'GridCAD/GridCAD'

class Grid
 include GridMetric
 include GridCAD
 def totalVolume
  vol = 0.0
  ncell.times { |cellId| vol += volume(cell(cellId)) }
  vol
 end
end

class TestSampleUnit < Test::Unit::TestCase

 def testEdgeProjection
  assert_not_nil              grid = Grid.new(3,0,0,2)
  assert_equal 0,             grid.addNode(0.0,0.0,0.0)
  assert_equal 1,             grid.addNode(0.5,0.1,0.1)
  assert_nil                  grid.projectNodeToEdge(1,1)
  assert_equal grid,          grid.addEdge(0,1,1,0.0,0.55)
  assert_equal grid,          grid.projectNodeToEdge(1,1)
  assert_equal [0.5,0.0,0.0], grid.nodeXYZ(1)
  assert_equal 0.0,           grid.nodeT(0,1)
  assert_equal 0.5,           grid.nodeT(1,1)
 end

 def testFaceProjection
  assert_not_nil              grid = Grid.new(4,0,2,0)
  assert_equal 0,             grid.addNode(0.0,0.0,0.1)
  assert_equal 1,             grid.addNode(1.0,0.0,0.0)
  assert_equal 2,             grid.addNode(0.0,1.0,0.0)
  assert_equal 3,             grid.addNode(1.0,1.0,0.0)
  assert_nil                  grid.projectNodeToFace(0,1)
  assert_equal grid,          grid.addFaceUV(0, 10.1, 20.1,
					     3, 11.0, 21.0,
					     1, 11.0, 20.0,
					     1)
  assert_equal grid,          grid.addFaceUV(0, 10.1, 20.1,
					     3, 11.0, 21.0,
					     2, 10.0, 21.0,
					     1)
  assert_equal grid,          grid.projectNodeToFace(0,1)
  assert_equal [0.0,0.0,0.0], grid.nodeXYZ(0)
  assert_equal [10.0,20.0],   grid.nodeUV(0,1)
 end

 def testProjectionToGeomerty
  assert_not_nil grid = Grid.new(5,1,1,1)
  assert_equal 0, grid.addNode(5.0,5.0,5.0)
  assert_equal 1, grid.addNode(0.0,0.1,0.1)
  assert_equal 2, grid.addNode(1.0,0.1,0.1)
  assert_equal 3, grid.addNode(0.0,1.0,0.1)
  assert_equal 4, grid.addNode(0.0,0.0,1.0)
  assert_equal grid, grid.addCell(1,2,3,4)
  assert_equal grid, grid.addFace(1,2,3,10)
  assert_equal grid, grid.addEdge(1,2,20,0.0,1.0)
  assert_equal grid, grid.setNGeomNode(1)
  5.times do |i| 
   assert_equal grid, grid.safeProjectNode(i), "could not project node #{i}"
  end
  assert_equal [5.0,5.0,5.0], grid.nodeXYZ(0)
  assert_equal [0.0,0.0,0.0], grid.nodeXYZ(1)
  assert_equal [1.0,0.0,0.0], grid.nodeXYZ(2)
  assert_equal [0.0,1.0,0.0], grid.nodeXYZ(3)
  assert_equal [0.0,0.0,1.0], grid.nodeXYZ(4)
  assert_equal 0.0,           grid.nodeT(1,20)
  assert_equal 1.0,           grid.nodeT(2,20)
  assert_equal [10.0,20.0],   grid.nodeUV(1,10)
  assert_equal [11.0,20.0],   grid.nodeUV(2,10)
  assert_equal [10.0,21.0],   grid.nodeUV(3,10)
 end

 def testSafeProjectionAndPositiveVolume
  assert_not_nil grid = Grid.new(4,1,1,1)
  assert_equal 0, grid.addNode(0.0,0.0,-0.5)
  assert_equal 1, grid.addNode(1.0,0.0,0.0)
  assert_equal 2, grid.addNode(0.0,1.0,0.0)
  assert_equal 3, grid.addNode(0.0,0.0,-0.1)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal grid, grid.addFace(0,1,2,10)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s  
  assert_nil   grid.safeProjectNode(0)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_equal [0.0,0.0,-0.5], grid.nodeXYZ(0)
  assert_equal grid, grid.addEdge(0,1,20,0.0,1.0)
  assert_nil   grid.safeProjectNode(0)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_equal [0.0,0.0,-0.5], grid.nodeXYZ(0)
 end

 def isoTet(pert = 0.0)
  grid = Grid.new(4,1,1,1)
  grid.addNode( pert,  0.000, 0.000 )
  grid.addNode( 1.000, 0.000, 0.000 )
  grid.addNode( 0.500, 0.866, 0.000 )
  grid.addNode( 0.500, 0.289, 0.823 ) 
  grid.addCell(0,1,2,3)
  grid.addFaceUV(0,10.0+pert,20.0,
		 1,11.0,20.0,
		 2,10.5,20.866,
		 10)
 end

 def testIsotropicTet
  assert_in_delta 1.000, isoTet.minAR, 1.0e-4
  assert_in_delta 0.975, isoTet(-0.2).minAR, 1.0e-4
 end

 def testOptimizeUVDispacement
  assert_not_nil grid = isoTet(-0.2)
  assert_equal grid, grid.optimizeUV(0,[1.0,0.0])
  assert_in_delta 0.999, grid.minAR, 1.0e-3
  assert_in_delta 0.999, grid.ar([0,1,2,3]), 1.0e-3
 end

 def testSmoothSurf
  assert_not_nil grid = isoTet(-0.2)
  assert_equal grid, grid.smoothNode(0)
  assert_in_delta 0.999, grid.minAR, 1.0e-3
 end

 def testSmooth
  assert_not_nil grid = isoTet(-0.2)
  assert_equal grid, grid.smooth
  assert_in_delta 0.999, grid.minAR, 1.0e-3
 end

end
