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

class TestGridCAD < Test::Unit::TestCase

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
  assert_nil   grid.project
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_equal [0.0,0.0,-0.5], grid.nodeXYZ(0)
  assert_equal grid, grid.addEdge(0,1,20,0.0,1.0)
  assert_nil   grid.safeProjectNode(0)
  assert_nil   grid.project
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_equal [0.0,0.0,-0.5], grid.nodeXYZ(0)
 end

 def isoTet(xpert = 0.0, zpert = 0.0, edge = nil)
  grid = Grid.new(4,1,1,1)
  grid.addNode( xpert, 0.000, 0.000 )
  grid.addNode( 1.000, 0.000, 0.000 )
  grid.addNode( 0.500, 0.866, 0.000 )
  grid.addNode( 0.500, 0.289, 0.823+zpert ) 
  grid.addCell(0,1,2,3)
  grid.addFaceUV(0,10.0+xpert,20.0,
		 1,11.0,20.0,
		 2,10.5,20.866,
		 10)
  grid.addEdge(0, 1, 20, 0.0+xpert, 1.0) if edge
  grid
 end

 def testIsotropicTet
  assert_in_delta 1.000, isoTet.minAR, 1.0e-4
  assert_in_delta 0.975, isoTet(-0.2).minAR, 1.0e-4
 end

 def testOptimizeTDispacement
  assert_not_nil grid = isoTet(-0.2,0.0,true)
  assert_equal grid, grid.optimizeT(0,1.0)
  assert_in_delta 0.999, grid.minAR, 1.0e-3
  assert_in_delta 0.0, grid.nodeT(0,20), 5.0e-2
 end

 def testSmoothEdge
  assert_not_nil grid = isoTet(-0.2,0.0,true)
  assert_equal grid, grid.smoothNode(0)
  assert_in_delta 0.999, grid.minAR, 1.0e-3
 end

 def testOptimizeUVDispacement
  assert_not_nil grid = isoTet(-0.2)
  assert_equal grid, grid.optimizeUV(0,[1.0,0.0])
  assert_in_delta 0.999, grid.minAR, 1.0e-3
  assert_in_delta 10.0, grid.nodeUV(0,10)[0], 5.0e-2
  assert_in_delta 20.0, grid.nodeUV(0,10)[1], 1.0e-15
 end

 def testSmoothSurf
  assert_not_nil grid = isoTet(-0.2)
  assert_equal grid, grid.smoothNode(0)
  assert_in_delta 0.999, grid.minAR, 1.0e-3
 end

 def testOptimizeXYZDispacement
  assert_not_nil grid = isoTet( 0.0, 1.0 )
  assert_equal grid, grid.optimizeXYZ(3,[0.0,0.0,-1.0])
  assert_in_delta 1.00, grid.minAR, 1.0e-2
  assert_equal grid, grid.optimizeXYZ(3,[0.0,0.0,1.0])
  assert_in_delta 1.00, grid.minAR, 1.0e-4
 end

 def testSmoothVol
  assert_not_nil grid = isoTet(0.0,1.0)
  assert_equal grid, grid.smoothNode(3)
  assert_in_delta 1.00, grid.minAR, 1.0e-1
  assert_equal grid, grid.smoothNode(3)
  assert_in_delta 1.00, grid.minAR, 1.0e-3
  assert_equal grid, grid.smoothNode(3)
  assert_in_delta 1.00, grid.minAR, 1.0e-4
 end

 def testSmooth
  assert_not_nil grid = isoTet(-4.0)
  assert_equal grid, grid.smooth
  assert_in_delta 1.0, grid.minAR, 1.0e-3
 end

 def gemGrid(nequ=4)
  grid = Grid.new(nequ+2+1,nequ*2,0,0)
  n = Array.new
  n.push grid.addNode( 1.0,0.0,0.0)
  n.push grid.addNode(-1.0,0.0,0.0)
  nequ.times do |i| 
   angle = 2.0*Math::PI*(i-1)/(nequ)
   n.push grid.addNode(0.0,Math.sin(angle),Math.cos(angle)) 
  end
  n.push 2
  center = grid.addNode(0.2,0.2,0.2)
  nequ.times do |i|
   grid.addCell(n[0],center,n[i+2],n[i+3])
   grid.addCell(center,n[1],n[i+2],n[i+3])
  end
  grid  
 end

 def testSmartLaplacianSmooth
  ngem =6
  assert_not_nil  grid=gemGrid(ngem)
  assert_in_delta 0.38, grid.minAR, 1.0e-3
  assert_equal grid, grid.smartLaplacian(ngem+2)
  assert_in_delta 0.812, grid.minAR, 1.0e-3
  assert_in_delta 0.0, grid.nodeXYZ(ngem+2)[0], 1.0e-15
  assert_in_delta 0.0, grid.nodeXYZ(ngem+2)[1], 1.0e-15
  assert_in_delta 0.0, grid.nodeXYZ(ngem+2)[2], 1.0e-15
 end

end
