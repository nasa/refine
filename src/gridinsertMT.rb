#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridMetric adj.c grid.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridInsert adj.c grid.h gridmetric.h gridStruct.h master_header.h'

require 'test/unit'
require 'Grid/Grid'
require 'GridMetric/GridMetric'
require 'GridInsert/GridInsert'

class Grid
 include GridMetric
 include GridInsert
end

class TestGridInsert < Test::Unit::TestCase

 def rightTet
  grid = Grid.new(10,12,0,0)
  grid.addCell( 
	       grid.addNode(0.0,0.0,0.0), 
	       grid.addNode(1.0,0.0,0.0), 
	       grid.addNode(0.0,1.0,0.0), 
	       grid.addNode(0.0,0.0,1.0) )
 end

 def gemGrid(nequ=4, a=nil, dent=nil, x0 = nil, gap = nil)
  a  = a  || 0.1
  x0 = x0 || 1.0
  grid = Grid.new(nequ+2+1,14,14,14)
  n = Array.new
  n.push grid.addNode(  x0,0.0,0.0)
  n.push grid.addNode(-1.0,0.0,0.0)
  nequ.times do |i| 
   angle = 2.0*Math::PI*(i-1)/(nequ)
   s = if (dent==i) then 0.9 else 1.0 end
   n.push grid.addNode(0.0,s*a*Math.sin(angle),s*a*Math.cos(angle)) 
  end
  n.push 2
  ngem = (gap)?(nequ-1):(nequ)
  ngem.times do |i|
   grid.addCell(n[0],n[1],n[i+2],n[i+3])
  end
  grid  
 end

 def testSplitEdge4
  assert_not_nil grid = gemGrid
  assert_equal grid, grid.splitEdge(0,1)
  assert_equal 7, grid.nnode
  assert_equal 8, grid.ncell
 end

 def testSplitEdge4onSameBC
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal grid, grid.addFaceUV(0,0.0,10.0,
				    1,1.0,11.0,
				    2,2.0,12.0,11)
  assert_equal grid, grid.addFaceUV(1,1.0,11.0,
				    0,0.0,10.0,
				    5,5.0,15.0,11)
  assert grid.rightHandedBoundary, "original boundary is not right handed"
  assert_equal grid, grid.splitEdge(0,1)
  assert_nil         grid.faceId(0,1,2)
  assert_nil         grid.faceId(0,1,5) 
  assert_equal 11,   grid.faceId(0,6,2)
  assert_equal 11,   grid.faceId(6,1,2)
  assert_equal 11,   grid.faceId(0,6,5) 
  assert_equal 11,   grid.faceId(6,1,5) 
  assert grid.rightHandedBoundary, "split boundary is not right handed"
  assert_equal [0.0,10.0], grid.nodeUV(0,11)
  assert_equal [1.0,11.0], grid.nodeUV(1,11)
  assert_equal [2.0,12.0], grid.nodeUV(2,11)
  assert_equal [5.0,15.0], grid.nodeUV(5,11)
  assert_equal [0.5,10.5], grid.nodeUV(6,11)
 end

 def testSplitEdge4onDifferentBC
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal grid, grid.addFaceUV(0,20.0,120.0,
				    1,21.0,121.0,
				    2,22.0,122.0,2)
  assert_equal grid, grid.addFaceUV(1,51.0,151.0,
				    0,50.0,150.0,
				    5,55.0,155.0,5)
  assert grid.rightHandedBoundary, "original boundary is not right handed"
  assert_equal grid, grid.splitEdge(0,1)
  assert_nil         grid.faceId(0,1,2)
  assert_nil         grid.faceId(0,1,5) 
  assert_equal 2,    grid.faceId(0,6,2)
  assert_equal 2,    grid.faceId(6,1,2)
  assert_equal 5,    grid.faceId(0,6,5) 
  assert_equal 5,    grid.faceId(6,1,5) 
  assert grid.rightHandedBoundary, "split boundary is not right handed"
  assert_equal [20.0,120.0], grid.nodeUV(0,2)
  assert_equal [50.0,150.0], grid.nodeUV(0,5)
  assert_equal [21.0,121.0], grid.nodeUV(1,2)
  assert_equal [51.0,151.0], grid.nodeUV(1,5)
  assert_equal [22.0,122.0], grid.nodeUV(2,2)
  assert_equal [55.0,155.0], grid.nodeUV(5,5)
  assert_equal [20.5,120.5], grid.nodeUV(6,2)
  assert_equal [50.5,150.5], grid.nodeUV(6,5)
 end

 def testSplitGeometryEdge4
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal grid, grid.addFace(0,1,2,2)
  assert_equal grid, grid.addFace(1,0,5,5)
  assert grid.rightHandedBoundary, "original boundary is not right handed"
  assert_equal grid, grid.addEdge(0,1,15, 0.0, 1.0)
  assert_equal grid, grid.splitEdge(0,1)
  assert_nil         grid.edgeId(0,1)
  assert_equal 15,   grid.edgeId(0,6)
  assert_equal 15,   grid.edgeId(6,1)
  assert_equal [0.0, 0.5, 1.0], grid.geomCurveT(15,0)
  assert_equal [1.0, 0.5, 0.0], grid.geomCurveT(15,1)
 end

 def testSplitWithoutGeometryEdge4
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal grid, grid.addFace(0,1,2,11)
  assert_equal grid, grid.addFace(1,0,5,11)
  assert grid.rightHandedBoundary, "original boundary is not right handed"
  assert_equal grid, grid.splitEdge(0,1)
  assert_nil         grid.edgeId(0,6)
  assert_nil         grid.edgeId(6,1)
 end

 def testAdaptToSpacing1
  assert_not_nil     grid = rightTet
  assert_equal grid, grid.resetSpacing
  assert_equal grid, grid.scaleSpacing(0,0.4)
  assert_equal grid, grid.scaleSpacing(1,0.4)
  assert_equal grid, grid.adapt
  assert_equal 2, grid.ncell
 end

 def testAdaptToSpacing3
  assert_not_nil     grid = rightTet
  assert_equal grid, grid.resetSpacing
  assert_equal grid, grid.scaleSpacing(0,0.2)
  assert_equal grid, grid.scaleSpacing(1,0.6)
  assert_equal grid, grid.scaleSpacing(2,0.6)
  assert_equal grid, grid.scaleSpacing(3,0.6)
  assert_equal grid, grid.adapt
  assert_equal 4, grid.ncell
 end

end
