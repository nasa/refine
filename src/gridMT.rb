#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h'

require 'test/unit'
require 'Grid/Grid'

class Grid
 def totalVolume
  vol = 0.0
  ncell.times { |cellId| vol += volume(cell(cellId)) }
  vol
 end
end

class TestSampleUnit < Test::Unit::TestCase

 def set_up
  @grid = Grid.new(4,1,0,0)
 end

 def testSafeAlloc
  grid = Grid.new(0,0,0,0)
  assert_equal 1, grid.maxnode
  assert_equal 1, grid.maxcell
  assert_equal 1, grid.maxface
  assert_equal 1, grid.maxedge
 end

  def testCreateGrid
  assert_equal 4, @grid.maxnode
  assert_equal 0, @grid.nnode
  assert_equal 1, @grid.maxcell
  assert_equal 0, @grid.ncell
  assert_equal 1, @grid.maxface
  assert_equal 0, @grid.nface
  assert_equal 1, @grid.maxedge
  assert_equal 0, @grid.nedge
 end

 def testAddCellDegAndCell
  assert_equal 0, @grid.ncell
  assert_equal nil, @grid.cell(0)
  assert_equal nil, @grid.cell(5)
  assert_equal @grid, @grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], @grid.cell(0)
  assert_equal 1, @grid.ncell
  (0..3).each { |n| assert_equal 1, @grid.cellDegree(n)}
 end

 def testAddCellRegFailure
  grid = Grid.new(3,1,0,0)
  assert_equal nil, grid.addCell(0,1,2,3)
  grid = Grid.new(4,1,0,0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal nil, grid.addCell(0,1,2,3)
 end
 
 def testRemoveCell
  assert_not_nil     grid = Grid.new(4,2,0,0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_nil         grid.removeCell(-1)
  assert_nil         grid.removeCell(1)
  assert_nil         grid.removeCell(25625)
  assert_equal 1,    grid.ncell
  assert_equal grid, grid.removeCell(0)
  assert_nil         grid.cell(0)
  assert_nil         grid.removeCell(0)
  assert_equal 0,    grid.ncell
  (0..3).each { |n| assert_equal 0, grid.cellDegree(n)}
 end

 def testReplaceCell
  grid = Grid.new(8,2,0,0)
  assert_equal grid, grid.addCell(0,1,2,3).addCell(4,5,6,7)
  assert_equal grid, grid.removeCell(0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], grid.cell(0)
  assert_equal [4, 5, 6, 7], grid.cell(1)
 end

 def testGetGem
  grid = Grid.new(5,3,0,0)
  assert_equal grid, grid.addCell(3,4,0,1).addCell(3,4,1,2).addCell(3,4,2,0)
  assert_equal [], grid.gem(5,6)
  assert_equal [0], grid.gem(0,1)
  assert_equal [1], grid.gem(1,2)
  assert_equal [2], grid.gem(0,2)
  assert_equal [2,0], grid.gem(3,0)
  assert_equal [2,1,0], grid.gem(3,4)
 end
 
  def testOrient
  assert_equal nil, @grid.orient(0,1,2,3,4,5)
  
  assert_equal [0, 1, 2, 3], @grid.orient(0,1,2,3,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(0,3,1,2,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(0,2,3,1,0,1)

  assert_equal [0, 1, 2, 3], @grid.orient(1,0,3,2,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(1,2,0,3,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(1,3,2,0,0,1)

  assert_equal [0, 1, 2, 3], @grid.orient(2,3,0,1,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(2,1,3,0,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(2,0,1,3,0,1)

  assert_equal [0, 1, 2, 3], @grid.orient(3,2,1,0,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(3,0,2,1,0,1)
  assert_equal [0, 1, 2, 3], @grid.orient(3,1,0,2,0,1)
 end

 def testEquator
  grid = Grid.new(6,4,0,0)
  assert_equal grid, grid.
   addCell(4,5,1,2).addCell(4,5,2,3).addCell(4,5,3,0).addCell(4,5,0,1)
  assert_equal [3,2,1,0], grid.gem(4,5)
  assert_equal 2, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(5)
  assert_equal [], grid.equator(0,2)
  assert_equal [], grid.equator(6,7)
  assert_equal [1,2,3,0,1], grid.equator(4,5)
 end

 def testEquatorGapInMiddle
  grid = Grid.new(6,3,0,0)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,2,3)
  assert_equal [1,0], grid.gem(4,5)
  assert_equal [1,2,3,1], grid.equator(4,5)
 end

 def testEquatorGapInEnd
  grid = Grid.new(6,3,0,0)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,3,1)
  assert_equal [3,1,2,3], grid.equator(4,5)
 end

 def testEquatorTwoGaps
  grid = Grid.new(6,3,0,0)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,3,0)
  assert_equal nil, grid.equator(4,5)
 end

 def testAddNode
  grid = Grid.new(1,0,0,0)
  assert_nil      grid.nodeXYZ(0)
  assert_equal 0, grid.addNode(1.0,2.0,3.0)
  assert_nil      grid.addNode(1.0,0.0,0.0)
  assert_equal [1.0,2.0,3.0], grid.nodeXYZ(0)
 end

 def testMetrics
  assert_equal @grid, @grid.
   addCell( @grid.addNode(0.0,0.0,0.0), @grid.addNode(1.0,0.0,0.0), 
	    @grid.addNode(0.0,1.0,0.0), @grid.addNode(0.0,0.0,1.0) )
  nodes = [0,1,2,3]
  assert_in_delta 1.0/6.0, @grid.volume(nodes), 1.0e-15
  assert_in_delta 1.0/6.0, @grid.minVolume, 1.0e-15
  assert_in_delta 0.732050807568877, @grid.ar(nodes), 1.0e-15
  assert_in_delta 0.732050807568877, @grid.minAR, 1.0e-15
 end

 def testNumberOfFaces
  assert_not_nil  grid = Grid.new(4,1,2,0)
  assert_equal 0, grid.nface 
  assert_equal 2, grid.maxface 
 end

 def testAddAndFindFace
  assert_not_nil     grid = Grid.new(4,1,2,0)
  assert_equal grid, grid.addFace(0, 1, 2, 10)
  assert_equal 0,    grid.findFace(0,1,2)
  assert_nil         grid.findFace(3,1,2)
 end

 def testAddAndRemoveFace
  assert_not_nil     grid = Grid.new(4,1,2,0)
  assert_nil         grid.removeFace(0)
  assert_nil         grid.removeFace(1)
  assert_equal grid, grid.addFace(0, 1, 2, 10)
  assert_nil         grid.removeFace(-1)
  assert_nil         grid.removeFace(1)
  assert_equal grid, grid.addFace(3, 1, 2, 11)
  assert_equal 2,    grid.nface 
  assert_nil         grid.addFace(0, 1, 3, 12)
  assert_equal 2,    grid.nface 
  assert_nil         grid.removeFace(3)
  assert_equal 2,    grid.nface 
 end

 def testFaceId
  assert_not_nil     grid = Grid.new(4,1,2,0)

  assert_nil         grid.faceId( 1, 2, 3 )

  assert_equal grid, grid.addFace(0, 1, 2, 10)
  assert_equal 10,   grid.faceId( 0, 1, 2 )
  assert_equal 10,   grid.faceId( 1, 2, 0 )
  assert_equal 10,   grid.faceId( 2, 0, 1 )
  assert_equal 10,   grid.faceId( 2, 1, 0 )
  assert_nil         grid.faceId( 1, 2, 3 )

  assert_equal grid, grid.addFace(3, 1, 2, 11)
  assert_equal 10,   grid.faceId( 0, 1, 2 )
  assert_equal 11,   grid.faceId( 1, 2, 3 )
 end

 def testNodeUV
  assert_not_nil     grid = Grid.new(4,1,2,0)
  assert_equal grid, grid.addFaceUV(0,20.0,120.0,
				    1,21.0,121.0,
				    2,22.0,122.0,2)
  assert_equal grid, grid.addFaceUV(0,30.0,130.0,
				    1,31.0,131.0,
				    3,33.0,133.0,3)
  assert_nil                 grid.nodeUV(2,3)
  assert_equal [20.0,120.0], grid.nodeUV(0,2)
  assert_equal [21.0,121.0], grid.nodeUV(1,2)
  assert_equal [22.0,122.0], grid.nodeUV(2,2)
  assert_equal [30.0,130.0], grid.nodeUV(0,3)
  assert_equal [31.0,131.0], grid.nodeUV(1,3)
  assert_equal [33.0,133.0], grid.nodeUV(3,3)
 end

 def testSplitEdge4
  assert_not_nil grid = gemGrid
  initalVolume = grid.totalVolume
  assert_equal grid, grid.splitEdge(0,1)
  assert_equal 7, grid.nnode
  assert_equal 8, grid.ncell
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
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

#test for not enough mem for swap and split

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

 def testNumberOfEdges
  assert_not_nil  grid = Grid.new(0,0,0,2)
  assert_equal 0, grid.nedge
  assert_equal 2, grid.maxedge
 end

 def testAddAndFindEdge
  assert_not_nil     grid = Grid.new(4,0,0,2)
  assert_nil         grid.findEdge(0,1)
  assert_equal grid, grid.addEdge(0, 1, 10, 0.0, 0.0)
  assert_equal 1,    grid.nedge
  assert_equal 0,    grid.findEdge(0,1)
 end

 def testAddAndRemoveEdge
  assert_not_nil     grid = Grid.new(4,0,0,2)
  assert_nil         grid.removeEdge(-1)
  assert_nil         grid.removeEdge(0)
  assert_nil         grid.removeEdge(1)
  assert_equal grid, grid.addEdge(0, 1, 10, 0.0, 0.0)
  assert_equal 10,   grid.edgeId(1, 0)
  assert_equal grid, grid.removeEdge(0)
  assert_nil         grid.edgeId(1, 0)
  assert_equal grid, grid.addEdge(3, 1, 11, 0.0, 0.0)
  assert_equal grid, grid.addEdge(0, 2, 12, 0.0, 0.0)
  assert_equal 11,   grid.edgeId(3, 1)
  assert_equal 2,    grid.nedge
  assert_nil         grid.addEdge(1, 2, 13, 0.0, 0.0)
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

 def testGetGeomCurve
  assert_not_nil     grid = Grid.new(4,0,0,4)
  assert_equal grid, grid.addEdge(0, 1, 10, 10.0, 11.0)
  assert_equal grid, grid.addEdge(1, 2, 11, 1.0, 2.0)
  assert_equal grid, grid.addEdge(2, 3, 11, 2.0, 3.0)
  assert_equal grid, grid.addEdge(3, 1, 12, 23.0, 21.0)
  assert_equal 0,         grid.geomCurveSize(15,0)
  assert_equal 2,         grid.geomCurveSize(10,0)
  assert_equal 2,         grid.geomCurveSize(10,1)
  assert_equal 3,         grid.geomCurveSize(11,1)
  assert_equal 3,         grid.geomCurveSize(11,3)
  assert_nil              grid.geomCurve(15,0)
  assert_equal [0, 1],    grid.geomCurve(10,0)
  assert_equal [1, 0],    grid.geomCurve(10,1)
  assert_equal [1, 2, 3], grid.geomCurve(11,1)
  assert_equal [3, 2, 1], grid.geomCurve(11,3)
  assert_equal [10.0, 11.0],    grid.geomCurveT(10,0)
  assert_equal [11.0, 10.0],    grid.geomCurveT(10,1)
  assert_equal [1.0, 2.0, 3.0], grid.geomCurveT(11,1)
 end

 def testFindCellWithFace
  assert_not_nil     grid = Grid.new(5,2,0,4)
  assert_equal grid, grid.addCell(0,1,4,3)
  assert_nil         grid.findCellWithFace(34)
  assert_nil         grid.findCellWithFace(0)
  assert_equal grid, grid.addFace(0,1,2,11)
  assert_nil         grid.findCellWithFace(0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal 1,    grid.findCellWithFace(0)
 end

 def testRightHandedFaces
  assert_not_nil grid = Grid.new(4,1,2,0)
  assert_equal grid, grid.
   addCell( grid.addNode(0.0,0.0,0.0), grid.addNode(1.0,0.0,0.0), 
	    grid.addNode(0.0,1.0,0.0), grid.addNode(0.0,0.0,1.0) )
  assert_equal grid,  grid.addFace(0,1,2,11)
  assert_equal true,  grid.rightHandedFace(0)
  assert_equal grid,  grid.addFace(0,2,3,11)
  assert_equal true,  grid.rightHandedBoundary
  assert_equal grid,  grid.removeFace(grid.findFace(0,1,2))
  assert_equal grid,  grid.addFace(0,2,1,11)
  assert_equal false, grid.rightHandedFace(0)
  assert_equal false, grid.rightHandedBoundary
 end
 
 # make register unique

 # allocating a new chunk of nodes, faces, cells

end
