#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Adj master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.h gridStruct.h master_header.h'

require 'test/unit'
require 'Adj/Adj'
require 'Grid/Grid'

class TestGrid < Test::Unit::TestCase

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

 def testAddCellAndCellDegree
  assert_equal 0, @grid.ncell
  assert_equal nil, @grid.cell(-1)
  assert_equal nil, @grid.cell(0)
  assert_equal nil, @grid.cell(5)
  assert_equal @grid, @grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], @grid.cell(0)
  assert_equal 1, @grid.ncell
  (0..3).each { |n| assert_equal 1, @grid.cellDegree(n)}
 end

 def testAddCellFailure
  assert_not_nil     grid = Grid.new(3,1,0,0)
  assert_nil         grid.addCell(0,1,2,3)
  assert_not_nil     grid = Grid.new(4,1,0,0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_nil         grid.addCell(0,1,2,3)
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

 def testInitNodeFrozenState
  @grid.addNode(0,0,0)
  @grid.addNode(1,0,0)
  @grid.removeNode(0)
  assert_equal true, @grid.nodeFrozen(0)
  assert_equal false, @grid.nodeFrozen(1)
  assert_equal true, @grid.nodeFrozen(2)
  assert_equal true, @grid.nodeFrozen(@grid.maxnode)
 end

 def testNodeFrozenState
  @grid.addNode(0,0,0)
  @grid.addNode(1,0,0)
  assert_equal 0,     @grid.nfrozen
  assert_equal false, @grid.nodeFrozen(0)
  assert_equal false, @grid.nodeFrozen(1)
  assert_nil          @grid.freezeNode(@grid.maxnode)
  assert_equal @grid, @grid.freezeNode(0)
  assert_equal 1,     @grid.nfrozen
  assert_equal true,  @grid.nodeFrozen(0)
  assert_equal false, @grid.nodeFrozen(1)
  assert_nil          @grid.thawNode(@grid.maxnode)
  assert_equal @grid, @grid.thawNode(0)
  assert_equal 0,     @grid.nfrozen
  assert_equal false, @grid.nodeFrozen(0)
  assert_equal false, @grid.nodeFrozen(1)
 end

 def testGlobalFreezeState
  grid = Grid.new(2,0,0,0)
  assert_equal 0,     grid.addNode(0,0,0)
  assert_equal 1,     grid.addNode(1,0,0)
  assert_equal grid,  grid.freezeAll
  assert_equal true,  grid.nodeFrozen(0)
  assert_equal true,  grid.nodeFrozen(1)
  assert_equal grid,  grid.thawAll
  assert_equal false, grid.nodeFrozen(0)
  assert_equal false, grid.nodeFrozen(1)
 end

 def testReplaceCell
  grid = Grid.new(8,2,0,0)
  assert_equal grid, grid.addCell(0,1,2,3).addCell(4,5,6,7)
  assert_equal grid, grid.removeCell(0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], grid.cell(0)
  assert_equal [4, 5, 6, 7], grid.cell(1)
 end

 def testReconnectCell
  assert_not_nil             grid = Grid.new(8,3,0,0)
  assert_equal grid,         grid.addCell(0,1,2,3).addCell(3,4,5,6)
  assert_equal [0, 1, 2, 3], grid.cell(0)
  assert_equal [3, 4, 5, 6], grid.cell(1)
  assert_equal 2,            grid.cellDegree(3)
  assert_equal 0,            grid.cellDegree(7)
  assert_nil                 grid.reconnectCell(-1,-1)
  assert_nil                 grid.reconnectCell(8,8)
  assert_equal grid,         grid.reconnectCell(3,7)
  assert_equal [0, 1, 2, 7], grid.cell(0)
  assert_equal [7, 4, 5, 6], grid.cell(1)
  assert_equal 0,            grid.cellDegree(3)
  assert_equal 2,            grid.cellDegree(7)
 end

 def testCellConnectivity
  assert_not_nil      grid = Grid.new(5,2,0,0)
  assert_equal grid,  grid.addCell(0,1,2,3)
  assert_equal grid,  grid.addCell(0,1,2,4)

  assert_equal true,  grid.cellEdge(0,1)
  assert_equal false, grid.cellEdge(3,4)
  assert_equal false, grid.cellEdge(-1,4)
  assert_equal false, grid.cellEdge(3,-1)

  assert_equal true,  grid.cellFace(0,1,2)
  assert_equal false, grid.cellFace(0,3,4)
  assert_equal false, grid.cellFace(-1,0,0)
  assert_equal false, grid.cellFace(0,-1,0)
  assert_equal false, grid.cellFace(0,0,-1)
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
  assert_not_nil          grid = Grid.new(6,3,0,0)
  assert_equal grid,      grid.addCell(4,5,1,2).addCell(4,5,3,1)
  assert_equal [3,1,2,3], grid.equator(4,5)
 end

 def testEquatorTwoGaps
  assert_not_nil     grid = Grid.new(6,3,0,0)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,3,0)
  assert_nil         grid.equator(4,5)
 end

 def testAddNode
  assert_not_nil              grid = Grid.new(1,0,0,0)
  assert_nil                  grid.nodeXYZ(0)
  assert_equal 0,             grid.addNode(1.0,2.0,3.0)
  assert_nil                  grid.addNode(1.0,0.0,0.0)
  assert_equal [1.0,2.0,3.0], grid.nodeXYZ(0)
 end

 def testSetNodeXYZ
  assert_not_nil                 grid = Grid.new(1,0,0,0)
  assert_equal 0,                grid.addNode(1.0,2.0,3.0)
  assert_equal [1.0,2.0,3.0],    grid.nodeXYZ(0)
  assert_equal grid,             grid.setNodeXYZ(0,[10.0,20.0,30.0])
  assert_equal [10.0,20.0,30.0], grid.nodeXYZ(0)
 end

 def testAddAndRemoveNode
  assert_not_nil      grid = Grid.new(4,0,0,0)
  assert_nil          grid.removeNode(5)
  assert_nil          grid.removeNode(2)
  assert_equal false, grid.validNode(-1)
  assert_equal false, grid.validNode(2)
  assert_equal false, grid.validNode(20)
  assert_equal 0,     grid.addNode(1.0,2.0,3.0)
  assert_equal 1,     grid.addNode(1.1,2.1,3.1)
  assert_equal 2,     grid.addNode(1.2,2.2,3.2)
  assert_equal 3,     grid.addNode(1.3,2.3,3.3)
  assert_equal true,  grid.validNode(2)
  assert_equal 4,     grid.nnode
  assert_equal grid,  grid.removeNode(2)
  assert_equal false, grid.validNode(2)
  assert_equal true,  grid.validNode(3)
  assert_not_nil      grid.nodeXYZ(3)
  assert_equal 3,     grid.nnode
  assert_nil          grid.removeNode(2)
  assert_nil          grid.nodeXYZ(2)
  assert_equal 3,     grid.nnode
  assert_equal 2,     grid.addNode(1.2,2.2,3.2)
 end

 def testNumberOfFaces
  assert_not_nil  grid = Grid.new(4,1,2,0)
  assert_equal 0, grid.nface 
  assert_equal 2, grid.maxface 
 end

 def testAddAndFindFace
  assert_not_nil           grid = Grid.new(4,1,2,0)
  assert_nil               grid.face(0)
  assert_equal grid,       grid.addFace(0, 1, 2, 10)
  assert_equal 0,          grid.findFace(0,1,2)
  assert_equal [0,1,2,10], grid.face(0)
  assert_nil               grid.findFace(3,1,2)
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

 def testReconnectFace
  assert_not_nil             grid = Grid.new(6,0,3,0)
  assert_equal grid,         grid.addFace(0,3,1,1)
  assert_equal grid,         grid.addFace(0,1,2,2)
  assert_equal grid,         grid.addFace(2,1,4,2)

  assert_equal [0, 3, 1, 1], grid.face(0)
  assert_equal [0, 1, 2, 2], grid.face(1)
  assert_equal [2, 1, 4, 2], grid.face(2)
  assert_nil                 grid.reconnectFace(2,-1,-1)
  assert_nil                 grid.reconnectFace(2, 8, 8)
  assert_equal grid,         grid.reconnectFace(2, 3, 5)
  assert_equal [0, 3, 1, 1], grid.face(0)
  assert_equal [0, 1, 2, 2], grid.face(1)
  assert_equal [2, 1, 4, 2], grid.face(2)
  assert_equal grid,         grid.reconnectFace(2, 1, 5)
  assert_equal [0, 3, 1, 1], grid.face(0)
  assert_equal [0, 5, 2, 2], grid.face(1)
  assert_equal [2, 5, 4, 2], grid.face(2)

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
  assert_equal grid,         grid.setNodeUV(0,2,8.0,9.0)
  assert_equal [8.0,9.0],    grid.nodeUV(0,2)
 end

 def testNumberOfGeomEdges
  assert_not_nil  grid = Grid.new(0,0,0,2)
  assert_equal 0, grid.nedge
  assert_equal 2, grid.maxedge
 end

 def testAddAndFindGeomEdge
  assert_not_nil         grid = Grid.new(4,0,0,2)
  assert_nil             grid.findEdge(0,1)
  assert_equal grid,     grid.addEdge(0, 1, 10, 0.0, 1.0)
  assert_equal 1,        grid.nedge
  assert_equal 0,        grid.findEdge(0,1)
  assert_nil             grid.edge(-1)
  assert_nil             grid.edge(5)
  assert_equal [0,1,10], grid.edge(0)
 end

 def testGeomEdgeTValues
  assert_not_nil     grid = Grid.new(4,0,0,2)
  assert_nil         grid.nodeT(0,10)
  assert_nil         grid.setNodeT(1,10,1.5) 
  assert_equal grid, grid.addEdge(0, 1, 10, 0.0, 1.0)
  assert_equal 0.0,  grid.nodeT(0,10)
  assert_equal 1.0,  grid.nodeT(1,10)
  assert_equal grid, grid.addEdge(1, 2, 10, 1.0, 2.0)
  assert_equal grid, grid.setNodeT(1,10,1.5) 
  assert_equal 1.5,  grid.nodeT(1,10) 
 end

 def testAddAndRemoveGeomEdge
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

 def testDeleteEdgesWithThawedNode
  assert_not_nil             grid = Grid.new(10,0,0,10)
  assert_equal 0,            grid.addNode(0,0,0)
  assert_equal 1,            grid.addNode(1,0,0)
  assert_equal 2,            grid.addNode(2,0,0)
  assert_equal 3,            grid.addNode(3,0,0)
  assert_equal grid,         grid.addEdge(0,1,1,0.0,1.0)
  assert_equal grid,         grid.addEdge(2,1,1,2.0,1.0)
  assert_equal grid,         grid.addEdge(2,3,1,2.0,3.0)
  assert_equal grid,         grid.addEdge(1,3,99,99.0,99.0)
  assert_equal 4,            grid.nedge
  assert_equal grid,         grid.freezeNode(0)
  assert_equal grid,         grid.freezeNode(1)
  assert_nil                 grid.deleteThawedEdgeSegments(-1)
  assert_nil                 grid.deleteThawedEdgeSegments(0)
  assert_equal grid,         grid.deleteThawedEdgeSegments(1)
  assert_equal 2,            grid.nedge
  assert_equal [0,1,1],  grid.edge(0)
  assert_nil             grid.edge(1)
  assert_nil             grid.edge(2)
  assert_equal [1,3,99], grid.edge(3)
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

 def testNGeomElements
  assert_not_nil grid = Grid.new(1,1,1,1)
  assert_equal 0,  grid.nGeomNode
  assert_equal 0,  grid.nGeomEdge
  assert_equal 0,  grid.nGeomFace
  assert_equal grid, grid.setNGeomNode(1)
  assert_equal grid, grid.setNGeomEdge(2)
  assert_equal grid, grid.setNGeomFace(3)
  assert_equal 1,    grid.nGeomNode
  assert_equal 2,    grid.nGeomEdge
  assert_equal 3,    grid.nGeomFace
 end

 def testGeometryNode
  assert_not_nil grid = Grid.new(3,0,0,0)
  assert_equal false, grid.geometryNode(-1)
  assert_equal false, grid.geometryNode(0)
  assert_equal grid,  grid.setNGeomNode(2)
  assert_equal 2,     grid.nGeomNode
  assert_equal false, grid.geometryNode(-1)
  assert_equal true,  grid.geometryNode(0)
  assert_equal true,  grid.geometryNode(1)
  assert_equal false, grid.geometryNode(2)
  assert_equal false, grid.geometryNode(grid.maxnode)
 end

 def testGeometryEdge
  assert_not_nil grid = Grid.new(3,0,0,1)
  assert_equal false, grid.geometryEdge(0)
  assert_equal grid,  grid.addEdge(0,1,10,0.0,1.0)
  assert_equal false, grid.geometryEdge(-1)
  assert_equal true,  grid.geometryEdge(0)
  assert_equal false, grid.geometryEdge(grid.maxnode)
 end

 def testGeometryFace
  assert_not_nil grid = Grid.new(3,0,0,1)
  assert_equal false, grid.geometryFace(0)
  assert_equal grid,  grid.addFace(0,1,2,10)
  assert_equal false, grid.geometryFace(-1)
  assert_equal true,  grid.geometryFace(0)
  assert_equal false, grid.geometryFace(grid.maxnode)
 end

 def testPack
  assert_not_nil grid = Grid.new(5,2,2,2)
  assert_equal 0, grid.addNode(9.0,0.0,9.0)
  assert_equal 1, grid.addNode(0.0,0.0,0.0)
  assert_equal grid, grid.addCell( 
				  0, 
				  grid.addNode(1.0,0.0,0.0), 
				  grid.addNode(0.0,1.0,0.0), 
				  grid.addNode(0.0,0.0,1.0) )
  assert_equal [0,2,3,4], grid.cell(0)
  assert_equal grid, grid.addCell(1,2,3,4)
  assert_equal [1,2,3,4], grid.cell(1)
  assert_equal [1], grid.gem(1,2)
  assert_equal grid, grid.addFace(2,3,4,1)
  assert_equal grid, grid.addFaceUV(1,1.0,11.0,
				    2,2.0,12.0,
				    3,3.0,13.0,
				    1)
  assert_equal [1.0,11.0], grid.nodeUV(1,1)
  assert_equal [2.0,12.0], grid.nodeUV(2,1)
  assert_equal [3.0,13.0], grid.nodeUV(3,1)
  assert_equal grid, grid.addEdge(2,3,1,22.0,23.0)
  assert_equal grid, grid.addEdge(1,2,1,21.0,22.0)
  assert_equal 21.0, grid.nodeT(1,1)
  assert_equal 22.0, grid.nodeT(2,1)
  assert_equal grid, grid.removeCell(0)
  assert_equal grid, grid.removeFace(0)
  assert_equal grid, grid.removeEdge(0)
  assert_equal grid, grid.removeNode(0)
  assert_equal grid, grid.freezeNode(1)

  assert_equal grid, grid.pack
  assert_equal grid, grid.pack
  assert_equal true, grid.nodeFrozen(0)
  assert_equal false, grid.nodeFrozen(1)
  assert_equal 1, grid.ncell
  assert_equal [0,1,2,3], grid.cell(0)
  assert_equal [0], grid.gem(0,1)
  assert_equal [1.0,11.0], grid.nodeUV(0,1)
  assert_equal [2.0,12.0], grid.nodeUV(1,1)
  assert_equal [3.0,13.0], grid.nodeUV(2,1)
  assert_equal 21.0, grid.nodeT(0,1)
  assert_equal 22.0, grid.nodeT(1,1)
  assert_equal 4, grid.addNode(9.0,0.0,9.0)
  assert_equal grid, grid.addCell(0,2,3,4)
  assert_equal [0,2,3,4], grid.cell(1)
  assert_equal grid, grid.addFace(1,2,3,2)
  assert_equal 1, grid.findFace(1,2,3)
  assert_equal grid, grid.addEdge(2,3,1,22.0,23.0)
  assert_equal 2, grid.maxedge
  assert_equal 1, grid.findEdge(2,3)
  assert_equal 1, grid.findEdge(3,2)
  assert_equal 0, grid.findEdge(0,1)
  assert_equal 0, grid.findEdge(1,0)

  assert_equal grid, grid.pack

  assert_equal grid, grid.removeCell(1)
  assert_equal grid, grid.removeFace(1)
  assert_equal grid, grid.removeEdge(1)
  assert_equal grid, grid.removeNode(4)

  assert_equal grid, grid.pack

 end

 def testRetreveGeomEdgeAndStoredEndPoints
  assert_not_nil          grid = Grid.new(4,1,1,2)
  assert_nil              grid.addGeomEdge(1,0,1)
  assert_equal grid,      grid.setNGeomEdge(1)
  assert_equal 1,         grid.nGeomEdge
  assert_equal( -1,       grid.geomEdgeStart(1))
  assert_equal( -1,       grid.geomEdgeEnd(1))
  assert_equal grid,      grid.addGeomEdge(1,0,1)
  assert_nil              grid.addGeomEdge(2,0,1)
  assert_equal grid,      grid.addEdge(0,2,1,0.0,2.0)
  assert_equal grid,      grid.addEdge(2,1,1,2.0,1.0)
  assert_equal 3,         grid.geomEdgeSize(1)
  assert_equal grid,      grid.addGeomEdge(1,1,0)
  assert_equal [1, 2, 0], grid.geomEdge(1)
  assert_equal( -1,        grid.geomEdgeStart(-1))
  assert_equal( -1,        grid.geomEdgeEnd(-1))
  assert_equal( -1,        grid.geomEdgeStart(0))
  assert_equal( -1,        grid.geomEdgeEnd(0))
  assert_equal(  1,        grid.geomEdgeStart(1))
  assert_equal(  0,        grid.geomEdgeEnd(1))
  assert_equal( -1,        grid.geomEdgeStart(2))
  assert_equal( -1,        grid.geomEdgeEnd(2))
 end

 def testFindFrozenEdgeEndPoints
  assert_not_nil             grid = Grid.new(10,0,0,10)
  assert_equal 0,            grid.addNode(0,0,0)
  assert_equal 1,            grid.addNode(1,0,0)
  assert_equal 2,            grid.addNode(2,0,0)
  assert_equal 3,            grid.addNode(3,0,0)
  assert_equal grid,         grid.setNGeomEdge(1)
  assert_equal grid,         grid.addGeomEdge(1,0,3)
  assert_equal grid,         grid.addEdge(0,1,1,0.0,1.0)
  assert_equal grid,         grid.addEdge(2,1,1,2.0,1.0)
  assert_equal grid,         grid.addEdge(2,3,1,2.0,3.0)
  assert_equal [0, 1, 2, 3], grid.geomEdge(1)

  assert_equal( -1,          grid.frozenEdgeEndPoint(2,0) )
  assert_equal( -1,          grid.frozenEdgeEndPoint(1,4) )

  assert_equal   0,          grid.frozenEdgeEndPoint(1,0)
  assert_equal   3,          grid.frozenEdgeEndPoint(1,3)

  assert_equal grid,         grid.freezeNode(0)
  assert_equal   0,          grid.frozenEdgeEndPoint(1,0)
  assert_equal   3,          grid.frozenEdgeEndPoint(1,3)

  assert_equal grid,         grid.freezeNode(1)
  assert_equal   1,          grid.frozenEdgeEndPoint(1,0)
  assert_equal   3,          grid.frozenEdgeEndPoint(1,3)

  assert_equal grid,         grid.freezeNode(2)
  assert_equal   2,          grid.frozenEdgeEndPoint(1,0)
  assert_equal   3,          grid.frozenEdgeEndPoint(1,3)

  assert_equal grid,         grid.freezeNode(3)
  assert_equal   3,          grid.frozenEdgeEndPoint(1,0)
  assert_equal   0,          grid.frozenEdgeEndPoint(1,3)

  assert_equal grid,         grid.thawNode(0)
  assert_equal   0,          grid.frozenEdgeEndPoint(1,0)
  assert_equal   1,          grid.frozenEdgeEndPoint(1,3)

 end

 def testSortNodesToGridExStandard
  assert_not_nil          grid = Grid.new(6,3,3,2)
  assert_equal 0,         grid.addNode( 1, 0, 0)
  assert_equal 1,         grid.addNode(-1, 0, 0)
  assert_equal 2,         grid.addNode( 0, 1, 0)
  assert_equal 3,         grid.addNode( 0, 0, 1)
  assert_equal 4,         grid.addNode( 0, 0, 0)
  assert_equal 5,         grid.addNode( 1, 0, 0)
  assert_equal [0,0,1],   grid.nodeXYZ(3)
  assert_equal grid,      grid.freezeNode(3)
  assert_equal grid,      grid.addCell( 1, 4, 2, 3)
  assert_equal grid,      grid.addCell( 0, 2, 4, 3)
  assert_equal grid,      grid.addCell( 0, 4, 5, 3)
  assert_equal [0,4,5,3], grid.cell(2)
  assert_equal grid,      grid.addFace( 1, 4, 2, 2)
  assert_equal grid,      grid.addFace( 0, 2, 4, 2)
  assert_equal grid,      grid.addFace( 0, 4, 5, 1)
  assert_equal 2,         grid.findFace(0, 4, 5)
  assert_equal grid,      grid.addEdge( 0, 4, 1, 0.0, 4.0)
  assert_equal 0,         grid.findEdge(0,4)
  assert_equal grid,      grid.addEdge( 1, 4, 1, 1.0, 4.0)
  assert_equal 1,         grid.findEdge(1,4)
  assert_equal grid,      grid.setNGeomNode(2)
  assert_equal grid,      grid.setNGeomEdge(1)
  assert_equal grid,      grid.addGeomEdge(1,0,1)
  assert_equal [0,4,1],   grid.geomEdge(1)

  assert_equal grid,      grid.sortNodeGridEx

  assert_equal [0,0,1],   grid.nodeXYZ(5)
  assert_equal false,     grid.nodeFrozen(3)
  assert_equal true,      grid.nodeFrozen(5)
  assert_equal [0,2,3,5], grid.cell(2)
  assert_equal [2,1,0],   grid.gem(2,5)
  assert_equal 0,         grid.findFace(0, 2, 3)
  assert_equal 0,         grid.findEdge(0,2)
  assert_equal 1,         grid.findEdge(1,2)
  assert_equal [0,2,1],   grid.geomEdge(1)

  assert_equal grid,      grid.sortNodeGridEx

  assert_equal [0,0,1],   grid.nodeXYZ(5)
  assert_equal [0,2,3,5], grid.cell(2)
  assert_equal [2,1,0],   grid.gem(2,5)
  assert_equal 0,         grid.findFace(0, 2, 3)
  assert_equal 0,         grid.findEdge(0,2)
  assert_equal 1,         grid.findEdge(1,2)
  assert_equal [0,2,1],   grid.geomEdge(1)
 end

 def testTec
  assert_not_nil     grid = Grid.new(3,0,1,0)
  assert_equal grid, grid.addFace( grid.addNode(0,0,0),
				   grid.addNode(1,0,0),
				   grid.addNode(0,1,0),
				   1)
  assert_equal grid, grid.writeTecplotSurfaceZone
 end


 # make register unique

 # allocating a new chunk of nodes, faces, cells

end
