#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Grid').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Sort/Sort'
require 'Grid/Grid'

class TestGrid < Test::Unit::TestCase

 EMPTY = -1

 def set_up
  @grid = Grid.new(4,1,0,0)
 end
 def setup ; set_up ; end

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
  assert_equal 0, @grid.nprism
  assert_equal 0, @grid.npyramid
  assert_equal 0, @grid.nquad
  assert_equal 0, @grid.partId
  assert_equal 0, @grid.globalnnode
  assert_equal 0, @grid.globalncell
  assert_equal 0, @grid.nUnusedNodeGlobal
  assert_equal 0, @grid.nUnusedCellGlobal
 end

 def testSetPartId
  assert_equal @grid, @grid.setPartId(17)
  assert_equal 17,    @grid.partId  
 end

 def testSetGlobalNNode
  assert_equal @grid, @grid.setGlobalNNode(293)
  assert_equal 293,   @grid.globalnnode
 end

 def testSetGlobalNCell
  assert_equal @grid, @grid.setGlobalNCell(12567)
  assert_equal 12567, @grid.globalncell
 end

 def testAddCellAndCellDegree
  assert_equal 0, @grid.ncell
  assert_equal nil, @grid.cell(-1)
  assert_equal nil, @grid.cell(0)
  assert_equal nil, @grid.cell(5)
  assert_equal 0, @grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], @grid.cell(0)
  assert_equal 1, @grid.ncell
  (0..3).each { |n| assert_equal 1, @grid.cellDegree(n)}
 end

 def testCellInitializedNotValid
  assert_nil @grid.cell(0)
 end
 
 def testRemoveCell
  assert_not_nil     grid = Grid.new(4,2,0,0)
  grid.addCell(0,1,2,3)
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

 def testAddingCellAllocatesMaxCell
  assert_not_nil     grid = Grid.new(5,1,0,0)
  assert_equal 0,    grid.addCell(0,1,2,3)
  assert_equal 1,    grid.maxcell
  assert_equal 1,    grid.addCell(0,1,2,4)
  assert_equal 5001, grid.maxcell # chunk size hard coded to 5000
 end
 
 def testGetAndSetCellGlobal
  grid = Grid.new(4,1,0,0)
  4.times { grid.addNode(1.0,2.0,3.0) }
  assert_equal EMPTY, grid.cellGlobal(-1)
  assert_equal EMPTY, grid.cellGlobal(0)
  assert_equal EMPTY, grid.cellGlobal(1)
  assert_nil          grid.setCellGlobal(-1,13)
  assert_nil          grid.setCellGlobal(0,13)
  assert_nil          grid.setCellGlobal(10,13)
  assert_equal 0,     grid.addCell(0,1,2,3)
  assert_equal EMPTY, grid.cellGlobal(0)
  assert_equal grid,  grid.setCellGlobal(0,15)
  assert_equal EMPTY, grid.cellGlobal(-1)
  assert_equal 15,    grid.cellGlobal(0)
  assert_equal EMPTY, grid.cellGlobal(1)
 end

 def testRemoveAndAddTracksGlobalCellId
  grid = Grid.new(5,2,0,0)
  5.times { grid.addNode(1.0,2.0,3.0) }
  grid.addCell(0,1,2,3)
  grid.setGlobalNCell(3512)
  grid.setCellGlobal(0,15)
  grid.removeCell(0)
  assert_equal [15], grid.getUnusedCellGlobal
  grid.addCell(0,1,2,3)
  assert_equal [], grid.getUnusedCellGlobal
  grid.addCell(0,2,1,4)
  assert_equal 15,   grid.cellGlobal(0)
  assert_equal 3512, grid.cellGlobal(1)
  assert_equal 3513, grid.globalncell
 end

 def testJoinUnusedGlobalCellId
  grid = Grid.new(5,2,0,0)
  5.times { grid.addNode(1.0,2.0,3.0) }
  grid.addCell(0,1,2,3)
  grid.setGlobalNCell(3512)
  grid.setCellGlobal(0,15)
  grid.removeCell(0)
  assert_equal 1, grid.nUnusedCellGlobal
  assert_equal [15], grid.getUnusedCellGlobal
  assert_equal grid, grid.joinUnusedCellGlobal(75)
  assert_equal [15,75], grid.getUnusedCellGlobal
  assert_equal grid, grid.joinUnusedCellGlobal(5)
  assert_equal [5,15,75], grid.getUnusedCellGlobal
  assert_equal grid, grid.joinUnusedCellGlobal(5)
  assert_equal grid, grid.joinUnusedCellGlobal(15)
  assert_equal grid, grid.joinUnusedCellGlobal(75)
  assert_equal [5,15,75], grid.getUnusedCellGlobal
  grid.addCell(0,1,2,3)
  assert_equal 5, grid.cellGlobal(0)
  assert_equal [15,75], grid.getUnusedCellGlobal
 end

 def testJoinUnusedGlobalCellIdPreviouslyUsed
  grid = Grid.new(4,3,0,0)
  4.times { grid.addNode(1.0,2.0,3.0) }
  grid.addCell(0,1,2,3)
  grid.setGlobalNCell(4512)
  grid.setCellGlobal(0,15)
  grid.addCell(0,1,2,3)
  grid.addCell(0,1,2,3)
  grid.removeCell(1)
  grid.removeCell(2)
  assert_equal [4512,4513], grid.getUnusedCellGlobal
  grid.addCell(0,1,2,3)
  grid.addCell(0,1,2,3)
  assert_equal [], grid.getUnusedCellGlobal
  grid.removeCell(2)
  assert_equal [4512], grid.getUnusedCellGlobal
  grid.removeCell(1)
  assert_equal [4512,4513], grid.getUnusedCellGlobal
 end

 def testEliminateUnusedGlobalCellId
  grid = Grid.new(5,5,0,0)
  5.times { grid.addNode(1.0,2.0,3.0) }
  grid.addCell(0,1,2,3)
  grid.setGlobalNCell(100)
  grid.setCellGlobal(0,99)
  4.times { grid.addCell(0,1,2,3) }
  grid.removeCell(1)
  grid.removeCell(3)
  assert_equal [100,102], grid.getUnusedCellGlobal
  assert_equal 104, grid.globalncell 
  assert_equal  99, grid.cellGlobal(0)
  assert_equal 101, grid.cellGlobal(2)
  assert_equal 103, grid.cellGlobal(4)

  assert_equal grid, grid.eliminateUnusedCellGlobal

  assert_equal [], grid.getUnusedCellGlobal
  assert_equal 102, grid.globalncell 
  assert_equal  99, grid.cellGlobal(0)
  assert_equal 100, grid.cellGlobal(2)
  assert_equal 101, grid.cellGlobal(4)  
 end

 def testEliminateUnusedGlobalCellIdOutOfOrder
  grid = Grid.new(5,5,0,0)
  5.times { grid.addNode(1.0,2.0,3.0) }
  grid.addCell(0,1,2,3)
  grid.setGlobalNCell(100)
  grid.setCellGlobal(0,99)
  4.times { grid.addCell(0,1,2,3) }
  5.times { |n| grid.setCellGlobal(n,103-n) }
  grid.removeCell(1)
  grid.removeCell(3)
  assert_equal [100,102], grid.getUnusedCellGlobal
  assert_equal 104, grid.globalncell 
  assert_equal 103, grid.cellGlobal(0)
  assert_equal 101, grid.cellGlobal(2)
  assert_equal  99, grid.cellGlobal(4)

  assert_equal grid, grid.eliminateUnusedCellGlobal

  assert_equal [], grid.getUnusedCellGlobal
  assert_equal 102, grid.globalncell 
  assert_equal 101, grid.cellGlobal(0)
  assert_equal 100, grid.cellGlobal(2)
  assert_equal  99, grid.cellGlobal(4)  
 end

 def testEliminateUnusedGlobalCellIdExtrema
  grid = Grid.new(4,3,0,0)
  4.times { grid.addNode(1.0,2.0,3.0) }
  grid.addCell(0,1,2,3)
  grid.setGlobalNCell(100)
  grid.setCellGlobal(0,99)
  2.times { grid.addCell(0,1,2,3) }
  grid.removeCell(0)
  grid.removeCell(2)
  assert_equal [ 99,101], grid.getUnusedCellGlobal
  assert_equal 102, grid.globalncell 
  assert_equal 100, grid.cellGlobal(1)

  assert_equal grid, grid.eliminateUnusedCellGlobal

  assert_equal [], grid.getUnusedCellGlobal
  assert_equal 100, grid.globalncell 
  assert_equal  99, grid.cellGlobal(1)
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

 def testThawedNodeAdd
  grid = Grid.new(2,0,0,0)
  assert_equal 0,     grid.addNode(0,0,0)
  assert_equal 1,     grid.addNode(1,0,0)
  assert_equal grid,  grid.freezeAll
  assert_equal true,  grid.nodeFrozen(0)
  assert_equal true,  grid.nodeFrozen(1)
  assert_equal grid,  grid.removeNode(1)
  assert_equal 1,     grid.addNode(1,0,0)
  assert_equal true,  grid.nodeFrozen(0)
  assert_equal false, grid.nodeFrozen(1)
 end

 def testGetAndSetNAux
  assert_equal 0, @grid.naux
  assert_nil @grid.setAux(0,0,0.5)
  assert_equal 0, @grid.aux(0,0)
  assert_equal @grid, @grid.setNAux(5)
  assert_equal 5, @grid.naux
  assert_equal @grid, @grid.setAux(0,0,0.5)
  assert_equal 0.5, @grid.aux(0,0)
 end

 def testReplaceCell
  grid = Grid.new(8,2,0,0)
  assert_equal 0, grid.addCell(0,1,2,3)
  assert_equal 1, grid.addCell(4,5,6,7)
  assert_equal grid, grid.removeCell(0)
  assert_equal 0, grid.addCell(0,1,3,2)
  assert_equal [0, 1, 3, 2], grid.cell(0)
  assert_equal [4, 5, 6, 7], grid.cell(1)
 end

 def testReconnectCell
  assert_not_nil             grid = Grid.new(8,3,0,0)
  grid.addCell(0,1,2,3)
  grid.addCell(3,4,5,6)
  assert_equal [0, 1, 2, 3], grid.cell(0)
  assert_equal [3, 4, 5, 6], grid.cell(1)
  assert_equal 2,            grid.cellDegree(3)
  assert_equal 0,            grid.cellDegree(7)
  assert_nil                 grid.reconnectAllCell(-1,-1)
  assert_nil                 grid.reconnectAllCell(8,8)
  assert_equal grid,         grid.reconnectAllCell(3,7)
  assert_equal [0, 1, 2, 7], grid.cell(0)
  assert_equal [7, 4, 5, 6], grid.cell(1)
  assert_equal 0,            grid.cellDegree(3)
  assert_equal 2,            grid.cellDegree(7)
 end

 def testCellConnectivity
  assert_not_nil      grid = Grid.new(5,2,0,0)
  grid.addCell(0,1,2,3)
  grid.addCell(0,1,2,4)

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

 def testFindCellWithFace
  assert_not_nil     grid = Grid.new(5,2,0,4)
  grid.addCell(0,1,4,3)
  assert_nil         grid.findCellWithFace(34)
  assert_nil         grid.findCellWithFace(0)
  grid.addFace(0,1,2,11)
  assert_nil         grid.findCellWithFace(0)
  grid.addCell(0,1,2,3)
  assert_equal 1,    grid.findCellWithFace(0)
 end

 def testFindCell
  assert_not_nil     grid = Grid.new(7,2,0,0)
  grid.addCell(0,1,2,3)
  assert_equal EMPTY, grid.findCell([-1,0,10,5])
  assert_equal 0,     grid.findCell([0,1,2,3])
  assert_equal 0,     grid.findCell([3,2,0,1])
  assert_equal EMPTY, grid.findCell([3,4,5,6])
  grid.addCell(3,4,5,6)
  assert_equal 1,     grid.findCell([3,4,5,6])
 end

 def testFindOtherCellWithThreeNodes
  assert_not_nil     grid = Grid.new(5,5,5,5)
  assert_equal(-1,   grid.findOtherCellWith3Nodes(0,1,2,0) )
  grid.addCell(0,1,2,3)
  assert_equal( 0,   grid.findOtherCellWith3Nodes(0,1,2,1) )
  assert_equal( 0,   grid.findOtherCellWith3Nodes(0,1,2,-1) )
  assert_equal(-1,   grid.findOtherCellWith3Nodes(0,1,2,0) )
  grid.addCell(0,1,2,4)
  assert_equal( 1,   grid.findOtherCellWith3Nodes(0,1,2,0) )
  assert_equal( 0,   grid.findOtherCellWith3Nodes(0,1,2,1) )
 end

 def testFindOtherCellWithThreeNodesOutOfRangeInputs
  assert_not_nil     grid = Grid.new(5,5,5,5)
  assert_equal(-1,   grid.findOtherCellWith3Nodes(-1,1,2,0) )
  assert_equal(-1,   grid.findOtherCellWith3Nodes(10,1,2,0) )
  assert_equal(-1,   grid.findOtherCellWith3Nodes(0,-1,2,0) )
  assert_equal(-1,   grid.findOtherCellWith3Nodes(0,10,2,0) )
  assert_equal(-1,   grid.findOtherCellWith3Nodes(0,1,-1,0) )
  assert_equal(-1,   grid.findOtherCellWith3Nodes(0,1,10,0) )
 end

 def testGetGem
  grid = Grid.new(5,3,0,0)
  grid.addCell(3,4,0,1)
  grid.addCell(3,4,1,2)
  grid.addCell(3,4,2,0)
  assert_equal [], grid.gem(5,6)
  assert_equal [0], grid.gem(0,1)
  assert_equal [1], grid.gem(1,2)
  assert_equal [2], grid.gem(0,2)
  assert_equal [2,0], grid.gem(3,0)
  assert_equal [2,1,0], grid.gem(3,4)
 end
 
 def testGemLocality
  grid = Grid.new(5,2,0,0).setPartId(5)
  5.times { |n| grid.addNode(1,2,3); grid.setNodePart(n,5)}
  grid.addCell(0,1,2,3)
  grid.addCell(0,2,1,4)
  grid.gem(0,1)
  assert_equal true, grid.gemIsAllLocal
  grid.setNodePart(4,2)
  assert_equal false, grid.gemIsAllLocal
  grid.setNodePart(4,5)
  assert_equal true, grid.gemIsAllLocal
  grid.setNodePart(3,2)
  assert_equal false, grid.gemIsAllLocal
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
  grid.addCell(4,5,1,2)
  grid.addCell(4,5,2,3)
  grid.addCell(4,5,3,0)
  grid.addCell(4,5,0,1)
  assert_equal [3,2,1,0], grid.gem(4,5)
  assert_equal 2, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(5)
  assert_equal [], grid.equator(0,2)
  assert_equal [], grid.equator(6,7)
  assert_equal [1,2,3,0,1], grid.equator(4,5)
 end

 def testEquatorGapInMiddle
  grid = Grid.new(6,3,0,0)
  grid.addCell(4,5,1,2)
  grid.addCell(4,5,2,3)
  assert_equal [1,0], grid.gem(4,5)
  assert_equal [1,2,3,1], grid.equator(4,5)
 end

 def testEquatorGapInEnd
  assert_not_nil          grid = Grid.new(6,3,0,0)
  grid.addCell(4,5,1,2)
  grid.addCell(4,5,3,1)
  assert_equal [3,1,2,3], grid.equator(4,5)
 end

 def testEquatorTwoGaps
  assert_not_nil     grid = Grid.new(6,3,0,0)
  grid.addCell(4,5,1,2)
  grid.addCell(4,5,3,0)
  assert_nil         grid.equator(4,5)
 end

 def testAddNode
  assert_not_nil              grid = Grid.new(1,0,0,0)
  assert_nil                  grid.nodeXYZ(0)
  assert_equal 0,             grid.addNode(1.0,2.0,3.0)
  assert_equal [1.0,2.0,3.0], grid.nodeXYZ(0)
 end

 def testAddingNodeAllocatesMaxNode
  assert_not_nil     grid = Grid.new(1,0,0,0)
  assert_equal 1,    grid.maxnode
  assert_equal 0,    grid.addNode( 1.0, 2.0, 3.0)
  assert_equal 1,    grid.maxnode
  assert_equal 1,    grid.addNode(11.0,12.0,13.0)
  assert_equal 5001, grid.maxnode # chunk size hard coded to 5000
  assert_equal false,         grid.nodeFrozen(1)
  assert_equal [1,0,0,1,0,1], grid.map(1)
 end
 
 def testAddingNodeResizesAdj
  assert_not_nil     grid = Grid.new(1,0,0,0)
  assert_equal 0,    grid.addNode( 1.0, 2.0, 3.0)
  assert_equal 1,    grid.addNode(11.0,12.0,13.0)
  assert_equal 2,    grid.addNode(21.0,22.0,23.0)
  assert_equal 3,    grid.addNode(31.0,32.0,33.0)
  assert_equal 0,    grid.addCell(0,1,2,3)
  assert_equal 0,    grid.addFace(0,1,2,57)
  assert_equal 0,    grid.addEdge(0,1,82,40.0,50.0)
 end

 def testSetNodeXYZ
  assert_not_nil                 grid = Grid.new(1,0,0,0)
  assert_equal 0,                grid.addNode(1.0,2.0,3.0)
  assert_equal [1.0,2.0,3.0],    grid.nodeXYZ(0)
  assert_equal grid,             grid.setNodeXYZ(0,[10.0,20.0,30.0])
  assert_equal [10.0,20.0,30.0], grid.nodeXYZ(0)
 end

 def testGetAndSetNodeGlobal
  assert_not_nil      grid = Grid.new(1,0,0,0)
  assert_equal EMPTY, grid.nodeGlobal(-1)
  assert_equal EMPTY, grid.nodeGlobal(0)
  assert_equal EMPTY, grid.nodeGlobal(1)
  assert_nil          grid.setNodeGlobal(0,3)
  assert_equal 0,     grid.addNode(1.0,2.0,3.0)
  assert_equal grid,  grid.setNodeGlobal(0,5)
  assert_equal EMPTY, grid.nodeGlobal(-1)
  assert_equal 5,     grid.nodeGlobal(0)
  assert_equal EMPTY, grid.nodeGlobal(1)
 end

 def testFindNodeLocalFromGlobal
  assert_not_nil      grid = Grid.new(3,0,0,0)
  3.times { grid.addNode(1.0,2.0,3.0) }
  assert_equal EMPTY, grid.global2local(-1)
  assert_equal EMPTY, grid.global2local(2)
  3.times { |node| grid.setNodeGlobal(node,node+100) }
  assert_equal EMPTY, grid.global2local(-1)
  assert_equal EMPTY, grid.global2local(2)
  assert_equal 0, grid.global2local(100)
  assert_equal 1, grid.global2local(101)
  assert_equal 2, grid.global2local(102)
 end

 def testFindNodeLocalFromGlobalWithRemovedNode
  assert_not_nil      grid = Grid.new(3,0,0,0)
  3.times { grid.addNode(1.0,2.0,3.0) }
  3.times { |node| grid.setNodeGlobal(node,node+100) }
  grid.removeNode(1)
  assert_equal 0, grid.global2local(100)
  assert_equal 2, grid.global2local(102)
 end

 def testFindNodeLocalFromGlobalWithAddedNode
  assert_not_nil      grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,100)
  grid.addNodeWithGlobal(1.0,2.0,3.0,110)
  assert_equal 0, grid.global2local(100)
  node = grid.addNodeWithGlobal(1.0,2.0,3.0,254)
  assert_equal node, grid.global2local(254)
  node = grid.addNodeWithGlobal(1.0,2.0,3.0,8)
  assert_equal node, grid.global2local(8)
  node = grid.addNodeWithGlobal(1.0,2.0,3.0,105)
  assert_equal node, grid.global2local(105)
 end

 def testFindNodeLocalAfterNodeRemove
  assert_not_nil      grid = Grid.new(3,0,0,0)
  3.times { grid.addNode(1.0,2.0,3.0) }
  3.times { |node| grid.setNodeGlobal(node,node+100) }
  assert_equal 1,     grid.global2local(101)
  assert_equal 2,     grid.global2local(102)
  grid.removeNode(1)
  assert_equal EMPTY, grid.global2local(101)
  assert_equal 2,     grid.global2local(102)
 end

 def testFindNodeLocalAfterNodeRemoveAndPack
  assert_not_nil      grid = Grid.new(3,0,0,0)
  3.times { grid.addNode(1.0,2.0,3.0) }
  3.times { |node| grid.setNodeGlobal(node,node+100) }
  assert_equal 2, grid.global2local(102)
  grid.removeNode(1).pack
  assert_equal 1, grid.global2local(102)
 end

 def testGetAndSetNodePart
  assert_not_nil      grid = Grid.new(1,0,0,0)
  assert_equal EMPTY, grid.nodePart(-1)
  assert_equal EMPTY, grid.nodePart(0)
  assert_equal EMPTY, grid.nodePart(1)
  assert_nil          grid.setNodePart(0,3)
  assert_equal 0,     grid.addNode(1.0,2.0,3.0)
  assert_equal grid,  grid.setNodePart(0,5)
  assert_equal EMPTY, grid.nodePart(-1)
  assert_equal 5,     grid.nodePart(0)
  assert_equal EMPTY, grid.nodePart(1)
 end

 def testInitializeMeshWithOnlyLocalNodes
  assert_not_nil grid = Grid.new(1,0,0,0)
  assert_equal false, grid.nodeLocal(0)
  assert_equal false, grid.nodeGhost(0)
  grid.addNode(1,2,3)
  assert_equal true,  grid.nodeLocal(0)
  assert_equal false, grid.nodeGhost(0)
  assert_equal false, grid.nodeLocal(1)
  assert_equal false, grid.nodeGhost(1)
 end

 def testLocalAndGhostNodes
  assert_not_nil grid = Grid.new(2,0,0,0)
  grid.addNode(1,2,3)
  grid.addNode(1,2,3)
  grid.setNodePart(0,2)
  grid.setNodePart(1,5)
  grid.setPartId(2)
  assert_equal true,  grid.nodeLocal(0)
  assert_equal false, grid.nodeGhost(0)
  grid.setPartId(5)
  assert_equal false, grid.nodeLocal(0)
  assert_equal true,  grid.nodeGhost(0)
 end

 def testAddAndRemoveNode
  assert_nil          @grid.removeNode(5)
  assert_nil          @grid.removeNode(2)
  assert_nil          @grid.removeNode(-1)
  assert_equal 0,     @grid.addNode(1.0,2.0,3.0)
  assert_equal 1,     @grid.addNode(1.1,2.1,3.1)
  assert_equal 2,     @grid.addNode(1.2,2.2,3.2)
  assert_equal 3,     @grid.nnode
  assert_equal @grid, @grid.removeNode(1)
  assert_nil          @grid.removeNode(1)
  assert_equal 2,     @grid.nnode
  assert_not_nil      @grid.nodeXYZ(2)
  assert_nil          @grid.nodeXYZ(1)
  assert_equal 1,     @grid.addNode(1.2,2.2,3.2)
 end

 def testValidNode
  assert_not_nil      grid = Grid.new(4,0,0,0)
  assert_equal false, grid.validNode(-1)
  assert_equal false, grid.validNode(0)
  assert_equal false, grid.validNode(20)
  assert_equal 0,     grid.addNode(1.0,2.0,3.0)
  assert_equal true,  grid.validNode(0)
  assert_not_nil      grid.removeNode(0)
  assert_equal false, grid.validNode(0)
 end

 def testRemoveAndAddTracksGlobaNodeId
  grid = Grid.new(2,0,0,0)
  grid.addNode(1.0,2.0,3.0)
  grid.setGlobalNNode(1148)
  grid.setNodeGlobal(0,17)
  grid.removeNode(0)
  assert_equal [17], grid.getUnusedNodeGlobal
  grid.addNode(1.0,2.0,3.0)
  assert_equal 1148, grid.globalnnode
  assert_equal [], grid.getUnusedNodeGlobal
  grid.addNode(1.0,2.0,3.0)
  assert_equal 17,   grid.nodeGlobal(0)
  assert_equal 1148, grid.nodeGlobal(1)
  assert_equal 1149, grid.globalnnode
 end

 def testEliminateUnusedGlobalNodeId
  grid = Grid.new(5,0,0,0)
  grid.setGlobalNNode(100)
  grid.addNodeWithGlobal(1.0,2.0,3.0,99)
  4.times { grid.addNode(0,1,2) }
  grid.removeNode(1)
  grid.removeNode(3)
  assert_equal [100,102], grid.getUnusedNodeGlobal
  assert_equal 104, grid.globalnnode 
  assert_equal  99, grid.nodeGlobal(0)
  assert_equal 101, grid.nodeGlobal(2)
  assert_equal 103, grid.nodeGlobal(4)

  assert_equal grid, grid.eliminateUnusedNodeGlobal

  assert_equal [], grid.getUnusedNodeGlobal
  assert_equal 102, grid.globalnnode 
  assert_equal  99, grid.nodeGlobal(0)
  assert_equal 100, grid.nodeGlobal(2)
  assert_equal 101, grid.nodeGlobal(4)  
  assert_equal 0, grid.global2local(99)
  assert_equal 2, grid.global2local(100)
  assert_equal 4, grid.global2local(101)
 end

 def testRenumberGlobalNodesIdentity
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,1)
  assert_not_nil grid.renumberGlobalNodes([0,1])
  assert_equal 0, grid.nodeGlobal(0)
  assert_equal 1, grid.nodeGlobal(1)
 end

 def testRenumberGlobalNodesMoveDown
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,1)
  grid.addNodeWithGlobal(1.0,2.0,3.0,0)
  assert_not_nil grid.renumberGlobalNodes([1])
  assert_equal 0, grid.nodeGlobal(0)
  assert_equal 1, grid.nodeGlobal(1)
 end

 def testRenumberGlobalNodesMoveAll
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,1)
  grid.addNodeWithGlobal(1.0,2.0,3.0,0)
  assert_not_nil grid.renumberGlobalNodes([1,0])
  assert_equal 0, grid.nodeGlobal(0)
  assert_equal 1, grid.nodeGlobal(1)
 end

 def testRenumberGlobalNodesMoveWithoutSwitchForOffProc
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,0)
  assert_not_nil grid.renumberGlobalNodes([1])
  assert_equal 1, grid.nodeGlobal(0)
 end

 def testRenumberGlobalNodesMoveForOffProc
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,4)
  assert_not_nil grid.renumberGlobalNodes([4])
  assert_equal 0, grid.nodeGlobal(0)
 end

 def testRenumberGlobalNodesMoveBigSmall
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,10)
  grid.addNodeWithGlobal(1.0,2.0,3.0,0)
  assert_not_nil grid.renumberGlobalNodes([10,0])
  assert_equal 0, grid.nodeGlobal(0)
  assert_equal 1, grid.nodeGlobal(1)
 end

 def testRenumberGlobalNodesShift
  grid = Grid.new(5,0,0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,0)
  grid.addNodeWithGlobal(1.0,2.0,3.0,1)
  grid.addNodeWithGlobal(1.0,2.0,3.0,2)
  grid.addNodeWithGlobal(1.0,2.0,3.0,3)
  grid.addNodeWithGlobal(1.0,2.0,3.0,4)
  assert_not_nil grid.renumberGlobalNodes([0,3,1])
  assert_equal 0, grid.nodeGlobal(0)
  assert_equal 2, grid.nodeGlobal(1)
  assert_equal 3, grid.nodeGlobal(2)
  assert_equal 1, grid.nodeGlobal(3)
  assert_equal 4, grid.nodeGlobal(4)
 end

 def testNumberOfFaces
  assert_not_nil  grid = Grid.new(4,1,2,0)
  assert_equal 0, grid.nface 
  assert_equal 2, grid.maxface 
 end

 def testAddAndFindFace
  assert_not_nil           grid = Grid.new(4,1,2,0)
  assert_nil               grid.face(0)
  assert_equal 0,          grid.addFace(0, 1, 2, 10)
  assert_equal 0,          grid.findFace(0,1,2)
  assert_equal [0,1,2,10], grid.face(0)
  assert_nil               grid.findFace(3,1,2)
 end

 def testAddAndRemoveFace
  assert_not_nil     grid = Grid.new(4,1,2,0)
  assert_nil         grid.removeFace(0)
  assert_nil         grid.removeFace(1)
  assert_equal 0,    grid.addFace(0, 1, 2, 10)
  assert_nil         grid.removeFace(-1)
  assert_nil         grid.removeFace(1)
  assert_equal 1,    grid.addFace(3, 1, 2, 11)
  assert_equal 2,    grid.nface 
  assert_nil         grid.removeFace(3)
  assert_nil         grid.removeFace(-1)
  assert_equal 2,    grid.nface 
 end

 def testAddFaceReallocsFaceArray
  assert_not_nil     grid = Grid.new(4,0,1,0)
  assert_equal 0,    grid.addFace(0, 1, 2, 10)
  assert_equal 1,    grid.maxface
  assert_equal 1,    grid.addFace(3, 1, 2, 11)
  assert_equal 5001, grid.maxface
 end

 def testDeleteFacesWithThawedNodes
  assert_not_nil             grid = Grid.new(10,10,10,10)
  assert_equal 0,            grid.addNode(0,0,0)
  assert_equal 1,            grid.addNode(1,0,0)
  assert_equal 2,            grid.addNode(0,1,0)
  assert_equal 3,            grid.addNode(1,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(2,1,3,1)
  grid.addFace(2,1,3,2)
  assert_equal(-1,           grid.nThawedFaces(-1)) 
  assert_equal(-1,           grid.nThawedFaces(0))
  assert_equal 2,            grid.nThawedFaces(1)
  assert_equal 1,            grid.nThawedFaces(2)
  assert_equal grid,         grid.freezeNode(0)
  assert_equal grid,         grid.freezeNode(1)
  assert_equal grid,         grid.freezeNode(2)
  assert_equal 1,            grid.nThawedFaces(1)
  assert_nil                 grid.deleteThawedFaces(-1)
  assert_nil                 grid.deleteThawedFaces(0)
  assert_equal grid,         grid.deleteThawedFaces(1)
  assert_equal 0,            grid.nThawedFaces(1)
  assert_equal 2,            grid.nface
  assert_equal [0,1,2,1],    grid.face(0)
  assert_nil                 grid.face(1)
  assert_equal [2,1,3,2],    grid.face(2)
 end

 def testFaceId
  assert_not_nil     grid = Grid.new(4,1,2,0)

  assert_nil         grid.faceId( 1, 2, 3 )

  grid.addFace(0, 1, 2, 10)
  assert_equal 10,   grid.faceId( 0, 1, 2 )
  assert_equal 10,   grid.faceId( 1, 2, 0 )
  assert_equal 10,   grid.faceId( 2, 0, 1 )
  assert_equal 10,   grid.faceId( 2, 1, 0 )
  assert_nil         grid.faceId( 1, 2, 3 )

  grid.addFace(3, 1, 2, 11)
  assert_equal 10,   grid.faceId( 0, 1, 2 )
  assert_equal 11,   grid.faceId( 1, 2, 3 )
 end

 def testReconnectFace
  assert_not_nil             grid = Grid.new(6,0,3,0)
  grid.addFace(0,3,1,1)
  grid.addFace(0,1,2,2)
  grid.addFace(2,1,4,2)

  assert_equal [0, 3, 1, 1], grid.face(0)
  assert_equal [0, 1, 2, 2], grid.face(1)
  assert_equal [2, 1, 4, 2], grid.face(2)
  assert_nil                 grid.reconnectAllFace(-1,-1)
  assert_nil                 grid.reconnectAllFace(8, 8)
  assert_equal grid,         grid.reconnectAllFace(3, 3)
  assert_equal grid,         grid.reconnectAllFace(3, 5)
  assert_equal [0, 5, 1, 1], grid.face(0)
  assert_equal [0, 1, 2, 2], grid.face(1)
  assert_equal [2, 1, 4, 2], grid.face(2)
 end

 def testNodeUV
  assert_not_nil     grid = Grid.new(4,1,2,0)
  grid.addFaceUV(0,20.0,120.0,
				    1,21.0,121.0,
				    2,22.0,122.0,2)
  grid.addFaceUV(0,30.0,130.0,
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

 def testNodeFaceIdDegree
  grid = Grid.new(3,0,20,0)
  assert_equal EMPTY, grid.nodeFaceIdDegree(0)
  assert_nil          grid.nodeFaceId(0)
  3.times{grid.addNode(1,2,3)}
  assert_equal 0,  grid.nodeFaceIdDegree(0)
  assert_equal [], grid.nodeFaceId(0)
  grid.addFace(0,1,2,4)
  assert_equal 1,   grid.nodeFaceIdDegree(0)
  assert_equal [4], grid.nodeFaceId(0)
  grid.addFace(0,1,2,4)
  assert_equal 1, grid.nodeFaceIdDegree(0)
  assert_equal [4], grid.nodeFaceId(0)
  grid.addFace(0,1,2,5)
  assert_equal 2, grid.nodeFaceIdDegree(0)
  assert_equal [4,5], grid.nodeFaceId(0)
  grid.addFace(0,1,2,4)
  assert_equal 2, grid.nodeFaceIdDegree(0)
  assert_equal [4,5], grid.nodeFaceId(0)
  grid.addFace(0,1,2,2)
  assert_equal 3, grid.nodeFaceIdDegree(0)
  assert_equal [2,4,5], grid.nodeFaceId(0)
 end

 def testNodeFaceIdDegreeRandom
  grid = Grid.new(3,0,20,0)
  assert_equal EMPTY, grid.nodeFaceIdDegree(0)
  3.times{grid.addNode(1,2,3)}
  grid.addFace(0,1,2,7)
  grid.addFace(0,1,2,3)
  grid.addFace(0,1,2,5)
  grid.addFace(0,1,2,4)
  grid.addFace(0,1,2,5)
  grid.addFace(0,1,2,2)
  grid.addFace(0,1,2,8)
  grid.addFace(0,1,2,1)
  assert_equal [1,2,3,4,5,7,8], grid.nodeFaceId(0)
 end

 def testNumberOfGeomEdges
  assert_not_nil  grid = Grid.new(0,0,0,2)
  assert_equal 0, grid.nedge
  assert_equal 2, grid.maxedge
 end

 def testAddAndFindGeomEdge
  assert_not_nil         grid = Grid.new(4,0,0,2)
  assert_nil             grid.findEdge(0,1)
  assert_equal 0,        grid.addEdge(0, 1, 10, 0.0, 1.0)
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
  assert_equal 0,    grid.addEdge(0, 1, 10, 0.0, 1.0)
  assert_equal 0.0,  grid.nodeT(0,10)
  assert_equal 1.0,  grid.nodeT(1,10)
  assert_equal 1,    grid.addEdge(1, 2, 10, 1.0, 2.0)
  assert_equal grid, grid.setNodeT(1,10,1.5) 
  assert_equal 1.5,  grid.nodeT(1,10) 
 end

 def testAddAndRemoveGeomEdge
  assert_not_nil     grid = Grid.new(4,0,0,2)
  assert_nil         grid.removeEdge(-1)
  assert_nil         grid.removeEdge(0)
  assert_nil         grid.removeEdge(1)
  assert_equal 0,    grid.addEdge(0, 1, 10, 0.0, 0.0)
  assert_equal 10,   grid.edgeId(1, 0)
  assert_equal grid, grid.removeEdge(0)
  assert_nil         grid.edgeId(1, 0)
  assert_equal 0,    grid.addEdge(3, 1, 11, 0.0, 0.0)
  assert_equal 1,    grid.addEdge(0, 2, 12, 0.0, 0.0)
  assert_equal 11,   grid.edgeId(3, 1)
  assert_equal 2,    grid.nedge
 end

 def testNodeEdgeIdDegree
  grid = Grid.new(2,0,10,0)
  assert_equal EMPTY, grid.nodeEdgeIdDegree(0)
  assert_nil          grid.nodeEdgeId(0)
  2.times{grid.addNode(1,2,3)}
  assert_equal 0,  grid.nodeEdgeIdDegree(0)
  assert_equal [], grid.nodeEdgeId(0)
  grid.addEdge(0,1,4,0.0,1.0)
  assert_equal 1,   grid.nodeEdgeIdDegree(0)
  assert_equal [4], grid.nodeEdgeId(0)
  grid.addEdge(0,1,4,0.0,1.0)
  assert_equal 1, grid.nodeEdgeIdDegree(0)
  assert_equal [4], grid.nodeEdgeId(0)
  grid.addEdge(0,1,5,0.0,1.0)
  assert_equal 2, grid.nodeEdgeIdDegree(0)
  assert_equal [4,5], grid.nodeEdgeId(0)
  grid.addEdge(0,1,4,0.0,1.0)
  assert_equal 2, grid.nodeEdgeIdDegree(0)
  assert_equal [4,5], grid.nodeEdgeId(0)
  grid.addEdge(0,1,2,0.0,1.0)
  assert_equal 3, grid.nodeEdgeIdDegree(0)
  assert_equal [2,4,5], grid.nodeEdgeId(0)
 end

 def testNodeEdgeIdDegreeRandom
  grid = Grid.new(3,0,10,0)
  assert_equal EMPTY, grid.nodeEdgeIdDegree(0)
  2.times{grid.addNode(1,2,3)}
  grid.addEdge(0,1,7,0.0,1.0)
  grid.addEdge(0,1,3,0.0,1.0)
  grid.addEdge(0,1,5,0.0,1.0)
  grid.addEdge(0,1,4,0.0,1.0)
  grid.addEdge(0,1,5,0.0,1.0)
  grid.addEdge(0,1,2,0.0,1.0)
  grid.addEdge(0,1,8,0.0,1.0)
  grid.addEdge(0,1,2,0.0,1.0)
  grid.addEdge(0,1,1,0.0,1.0)
  assert_equal [1,2,3,4,5,7,8], grid.nodeEdgeId(0)
 end

 def testAddEdgeRealloc
  assert_not_nil     grid = Grid.new(4,0,0,1)
  assert_equal 0,    grid.addEdge(0, 1, 10, 0.0, 1.0)
  assert_equal 1,    grid.maxedge
  assert_equal 1,    grid.addEdge(1, 2, 11, 1.0, 2.0)
  assert_equal 5001, grid.maxedge
  assert_equal 11,   grid.edgeId(1, 2)
 end

 def testDeleteEdgesWithThawedNode
  assert_not_nil             grid = Grid.new(10,0,0,10)
  assert_equal 0,            grid.addNode(0,0,0)
  assert_equal 1,            grid.addNode(1,0,0)
  assert_equal 2,            grid.addNode(2,0,0)
  assert_equal 3,            grid.addNode(3,0,0)
  assert_equal 0,            grid.addEdge(0,1,1,0.0,1.0)
  assert_equal 1,            grid.addEdge(2,1,1,2.0,1.0)
  assert_equal 2,            grid.addEdge(2,3,1,2.0,3.0)
  assert_equal 3,            grid.addEdge(1,3,99,99.0,99.0)
  assert_equal 4,            grid.nedge
  assert_equal(-1,           grid.nThawedEdgeSegments(-1))
  assert_equal(-1,           grid.nThawedEdgeSegments(0))
  assert_equal 3,            grid.nThawedEdgeSegments(1)
  assert_equal 1,            grid.nThawedEdgeSegments(99)
  assert_equal grid,         grid.freezeNode(0)
  assert_equal grid,         grid.freezeNode(1)
  assert_equal 2,            grid.nThawedEdgeSegments(1)
  assert_nil                 grid.deleteThawedEdgeSegments(-1)
  assert_nil                 grid.deleteThawedEdgeSegments(0)
  assert_equal grid,         grid.deleteThawedEdgeSegments(1)
  assert_equal 2,            grid.nedge
  assert_equal 0,            grid.nThawedEdgeSegments(1)
  assert_equal [0,1,1],  grid.edge(0)
  assert_nil             grid.edge(1)
  assert_nil             grid.edge(2)
  assert_equal [1,3,99], grid.edge(3)
 end

 def testGetGeomCurve
  assert_not_nil     grid = Grid.new(4,0,0,4)
  assert_equal 0,    grid.addEdge(0, 1, 10, 10.0, 11.0)
  assert_equal 1,    grid.addEdge(1, 2, 11, 1.0, 2.0)
  assert_equal 2,    grid.addEdge(2, 3, 11, 2.0, 3.0)
  assert_equal 3,    grid.addEdge(3, 1, 12, 23.0, 21.0)
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

 def testGetGeomCurveLoop
  assert_not_nil     grid = Grid.new(4,0,0,4)
  assert_equal 0,    grid.addEdge(0, 1, 10, 1.0, 2.0)
  assert_equal 1,    grid.addEdge(1, 2, 10, 1.0, 2.0)
  assert_equal 2,    grid.addEdge(2, 0, 10, 2.0, 0.0)
  assert_equal 4,            grid.geomCurveSize(10,0)
  assert_equal [0, 2, 1, 0], grid.geomCurve(10,0)
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
  grid.addEdge(0,1,10,0.0,1.0)
  assert_equal false, grid.geometryEdge(-1)
  assert_equal true,  grid.geometryEdge(0)
  assert_equal false, grid.geometryEdge(grid.maxnode)
 end

 def testGeometryFace
  assert_not_nil grid = Grid.new(3,0,0,1)
  assert_equal false, grid.geometryFace(0)
  grid.addFace(0,1,2,10)
  assert_equal false, grid.geometryFace(-1)
  assert_equal true,  grid.geometryFace(0)
  assert_equal false, grid.geometryFace(grid.maxnode)
 end

 def testPack
  assert_not_nil grid = Grid.new(5,2,2,2)
  assert_equal 0, grid.addNode(9.0,0.0,9.0)
  assert_equal 1, grid.addNode(0.0,0.0,0.0)
  grid.addCell( 
	       0, 
	       grid.addNode(1.0,0.0,0.0), 
	       grid.addNode(0.0,1.0,0.0), 
	       grid.addNode(0.0,0.0,1.0) )
  assert_equal [0,2,3,4], grid.cell(0)
  grid.addCell(1,2,3,4)
  assert_equal [1,2,3,4], grid.cell(1)
  assert_equal [1], grid.gem(1,2)
  grid.addFace(2,3,4,1)
  grid.addFaceUV(1,1.0,11.0,
				    2,2.0,12.0,
				    3,3.0,13.0,
				    1)
  assert_equal [1.0,11.0], grid.nodeUV(1,1)
  assert_equal [2.0,12.0], grid.nodeUV(2,1)
  assert_equal [3.0,13.0], grid.nodeUV(3,1)
  grid.addEdge(2,3,1,22.0,23.0)
  grid.addEdge(1,2,1,21.0,22.0)
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
  grid.addCell(0,2,3,4)
  assert_equal [0,2,3,4], grid.cell(1)
  grid.addFace(1,2,3,2)
  assert_equal 1, grid.findFace(1,2,3)
  grid.addEdge(2,3,1,22.0,23.0)
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
  grid.addEdge(0,2,1,0.0,2.0)
  grid.addEdge(2,1,1,2.0,1.0)
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
  grid.addEdge(0,1,1,0.0,1.0)
  grid.addEdge(2,1,1,2.0,1.0)
  grid.addEdge(2,3,1,2.0,3.0)
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
  grid.addCell( 1, 4, 2, 3)
  grid.addCell( 0, 2, 4, 3)
  grid.addCell( 0, 4, 5, 3)
  assert_equal [0,4,5,3], grid.cell(2)
  grid.addFace( 1, 4, 2, 2)
  grid.addFace( 0, 2, 4, 2)
  grid.addFace( 0, 4, 5, 1)
  assert_equal 2,         grid.findFace(0, 4, 5)
  assert_equal 0,         grid.addEdge( 0, 4, 1, 0.0, 4.0)
  assert_equal 0,         grid.findEdge(0,4)
  assert_equal 1,         grid.addEdge( 1, 4, 1, 1.0, 4.0)
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

 def XtestTec
  assert_not_nil     grid = Grid.new(3,0,1,0)
  grid.addFace(grid.addNode(0,0,0),
	       grid.addNode(1,0,0),
	       grid.addNode(0,1,0),
	       1)
  assert_equal grid, grid.writeTecplotSurfaceZone
 end

 def testInsertPrism
  grid = Grid.new(6,0,0,0)
  assert_equal 0,             grid.nprism
  assert_equal grid,          grid.addPrism(0,1,2,3,4,5)
  assert_equal 1,             grid.nprism
  assert_nil                  grid.prism(-1) 
  assert_nil                  grid.prism(1)
  assert_equal [0,1,2,3,4,5], grid.prism(0)
 end

 def testInsertPyramid
  grid = Grid.new(5,0,0,0)
  assert_equal 0,           grid.npyramid
  assert_equal grid,        grid.addPyramid(0,1,2,3,4)
  assert_equal 1,           grid.npyramid
  assert_nil                grid.pyramid(-1) 
  assert_nil                grid.pyramid(1)
  assert_equal [0,1,2,3,4], grid.pyramid(0)
 end

 def testInsertQuad
  grid = Grid.new(4,0,0,0)
  assert_equal 0,            grid.nquad
  assert_equal grid,         grid.addQuad(0,1,2,3,33)
  assert_equal 1,            grid.nquad
  assert_nil                 grid.quad(-1) 
  assert_nil                 grid.quad(1)
  assert_equal [0,1,2,3,33], grid.quad(0)
 end

 def testMapMatrix
  defaultMap = [1,0,0,1,0,1]
  assert_nil      @grid.map(0)
  assert_equal 0, @grid.addNode(0,0,0)
  assert_equal    defaultMap, @grid.map(0)
  assert_nil      @grid.setMap(-1, 1,0,0, 1,0, 1)
  assert_nil      @grid.setMap(@grid.nnode, 1,0,0, 1,0, 1)
  assert_nil      @grid.map(-1)
  assert_nil      @grid.map(@grid.nnode)
  assert_equal    defaultMap, @grid.map(0)
  assert_equal    @grid, @grid.setMap(0, 1,2,3,4,5,6)
  assert_equal [1,2,3,4,5,6], @grid.map(0)
  assert_equal    @grid, @grid.removeNode(0)
  assert_nil      @grid.map(0)
  assert_equal 0, @grid.addNode(0,0,0)
  assert_equal    defaultMap, @grid.map(0) 
 end

 # make register unique

 # allocating a new chunk of nodes, faces, cells

end
