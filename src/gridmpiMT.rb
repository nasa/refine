#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('GridMPI').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Queue/Queue'
require 'Sort/Sort'
require 'Grid/Grid'
require 'GridMetric/GridMetric'
require 'GridInsert/GridInsert'
require 'GridMPI/GridMPI'

class Grid
 include GridMetric
 include GridInsert
 include GridMPI
end

class TestGridMPI < Test::Unit::TestCase

 EMPTY = -1

 def rightTet
  grid = Grid.new(5,2,4,0)
  grid.addCell( 
	       grid.addNode(0,0,0), 
	       grid.addNode(1,0,0), 
	       grid.addNode(0,1,0), 
	       grid.addNode(0,0,1) )
  grid.addFaceUV(0,0,0,
		 3,0,1,
		 1,1,0,
		 10)
  grid.addFaceUV(0,0,0,
		 2,1,0,
		 3,0,1,
		 11)
  grid.identityNodeGlobal(100).identityCellGlobal(200)
  grid.setGlobalNNode(104).setGlobalNCell(201)
  grid
 end 

 def testCopyLocalNodeNumberingToGlobalNumbering
  plus = 544
  grid = Grid.new(10,0,0,0)
  10.times { grid.addNode(1,2,3) }
  assert_equal grid,grid.identityNodeGlobal(plus)
  10.times { |node| assert_equal node+plus, grid.nodeGlobal(node) }
 end

 def testCopyLocalCellNumberingToGlobalNumbering
  plus = 654
  grid = Grid.new(4,10,0,0)
  4.times { grid.addNode(1,2,3) }
  10.times { grid.addCell(0,1,2,3) }
  assert_equal grid, grid.identityCellGlobal(plus)
  10.times { |cell| assert_equal cell+plus, grid.cellGlobal(cell) }
 end

 def testMarkAllNodesAsLocal
  partId = 5
  grid = Grid.new(10,0,0,0).setPartId(partId)
  10.times { grid.addNode(1,2,3) }
  assert_equal grid, grid.setAllLocal
  10.times do |node| 
   assert_equal partId, grid.nodePart(node), "node #{node} part failed" 
  end
 end

 def testSetGhostNode
  partId = 7
  grid = Grid.new(1,0,0,0).setPartId(partId)
  grid.addNode(1,2,3)
  grid.setAllLocal
  assert_equal grid, grid.setGhost(0)
  assert_equal EMPTY, grid.nodePart(0)
 end

 def testLoadQueue
  q = Queue.new
  p1 = rightTet.setPartId(1).setAllLocal
  p1.setGhost(3)

  assert_equal 4, p1.parallelEdgeSplit(q,0,3)
  assert_equal 2, p1.ncell
  assert_equal 4, p1.nface
  assert_equal 105, p1.globalnnode
  assert_equal 202, p1.globalncell
  assert_equal p1.partId, p1.nodePart(4)
  assert_equal 104, p1.nodeGlobal(4)
  assert_equal 200, p1.cellGlobal(0)
  assert_equal 201, p1.cellGlobal(1)
  assert_equal 2, q.transactions
  assert_equal 1, q.removedCells(1)
  assert_equal [100,101,102,103], q.removedCellNodes(0)
  assert_equal 2, q.addedCells(1)
  assert_equal [104,101,102,103,200], q.addedCellNodes(0)
  assert_equal [100,101,102,104,201], q.addedCellNodes(1)
  h=0.5
  assert_equal [ 0,0,h, 1,0,0, 0,1,0, 0,0,1 ], q.addedCellXYZs(0)
  assert_equal [ 0,0,0, 1,0,0, 0,1,0, 0,0,h ], q.addedCellXYZs(1)
  assert_equal 2, q.removedFaces(1)
  assert_equal 4, q.addedFaces(1)
 end

 def testUnPackQueue
  q = Queue.new
  p1 = rightTet.setPartId(1).setAllLocal
  p2 = rightTet.setPartId(2).setAllLocal
  p2.setGhost(0).setGhost(1).setGhost(2)
  p1.setGhost(3)

  assert_equal EMPTY, p2.parallelEdgeSplit(q,0,1), "split a ghost edge"
  assert_equal 1, p1.ncell
  assert_equal 1, q.transactions

  assert_equal 4, p1.parallelEdgeSplit(q,0,3)

  assert_equal p2, p2.applyQueue(q)
  assert_equal 0, p2.ncell
 end

end
