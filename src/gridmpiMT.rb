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
  grid.addFace(0,3,1,10)
  grid.addFace(0,2,3,11)
  grid.identityGlobal(100)
  grid
 end 

 def testCopyLocalNodeNumberingToGlobalNumbering
  plus = 544
  grid = Grid.new(10,0,0,0)
  10.times { grid.addNode(1,2,3) }
  assert_equal grid,grid.identityGlobal(544)
  10.times { |node| assert_equal node+544, grid.nodeGlobal(node) }
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

 def testSplitEdgeAcrossProc
  q = Queue.new
  p1 = rightTet.setPartId(1).setAllLocal
  p2 = rightTet.setPartId(2).setAllLocal
  p2.setGhost(0)
  p2.setGhost(1)
  p2.setGhost(2)
  p1.setGhost(3)

  assert_equal EMPTY, p2.parallelEdgeSplit(q,0,1), "split a ghost edge"
  assert_equal 1, p1.ncell
  assert_equal 1, q.transactions

  assert_equal 4, p1.parallelEdgeSplit(q,0,3)
  assert_equal 2, p1.ncell
  assert_equal 4, p1.nface
  assert_equal p1.partId, p1.nodePart(4)
  assert_equal 2, q.transactions
  assert_equal 1, q.removedCells(1)
  assert_equal [0,1,2,3], q.removedCellNodes(0)
 end

end
