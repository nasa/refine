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
require 'Grid/Grid'
require 'GridMPI/GridMPI'

class Grid
 include GridMPI
end

class TestGridMPI < Test::Unit::TestCase

 EMPTY = -1

 def rightTet
  grid = Grid.new(5,2,0,0)
  grid.addCell( 
	       grid.addNode(0,0,0), 
	       grid.addNode(1,0,0), 
	       grid.addNode(0,1,0), 
	       grid.addNode(0,0,1) )
  grid.identityGlobal
  grid
 end 

 def testCopyLocalNodeNumberingToGlobalNumbering
  grid = Grid.new(10,0,0,0)
  10.times { grid.addNode(1,2,3) }
  assert_equal grid,grid.identityGlobal
  10.times { |node| assert_equal node, grid.nodeGlobal(node) }
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
  p1 = rightTet.setPartId(1).setAllLocal
  p2 = rightTet.setPartId(2).setAllLocal
  p2.setGhost(0)
  p2.setGhost(1)
  p2.setGhost(2)
  p1.setGhost(3)
 end

end
