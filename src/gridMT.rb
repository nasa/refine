#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit unless system 'ruby extconf.rb Grid'
exit unless system 'make'

require 'test/unit'
require 'Grid'

class TestSampleUnit < Test::Unit::TestCase

 def set_up
  @grid = Grid.new(4,1)
 end

 def testCreateGrid
  assert_equal 4, @grid.maxnode
  assert_equal 0, @grid.nnode
  assert_equal 1, @grid.maxcell
  assert_equal 0, @grid.ncell
 end

 def testNodeCellDeg
  assert_equal 0, @grid.nodeDeg(0)
  assert_not_nil @grid.registerNodeCell( 0, 299 )
  assert_equal 1, @grid.nodeDeg(0)
 end

 def testNodeCellIterator
  assert_equal false, @grid.validNodeCell
  assert_equal false, @grid.moreNodeCell
  
  @grid.firstNodeCell(0);
  assert_equal false, @grid.validNodeCell
  
  assert_not_nil @grid.registerNodeCell( 2, 299 )
  @grid.firstNodeCell(2);
  assert_equal 299, @grid.currentNodeCell
  assert_equal true,  @grid.validNodeCell
  assert_equal false, @grid.moreNodeCell
  @grid.nextNodeCell
  assert_equal false, @grid.validNodeCell
  
  assert_not_nil @grid.registerNodeCell( 3, 398 )
  assert_not_nil @grid.registerNodeCell( 3, 399 )
  @grid.firstNodeCell(3);
    assert_equal true,  @grid.validNodeCell
  assert_equal true,  @grid.moreNodeCell
  
  100.times {@grid.nextNodeCell} # abusive use of next
  assert_nil @grid.firstNodeCell(300000);
end
 
 def testAddAndRemoveNodeCell
  assert_equal false, @grid.cellExists(1,0)
  assert_equal nil,   @grid.removeNodeCell(1,0)
  assert_equal @grid, @grid.registerNodeCell(1,0)
  assert_equal true,  @grid.cellExists(1,0)
  assert_equal @grid, @grid.removeNodeCell(1,0)
  assert_equal false, @grid.cellExists(1,0)
  assert_equal nil,   @grid.removeNodeCell(1,0)
 end
 
 def testRegisterNodeCellLimit
  assert_equal nil, @grid.registerNodeCell(1000,0)
 end

 def testMultipleNodeCellExists
  assert_equal false, @grid.cellExists(1,198)
  @grid.registerNodeCell(1,198)
  @grid.registerNodeCell(2,198)
  @grid.registerNodeCell(1,199)
  
  assert_equal true,  @grid.cellExists(1,198)
  assert_equal true,  @grid.cellExists(1,199)
  @grid.removeNodeCell(1,198)
  assert_equal false, @grid.cellExists(1,198)
  assert_equal true,  @grid.cellExists(1,199)
  @grid.registerNodeCell(1,198)
  assert_equal true,  @grid.cellExists(1,198)
  assert_equal true,  @grid.cellExists(1,199)
 end
 
 def testAddCellDegAndCell
  assert_equal 0, @grid.ncell
  assert_equal nil, @grid.cell(0)
  assert_equal nil, @grid.cell(5)
  assert_equal @grid, @grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], @grid.cell(0)
  assert_equal 1, @grid.ncell
  (0..3).each { |n| assert_equal 1, @grid.nodeDeg(n)}
 end

 def testAddCellRegFailure
  grid = Grid.new(3,1)
  assert_equal nil, grid.addCell(0,1,2,3)
  grid = Grid.new(4,1)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal nil, grid.addCell(0,1,2,3)
 end
 
 def testRemoveCell
  assert_equal @grid, @grid.addCell(0,1,2,3)
  assert_equal nil, @grid.removeCell(25625)
  assert_equal 1, @grid.ncell
  assert_equal @grid, @grid.removeCell(0)
  assert_equal nil, @grid.cell(0)
  assert_equal nil, @grid.removeCell(0)
  assert_equal 0, @grid.ncell
  (0..3).each { |n| assert_equal 0, @grid.nodeDeg(n)}
 end

 def testReplaceCell
  grid = Grid.new(8,2)
  assert_equal grid, grid.addCell(0,1,2,3).addCell(4,5,6,7)
  assert_equal grid, grid.removeCell(0)
  assert_equal grid, grid.addCell(0,1,2,3)
  assert_equal [0, 1, 2, 3], grid.cell(0)
  assert_equal [4, 5, 6, 7], grid.cell(1)
 end

 def testGetGem
  grid = Grid.new(5,3)
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
  grid = Grid.new(6,4)
  assert_equal grid, grid.
   addCell(4,5,1,2).addCell(4,5,2,3).addCell(4,5,3,0).addCell(4,5,0,1)
  assert_equal [3,2,1,0], grid.gem(4,5)
  assert_equal 2, grid.nodeDeg(0)
  assert_equal 4, grid.nodeDeg(5)
  assert_equal [], grid.equator(0,2)
  assert_equal [], grid.equator(6,7)
  assert_equal [1,2,3,0,1], grid.equator(4,5)
 end

 def testEquatorGapInMiddle
  grid = Grid.new(6,3)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,2,3)
  assert_equal [1,0], grid.gem(4,5)
  assert_equal [1,2,3,1], grid.equator(4,5)
 end

 def testEquatorGapInEnd
  grid = Grid.new(6,3)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,3,1)
  assert_equal [3,1,2,3], grid.equator(4,5)
 end

 def testEquatorTwoGaps
  grid = Grid.new(6,3)
  assert_equal grid, grid.addCell(4,5,1,2).addCell(4,5,3,0)
  assert_equal nil, grid.equator(4,5)
 end

 def testSwap4for4
  grid = Grid.new(6,4)
  assert_equal grid, grid.
   addCell(4,5,0,1).addCell(4,5,1,2).addCell(4,5,2,3).addCell(4,5,3,0)
  assert_equal grid, grid.swap(4,5)
  assert_equal 2, grid.nodeDeg(4)
  assert_equal 2, grid.nodeDeg(5)
  assert_equal 4, grid.nodeDeg(0)
  assert_equal 4, grid.nodeDeg(2)
  assert_equal 2, grid.nodeDeg(1)
  assert_equal 2, grid.nodeDeg(3)
 end

 def testAddNode
  grid = Grid.new(1,1)
  assert_equal 0, grid.addNode(1.0,0.0,0.0)
  assert_equal( -1, grid.addNode(1.0,0.0,0.0))
 end

 def testMetrics
  assert_equal @grid, @grid.
   addCell( @grid.addNode(0.0,0.0,0.0), @grid.addNode(1.0,0.0,0.0), 
	    @grid.addNode(0.0,1.0,0.0), @grid.addNode(0.0,0.0,1.0) )
  nodes = [0,1,2,3]
  assert_equal( 1.0/6.0, @grid.volume(nodes) )
  assert_in_delta 1.366025403784439, @grid.ar(nodes), 1.0e-15
 end

 def testGemGrid
  assert_not_nil grid=gemGrid(nequ=4, a=0.1)
  assert nequ+2, grid.nnode
  assert nequ, grid.ncell
  grid.ncell.times do |cellId|
   puts cellId, grid.volume(grid.cell(cellId)), grid.ar(grid.cell(cellId))
  end
  grid.swap(0,1)
  grid.ncell.times do |cellId|
   puts cellId, grid.volume(grid.cell(cellId)), grid.ar(grid.cell(cellId))
  end
 end

 def gemGrid(nequ=4, a=0.1)
  grid = Grid.new(nequ+2,nequ)
  n = Array.new
  n.push grid.addNode(1.0,0.0,0.0)
  n.push grid.addNode(-1.0,0.0,0.0)
  nequ.times do |i| 
   angle = 2.0*Math::PI*(i-1)/(nequ)
   n.push grid.addNode(0.0,a*Math.sin(angle),a*Math.cos(angle)) 
  end
  n.push 2
  nequ.times do |i|
   grid.addCell(n[0],n[1],n[i+2],n[i+3])
  end
  grid  
 end


 def XtestMaxSize
  nnode = 6000000
  grid = Grid.new(nnode,nnode*6)
  1.upto(nnode*6) {grid.addCell(3,4,0,1)}
 end


 # make register unique

 # allocating a new chunk of celllist

end
