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
    @grid = Grid.new(4,1,0)
  end

  def testCreateGrid
    assert_equal 4, @grid.nnode
    assert_equal 1, @grid.maxcell
    assert_equal 0, @grid.ncell
  end

  def testNodeCellDeg
    assert_equal 0, @grid.nodeDeg(0)
    @grid.registerNodeCell( 0, 299 )
    assert_equal 1, @grid.nodeDeg(0)
  end

  def testNodeCellIterator
    assert_equal false, @grid.validNodeCell
    assert_equal false, @grid.moreNodeCell

    @grid.firstNodeCell(0);
    assert_equal false, @grid.validNodeCell

    @grid.registerNodeCell( 2, 299 )
    @grid.firstNodeCell(2);
    assert_equal 299, @grid.currentNodeCell
    assert_equal true,  @grid.validNodeCell
    assert_equal false, @grid.moreNodeCell
    @grid.nextNodeCell
    assert_equal false, @grid.validNodeCell

    @grid.registerNodeCell( 3, 398 )
    @grid.registerNodeCell( 3, 399 )
    @grid.firstNodeCell(3);
    assert_equal true,  @grid.validNodeCell
    assert_equal true,  @grid.moreNodeCell

    100.times {@grid.nextNodeCell} # abusive use of next
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

  def testEfficientStorage
    grid = Grid.new(1,1,1)
    assert_equal nil, grid.registerNodeCell(0,0)
    grid = Grid.new(1,1,2)
    assert_equal nil, grid.registerNodeCell(0,0)
    grid = Grid.new(1,1,3)
    assert_equal grid, grid.registerNodeCell(0,0)
    grid = Grid.new(1,1,4)
    assert_equal grid, grid.registerNodeCell(0,0)
    assert_equal grid, grid.registerNodeCell(0,1)
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

 def testAddCell
   assert_equal 0, @grid.ncell
   @grid.addCell(0,1,2,3)
   assert_equal 1, @grid.ncell
   (0..3).each { |n| assert_equal 1, @grid.nodeDeg(n)}
 end

 def testGetGem1
   grid = @grid
   grid.addCell(0,1,2,3)
   gem = grid.getGem(0,1)
   assert_equal [0], gem
 end

 def XtestGetGem2
   grid = Grid.new(5,2,20)
   grid.addCell(0,1,2,3)
   grid.addCell(0,1,2,4)
   grid.dump
   gem = grid.getGem(0,1)
   assert_equal [0, 1], gem
 end

# make register unique
# non-contiguos cellist for access and registering
# test that new list terminator is contiguous
# packing

# adding a cell bigger than maxcell

# get rid of long, use int

# allocating a new chunk of celllist

end
