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
    assert_equal 1, @grid.ncell
  end

  def testNodeDeg
    assert_equal 0, @grid.nodeDeg(0)
    @grid.registerNodeCell( 0, 299 )
    assert_equal 1, @grid.nodeDeg(0)
  end

  def testCellIterator
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

  def testAddAndRemoveCell
    assert_equal false, @grid.cellExists(1,0)
    assert_equal nil,   @grid.removeNodeCell(1,0)
    assert_equal @grid, @grid.registerNodeCell(1,0)
    assert_equal true,  @grid.cellExists(1,0)
    assert_equal @grid, @grid.removeNodeCell(1,0)
    assert_equal false, @grid.cellExists(1,0)
    assert_equal nil,   @grid.removeNodeCell(1,0)
  end

  def testEfficientStorage
    localGrid = Grid.new(1,1,1)
    assert_equal nil, localGrid.registerNodeCell(0,0)
    localGrid = Grid.new(1,1,2)
    assert_equal nil, localGrid.registerNodeCell(0,0)
    localGrid = Grid.new(1,1,3)
    assert_equal localGrid, localGrid.registerNodeCell(0,0)
    localGrid = Grid.new(1,1,4)
    assert_equal localGrid, localGrid.registerNodeCell(0,0)
    assert_equal localGrid, localGrid.registerNodeCell(0,1)
  end

 def testMultipleCellExists
   assert_equal false, @grid.cellExists(1,198)
   @grid.registerNodeCell(1,198)
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

# make register unique
# non-contiguos cellist for access and registering
# test that new list terminator is contiguous
# packing

# allocating a new chunk of celllist

end
