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
    0.upto(@grid.nnode-1) { |i| assert_equal 0, @grid.nodeDeg(i) }
    @grid.registerNodeCell( 2, 299 )
    assert_equal 1, @grid.nodeDeg(2)
  end

  def testCellIterator
    assert_equal false, @grid.validNodeCell
    assert_equal false, @grid.moreNodeCell

    @grid.firstNodeCell(0);
    assert_equal false, @grid.validNodeCell

    @grid.registerNodeCell( 2, 299 )
    @grid.firstNodeCell(2);
    assert_equal 299, @grid.currentNodeCell
    assert_equal true, @grid.validNodeCell
    assert_equal false, @grid.moreNodeCell
    @grid.nextNodeCell
    assert_equal false, @grid.validNodeCell

    @grid.registerNodeCell( 3, 398 )
    @grid.registerNodeCell( 3, 399 )
    @grid.firstNodeCell(3);
    assert_equal true, @grid.validNodeCell
    assert_equal true, @grid.moreNodeCell

    100.times {@grid.nextNodeCell} # abusive use of next
  end

end
