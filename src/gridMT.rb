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

    @grid.firstNodeCell(0);
    assert_equal false, @grid.validNodeCell
  end

end
