#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

system 'ruby extconf.rb Grid'
system 'make'

require 'test/unit'
require 'Grid'

class TestSampleUnit < Test::Unit::TestCase

  def testCreateGrid
    @grid = Grid.new(4,1,0)
    assert_equal 4, @grid.nnode
    assert_equal 1, @grid.ncell
  end

end
