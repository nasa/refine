#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for gridedger c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('GridEdger').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMath/GridMath'
require 'GridMetric/GridMetric'
require 'GridEdger/GridEdger'

class Grid
 include GridMetric
end

class TestGridEdger < Test::Unit::TestCase

 EMPTY = (-1)
 
 def test_initialize
  grid = Grid.new(0,0,0,0)
  edgeId = 5
  assert_not_nil ge = GridEdger.new(grid,edgeId)
  assert_equal( edgeId, ge.edgeId )
 end

end
