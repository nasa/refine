#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridCAD adj.c grid.c gridStruct.h'

require 'test/unit'
require 'Grid/Grid'
require 'GridCAD/GridCAD'

class Grid
 include GridCAD
 def totalVolume
  vol = 0.0
  ncell.times { |cellId| vol += volume(cell(cellId)) }
  vol
 end
end

class TestSampleUnit < Test::Unit::TestCase

 def testProjection
  assert_not_nil grid = Grid.new(0,0,0,0)
  assert_nil grid.safeProject(0)
 end

end
