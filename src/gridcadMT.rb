#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridCAD FAKEGeom adj.c grid.c gridStruct.h master_header.h'

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

 def testEdgeProjection
  assert_not_nil              grid = Grid.new(3,0,0,2)
  assert_equal 0,             grid.addNode(0.0,0.0,0.0)
  assert_equal 1,             grid.addNode(0.5,0.1,0.1)
  assert_nil                  grid.projectNodeToEdge(1,1)
  assert_equal grid,          grid.addEdge(0,1,1,0.0,0.55)
  assert_equal grid,          grid.projectNodeToEdge(1,1)
  assert_equal [0.5,0.0,0.0], grid.nodeXYZ(1)
  assert_equal 0.0,           grid.nodeT(0,1)
  assert_equal 0.5,           grid.nodeT(1,1)
 end


 def testSafeProjection
  assert_not_nil grid = Grid.new(3,0,0,2)
  assert_equal 0, grid.addNode(0.0,0.0,0.0)
  assert_equal 1, grid.addNode(0.5,0.0,0.0)
  assert_equal 2, grid.addNode(1.0,0.0,0.0)
  assert_nil grid.safeProject(0)
 end

end
