#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb Layer grid.h master_header.h'

require 'test/unit'
require 'Grid/Grid'
require 'Layer/Layer'

class TestLayer < Test::Unit::TestCase

 def testInit
  assert_not_nil  grid = Grid.new(2,0,0,0)
  assert_not_nil  layer = Layer.new(grid)
  assert_equal 2, layer.maxnode
  assert_equal 0, layer.nfront
 end

 def testInitGC
  assert_not_nil  grid = Grid.new(2,0,0,0)
  assert_not_nil  layer = Layer.new(grid)
  assert_equal 2, layer.maxnode
  assert_not_nil  grid = String.new("hello")
  assert_nil      GC.start
  assert_equal 2, layer.maxnode
 end

 def testMakeFront
  assert_not_nil        grid = Grid.new(4,0,2,0)
  assert_equal grid,    grid.addFace(0,1,2,1)
  assert_equal grid,    grid.addFace(0,1,3,2)
  assert_equal 2,       grid.nface
  assert_not_nil        layer = Layer.new(grid)
  assert_equal 0,       layer.nfront
  assert_nil            layer.front(0)
  assert_equal layer,   layer.makeFront([1,2])
  assert_equal 2,       layer.nfront
  assert_equal [0,1,2], layer.front(0)
  assert_equal [0,1,3], layer.front(1)
 end

end
