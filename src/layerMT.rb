#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb Layer adj.c grid.h master_header.h'

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

 def testMakeNormals
  assert_not_nil        grid = Grid.new(5,0,2,0)
  assert_equal grid,    grid.addFace(1,2,3,1)
  assert_equal grid,    grid.addFace(1,2,4,2)
  assert_not_nil        layer = Layer.new(grid)
  assert_nil            layer.makeNormal
  assert_equal layer,   layer.makeFront([1,2])
  assert_equal 0,       layer.nnormal
  assert_nil            layer.frontNormal(0)
  assert_equal 0,       layer.normalRoot(0)
  assert_equal layer,   layer.makeNormal
  assert_equal 4,       layer.nnormal
  assert_equal [0,1,2], layer.frontNormal(0)
  assert_equal [0,1,3], layer.frontNormal(1)
  assert_equal 1,       layer.normalRoot(0)
  assert_equal 2,       layer.normalRoot(1)
  assert_equal 3,       layer.normalRoot(2)
  assert_equal 4,       layer.normalRoot(3)
 end

 def testConstrainNormals
  assert_not_nil        grid = Grid.new(6,0,3,0)
  assert_equal grid,    grid.addFace(1,2,3,1)
  assert_equal grid,    grid.addFace(1,2,4,2)
  assert_equal grid,    grid.addFace(3,4,5,3)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.makeFront([1,2])
  assert_equal 2,       layer.nfront
  assert_nil            layer.constrainNormal(3)
  assert_equal 0,       layer.constrained(0)
  assert_equal layer,   layer.makeNormal
  assert_equal 4,       layer.nnormal
  assert_equal 0,       layer.constrained(0)
  assert_equal 0,       layer.constrained(1)
  assert_equal 0,       layer.constrained(2)
  assert_equal 0,       layer.constrained(3)
  assert_equal layer,   layer.constrainNormal(3)
  assert_equal 0,       layer.constrained(0)
  assert_equal 0,       layer.constrained(1)
  assert_equal 3,       layer.constrained(2)
  assert_equal 3,       layer.constrained(3)
 end

 def testNormalFrontNeighbors
  assert_not_nil        grid = Grid.new(4,0,3,0)
  assert_equal grid,    grid.addFace(0,1,2,1)
  assert_equal grid,    grid.addFace(0,1,3,2)
  assert_equal grid,    grid.addFace(1,2,3,3)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.makeFront([1,2,3])
  assert_equal 3,       layer.nfront
  assert_equal layer,   layer.makeNormal
  assert_equal 4,       layer.nnormal
  assert_equal [0,1,2], layer.frontNormal(0)
  assert_equal [0,1,3], layer.frontNormal(1)
  assert_equal [1,2,3], layer.frontNormal(2)
  assert_equal 2,       layer.normalDeg(0)
  assert_equal 3,       layer.normalDeg(1)
  assert_equal [1,0],   layer.normalFronts(0)
  assert_equal [2,1,0], layer.normalFronts(1)
  assert_equal [2,0],   layer.normalFronts(2)
  assert_equal [2,1],   layer.normalFronts(3)
 end

end
