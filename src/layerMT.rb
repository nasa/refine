#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridMetric adj.c grid.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb Layer adj.c grid.h master_header.h'

require 'test/unit'
require 'Grid/Grid'
require 'GridMetric/GridMetric'
require 'Layer/Layer'

class Grid
 include GridMetric
end

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
  assert_nil            layer.frontNormals(0)
  assert_equal 0,       layer.normalRoot(0)
  assert_equal layer,   layer.makeNormal
  assert_equal 4,       layer.nnormal
  assert_equal [0,1,2], layer.frontNormals(0)
  assert_equal [0,1,3], layer.frontNormals(1)
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
  assert_equal 0,       layer.normalDeg(0)
  assert_nil            layer.normalFronts(0)
  assert_equal layer,   layer.makeFront([1,2,3])
  assert_equal 3,       layer.nfront
  assert_equal 0,       layer.normalDeg(0)
  assert_nil            layer.normalFronts(0)
  assert_equal layer,   layer.makeNormal
  assert_equal 4,       layer.nnormal
  assert_equal [0,1,2], layer.frontNormals(0)
  assert_equal [0,1,3], layer.frontNormals(1)
  assert_equal [1,2,3], layer.frontNormals(2)
  assert_equal 2,       layer.normalDeg(0)
  assert_equal 3,       layer.normalDeg(1)
  assert_equal [1,0],   layer.normalFronts(0)
  assert_equal [2,1,0], layer.normalFronts(1)
  assert_equal [2,0],   layer.normalFronts(2)
  assert_equal [2,1],   layer.normalFronts(3)
 end

 def testNormalDirection1Face
  assert_not_nil          grid = Grid.new(3,0,1,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal grid,      grid.addFace(0,1,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_nil              layer.frontDirection(0)
  assert_nil              layer.normalDirection(0)
  assert_equal layer,     layer.makeFront([1])
  assert_equal 1,         layer.nfront
  assert_nil              layer.normalDirection(0)
  assert_equal layer,     layer.makeNormal
  assert_equal 3,         layer.nnormal
  direction = [0.0,0.0,1.0]
  assert_equal direction, layer.frontDirection(0)
  assert_equal direction, layer.normalDirection(0)
  assert_equal direction, layer.normalDirection(1)
  assert_equal direction, layer.normalDirection(2)
 end

 def testNormalDirection2Face
  assert_not_nil          grid = Grid.new(4,0,2,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  assert_equal grid,      grid.addFace(0,1,2,1)
  assert_equal grid,      grid.addFace(0,3,1,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.makeFront([1])
  assert_equal 2,         layer.nfront
  assert_nil              layer.normalDirection(0)
  assert_equal layer,     layer.makeNormal
  assert_equal 4,         layer.nnormal
  assert_equal [0.0,0.0,1.0], layer.frontDirection(0)
  assert_equal [0.0,1.0,0.0], layer.frontDirection(1)
  halfSqrt2 = 0.5 * Math::sqrt(2)
  direction = [0.0, halfSqrt2, halfSqrt2]
  assert_equal direction, layer.normalDirection(0)
  assert_equal direction, layer.normalDirection(1)
  assert_equal [0.0,0.0,1.0], layer.normalDirection(2)
  assert_equal [0.0,1.0,0.0], layer.normalDirection(3)
 end

 def testAdvanceLayer
  assert_not_nil          grid = Grid.new(7,4,1,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  assert_equal grid,      grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  assert_equal grid,      grid.addFace(0,1,2,1)
  assert_equal true,      grid.rightHandedFace(0)
  assert_equal true,      grid.rightHandedBoundary
  assert       0<         grid.minVolume, "negative volumes"
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.makeFront([1])
  assert_equal 1,         layer.nfront
  assert_equal layer,     layer.makeNormal
  assert_equal 3,         layer.nnormal
  assert_equal layer,     layer.advance(0.1)
  assert_equal 7,         grid.nnode
  assert_equal [4,5,6,3], grid.cell(0)
  assert_equal 4,         grid.ncell
  assert_equal true,      grid.rightHandedFace(0)
  assert_equal true,      grid.rightHandedBoundary
  assert       0<         grid.minVolume, "negative volumes"
 end

end
