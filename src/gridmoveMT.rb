#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('GridMove').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMove/GridMove'

class TestGridMove < Test::Unit::TestCase

 EMPTY = (-1)
 
 def equalTetWithFaceBase
  grid = Grid.new(4,1,4,0)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.866,0)
  grid.addNode(0.5,0.35,0.8)
  grid.addCell(0,1,2,3)
  grid.addFace(0,1,2,10)
  grid.addFace(0,3,1,10)
  grid.addFace(1,3,2,10)
  grid.addFace(0,2,3,10)
  grid
 end

 def testSpecifiedDisplacement
  grid = equalTetWithFaceBase
  assert_not_nil gm = GridMove.new(grid)
  zero = [0.0,0.0,0.0]
  up   = [0.0,0.0,1.0]
  assert_nil gm.displacement(-1)
  assert_nil gm.displacement(4)
  assert_equal false, gm.specified(-1)
  assert_equal false, gm.specified(4)
  assert_equal false, gm.specified(0)
  assert_equal zero,  gm.displacement(0)
  assert_equal gm,    gm.displace(0,up)
  assert_equal true,  gm.specified(0)
  assert_equal up,    gm.displacement(0)
 end

 def testInitGC
  grid = equalTetWithFaceBase
  gm = GridMove.new(grid)

  assert_not_nil  grid = String.new("hello")
  assert_nil      GC.start

  assert_nil gm.displacement(-1)
  assert_nil gm.displacement(4)
  zero = [0.0,0.0,0.0]
  assert_equal zero, gm.displacement(0)
 end

 def testNodeRealloc
  grid = Grid.new(1,0,0,0)
  grid.addNode(0.0,0.0,0.0)
  gm = GridMove.new(grid)
  1000.times{grid.addNode(0.0,0.0,0.0)}

  zero = [0.0,0.0,0.0]
  assert_equal false, gm.specified(1000)
  assert_equal zero,  gm.displacement(1000)
 end

 def testNodePack
  grid = Grid.new(4,0,0,0)
  4.times { grid.addNode(0.0,0.0,0.0) }
  gm = GridMove.new(grid)
  zero = [0.0,0.0,0.0]
  d = [1.1,2.2,3.3]
  assert_equal gm, gm.displace(3,d)
  grid.removeNode(0)
  grid.pack
  assert_equal false, gm.specified(3)
  assert_equal zero,  gm.displacement(3)
  assert_equal true,  gm.specified(2)
  assert_equal d,     gm.displacement(2)
 end

 def testSprings
  grid = Grid.new(5,2,0,0)
  gm = GridMove.new(grid)
  5.times { grid.addNode(0.0,0.0,0.0) }
  grid.addCell(0,1,2,3)
  assert_equal [0,1, 0,2, 0,3, 1,2, 1,3, 2,3], gm.springs  
  grid.addCell(1,0,2,4)
  assert_equal [0,1, 0,2, 0,3, 1,2, 1,3, 2,3, 1,4, 0,4, 2,4], gm.springs  
 end

 def testSpringRelaxationUp
  grid = equalTetWithFaceBase
  assert_not_nil gm = GridMove.new(grid)
  zero = [0.0,0.0,0.0]
  up  = [0.0,0.0,1.0]
  3.times{|n| gm.displace(n,up)}
  assert_equal zero, gm.displacement(3)
  assert_equal gm, gm.springRelaxation(1,1)
  assert_equal up, gm.displacement(0)
  assert_equal up, gm.displacement(1)
  assert_equal up, gm.displacement(2)
  delta = 1.0e-15
  assert_in_delta up[0], gm.displacement(3)[0], delta
  assert_in_delta up[1], gm.displacement(3)[1], delta
  assert_in_delta up[2], gm.displacement(3)[2], delta
 end

 def testSpringRelaxationUpSteps
  grid = equalTetWithFaceBase
  assert_not_nil gm = GridMove.new(grid)
  zero = [0.0,0.0,0.0]
  up  = [0.0,0.0,1.0]
  3.times{|n| gm.displace(n,up)}
  assert_equal zero, gm.displacement(3)
  assert_equal gm, gm.springRelaxation(10,1)
  assert_equal up, gm.displacement(0)
  assert_equal up, gm.displacement(1)
  assert_equal up, gm.displacement(2)
  delta = 1.0e-15
  assert_in_delta up[0], gm.displacement(3)[0], delta
  assert_in_delta up[1], gm.displacement(3)[1], delta
  assert_in_delta up[2], gm.displacement(3)[2], delta
 end

 def testSpringRelaxationUpSubiters
  grid = equalTetWithFaceBase
  assert_not_nil gm = GridMove.new(grid)
  zero = [0.0,0.0,0.0]
  up  = [0.0,0.0,1.0]
  3.times{|n| gm.displace(n,up)}
  assert_equal zero, gm.displacement(3)
  assert_equal gm, gm.springRelaxation(1,10)
  assert_equal up, gm.displacement(0)
  assert_equal up, gm.displacement(1)
  assert_equal up, gm.displacement(2)
  delta = 1.0e-15
  assert_in_delta up[0], gm.displacement(3)[0], delta
  assert_in_delta up[1], gm.displacement(3)[1], delta
  assert_in_delta up[2], gm.displacement(3)[2], delta
 end

 def testSpringRelaxationRotate
  grid = equalTetWithFaceBase
  assert_not_nil gm = GridMove.new(grid)
  2.times{|n| gm.displace(n,[0.0,0.0,0.0])}
  gm.displace(2,[0.5,0,-0.866])
  gm.springRelaxation(1,1)

  puts gm.displacement(3)
  ans = [0.5,0.8,-0.35]
  delta = 1.0e-15
  assert_in_delta ans[0], gm.displacement(3)[0], delta
  assert_in_delta ans[1], gm.displacement(3)[1], delta
  assert_in_delta ans[2], gm.displacement(3)[2], delta
 end

end
