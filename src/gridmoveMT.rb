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
  grid = Grid.new(4,1,1,0)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)
  grid.addNode(0.5,0.35,0.8)
  grid.addCell(0,1,2,3)
  grid.addFace(0,1,2,10)
  grid
 end

 def testDisplacement
  grid = equalTetWithFaceBase
  assert_not_nil gm = GridMove.new(grid)
  zero = [0.0,0.0,0.0]
  up   = [0.0,0.0,1.0]
  assert_nil gm.displacement(-1)
  assert_nil gm.displacement(4)
  assert_equal zero, gm.displacement(0)
  assert_equal gm, gm.displace(0,up)
  assert_equal up, gm.displacement(0)
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
  grid = equalTetWithFaceBase
  gm = GridMove.new(grid)
  1000.times{grid.addNode(0.0,0.0,0.0)}

  zero = [0.0,0.0,0.0]
  assert_equal zero, gm.displacement(1000)
 end

end
