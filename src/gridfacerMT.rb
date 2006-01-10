#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for gridfacer c lib

rebuild = false || true
if rebuild
 Dir.chdir ENV['srcdir'] if ENV['srcdir']
 require 'RubyExtensionBuilder'
 RubyExtensionBuilder.new('GridFacer').build
end

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMath/GridMath'
require 'GridMetric/GridMetric'
require 'GridCAD/GridCAD'
require 'GridInsert/GridInsert'
require 'GridFacer/GridFacer'

class Grid
 include GridMetric
 include GridCAD
 include GridInsert
end

class TestGridFacer < Test::Unit::TestCase

 EMPTY = (-1)
 
 def test_initialize
  grid = Grid.new(0,0,0,0)
  faceId = 5
  assert_not_nil gf = GridFacer.new(grid,faceId)
  assert_equal( faceId, gf.faceId )
  assert_equal( 0, gf.edges )
 end

 def test_edges_from_different_faces
  grid = Grid.new(4,0,2,0)
  4.times{ grid.addNode(0,0,0) }
  grid.addFace(0,1,2,1)
  grid.addFace(2,1,3,2)
  faceId = 1
  assert_not_nil gf = GridFacer.new(grid,faceId)
  assert_equal( 3, gf.edges )
 end

 def test_edges_from_same_faces
  grid = Grid.new(4,0,2,0)
  4.times{ grid.addNode(0,0,0) }
  grid.addFace(0,1,2,1)
  grid.addFace(2,1,3,1)
  faceId = 1
  assert_not_nil gf = GridFacer.new(grid,faceId)
  assert_equal( 5, gf.edges )
 end

end
