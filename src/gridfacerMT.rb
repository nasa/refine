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
  assert_not_nil ge = GridFacer.new(grid,faceId)
  assert_equal( faceId, ge.faceId )
 end

end
