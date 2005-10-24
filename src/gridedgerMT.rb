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
require 'GridCAD/GridCAD'
require 'GridEdger/GridEdger'

class Grid
 include GridMetric
 include GridCAD
end

class TestGridEdger < Test::Unit::TestCase

 EMPTY = (-1)
 
 def one_edge_grid
  grid = Grid.new(4,0,2,1)
  grid.addNode(0.0,0.0,0.0)
  grid.addNode(1.0,0.0,0.0)
  grid.addNode(0.5,0.5,0.0)
  grid.addNode(0.5,-0.5,0.0)

  edgeId = 1
  grid.addEdge(0,1,edgeId,0.0,1.0)

  #  2
  # 0 1
  #  3

  grid.addFaceUV(0, 10.0, 20.0,
                 1, 11.0, 20.0,
                 2, 10.5, 20.5,
                 1)

  grid.addFaceUV(0, 10.0, 20.0,
                 3, 10.5, 19.5,
                 1, 11.0, 20.0,
                 2)

  grid.setNGeomNode(2)
  grid.setNGeomEdge(1)
  grid.setNGeomFace(2)

  grid.addGeomEdge(edgeId,0,1)

  grid.scaleSpacing(0,0.5)
  grid.scaleSpacing(1,0.5)

  grid
 end

 def test_single_edge_setup
  grid = one_edge_grid

  assert grid.geometryNode(0)
  assert grid.geometryNode(1)
  assert grid.geometryFace(2)
  assert grid.geometryFace(3)

  assert_equal [0, 1], grid.geomEdge(1)

  assert_equal [0, 1],    grid.geomCurve(1,0)
  assert_equal [0.0, 1.0], grid.geomCurveT(1,0)
 end

 def two_edge_grid
  grid = Grid.new(5,0,4,2)
  grid.addNode(0.0,0.0,0.0)
  grid.addNode(3.0,0.0,0.0)
  grid.addNode(1.0,0.0,0.0)
  grid.addNode(0.5,0.5,0.0)
  grid.addNode(0.5,-0.5,0.0)

  edgeId = 1
  grid.addEdge(0,2,edgeId,0.0,1.0)
  grid.addEdge(2,1,edgeId,1.0,3.0)

  #  3
  # 0 2  1
  #  4

  grid.addFaceUV(0, 10.0, 20.0,
                 2, 11.0, 20.0,
                 3, 10.5, 20.5,
                 1)
  grid.addFaceUV(2, 11.0, 20.0,
                 1, 13.0, 20.0,
                 3, 10.5, 20.5,
                 1)

  grid.addFaceUV(0, 10.0, 20.0,
                 4, 10.5, 19.5,
                 2, 11.0, 20.0,
                 2)
  grid.addFaceUV(2, 11.0, 20.0,
                 4, 10.5, 19.5,
                 1, 13.0, 20.0,
                 2)

  grid.setNGeomNode(2)
  grid.setNGeomEdge(1)
  grid.setNGeomFace(2)

  grid.addGeomEdge(edgeId,0,1)

  grid.scaleSpacing(2,0.5)
  grid.scaleSpacing(1,0.25)

  grid
 end

 def test_two_edge_setup
  grid = two_edge_grid

  assert grid.geometryNode(0)
  assert grid.geometryNode(1)
  assert grid.geometryEdge(2)
  assert grid.geometryFace(3)
  assert grid.geometryFace(4)

  assert_equal [0, 2, 1], grid.geomEdge(1)

  assert_equal [0, 2, 1],    grid.geomCurve(1,0)
  assert_equal [1, 2, 0],    grid.geomCurve(1,1)
  assert_equal [0.0, 1.0, 3.0], grid.geomCurveT(1,0)
 end

 def test_initialize
  grid = Grid.new(0,0,0,0)
  edgeId = 5
  assert_not_nil ge = GridEdger.new(grid,edgeId)
  assert_equal( edgeId, ge.edgeId )
 end

 def test_get_t_from_segment_single_edge_out_of_range
  ge = GridEdger.new(one_edge_grid,1)
  assert_nil ge.segmentT(-1.0)
  assert_nil ge.segmentT(1.1)
 end

 def test_get_t_from_segment_single_edge
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-14
  assert_in_delta( 0.1, ge.segmentT(0.1), tol )
  assert_in_delta( 0.5, ge.segmentT(0.5), tol )
  assert_in_delta( 0.9, ge.segmentT(0.9), tol )
 end

 def test_get_t_from_segment_single_edge_endpoints
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-14
  assert_in_delta( 0.0, ge.segmentT(0.0), tol )
  assert_in_delta( 1.0, ge.segmentT(1.0), tol )
 end

 def test_get_t_from_segment_two_edges_out_of_range
  ge = GridEdger.new(two_edge_grid,1)
  assert_nil ge.segmentT(-1.0)
  assert_nil ge.segmentT(2.1)
 end

 def test_get_t_from_segment_two_edges
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-14
  assert_in_delta( 0.1, ge.segmentT(0.1), tol )
  assert_in_delta( 0.5, ge.segmentT(0.5), tol )
  assert_in_delta( 0.9, ge.segmentT(0.9), tol )
  assert_in_delta( 1.2, ge.segmentT(1.1), tol )
  assert_in_delta( 2.0, ge.segmentT(1.5), tol )
  assert_in_delta( 2.8, ge.segmentT(1.9), tol )
 end

 def test_get_t_from_segment_two_edges_endpoints
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-14
  assert_in_delta( 0.0, ge.segmentT(0.0), tol )
  assert_in_delta( 1.0, ge.segmentT(1.0), tol )
  assert_in_delta( 3.0, ge.segmentT(2.0), tol )
 end


 def test_get_m_from_two_segment_edge
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-14

  diag =  1.0
  m = [ diag, 0.0, 0.0, diag, 0.0, diag ];
  result = ge.segmentMap(0.0)
  6.times { |i| assert_in_delta( m[i], result[i], tol, "element #{i} error" ) }

  diag =  4.0
  m = [ diag, 0.0, 0.0, diag, 0.0, diag ];
  result = ge.segmentMap(1.0)
  6.times { |i| assert_in_delta( m[i], result[i], tol, "element #{i} error" ) }
  diag = 16.0

  m = [ diag, 0.0, 0.0, diag, 0.0, diag ];
  result = ge.segmentMap(2.0)
  6.times { |i| assert_in_delta( m[i], result[i], tol, "element #{i} error" ) }

  diag =  2.5
  m = [ diag, 0.0, 0.0, diag, 0.0, diag ];
  result = ge.segmentMap(0.5)
  6.times { |i| assert_in_delta( m[i], result[i], tol, "element #{i} error" ) }

  diag =  10.0
  m = [ diag, 0.0, 0.0, diag, 0.0, diag ];
  result = ge.segmentMap(1.5)
  6.times { |i| assert_in_delta( m[i], result[i], tol, "element #{i} error" ) }
 end

 def test_find_next_s_value_for_length_for_one_hlf_length_edge
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-14
  start = 0.0
  length = 2.0
  assert_in_delta( 1.0, ge.lengthToS(start,length), tol )
  start = 0.0
  length = 1.0
  assert_in_delta( 0.5, ge.lengthToS(start,length), tol )
 end

end
