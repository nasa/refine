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
  grid = Grid.new(5,2,2,1)
  grid.addNode(0.0,0.0,0.0)
  grid.addNode(1.0,0.0,0.0)
  grid.addNode(0.5,0.5,0.0)
  grid.addNode(0.5,-0.5,0.0)
  grid.addNode(0.5,0.0,0.5)

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

  grid.addCell( 0, 1, 2, 4)
  grid.addCell( 1, 0, 3, 4)

  grid.setNGeomNode(2)
  grid.setNGeomEdge(1)
  grid.setNGeomFace(2)

  grid.addGeomEdge(edgeId,0,1)

  grid.nnode.times do |node|
   grid.scaleSpacing(node,0.5)
  end

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

 def test_initialize_nodes_to_zero
  ge = GridEdger.new(Grid.new(0,0,0,0),5)
  assert_equal( 0, ge.nodes )
 end

 def test_empty_initial_edger_doesnt_return_s
  ge = GridEdger.new(Grid.new(0,0,0,0),5)
  assert_nil( ge.nodeS(0) )
 end

 def test_get_discrete_segment_and_ratio_from_s_out_of_range
  ge = GridEdger.new(two_edge_grid,1)
  assert_nil ge.discreteSegmentAndRatio(-1.0)
  assert_nil ge.discreteSegmentAndRatio( 3.0)
 end

 def test_get_discrete_segment_and_ratio_from_s_in_range
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-14

  result = ge.discreteSegmentAndRatio(0.0)
  truth = [ 0, 0.0 ]
  assert_equal( truth[0], result[0] )
  assert_in_delta( truth[1], result[1], tol )

  result = ge.discreteSegmentAndRatio(0.5)
  truth = [ 0, 0.5 ]
  assert_equal( truth[0], result[0] )
  assert_in_delta( truth[1], result[1], tol )

  result = ge.discreteSegmentAndRatio(1.0)
  truth = [ 1, 0.0 ]
  assert_equal( truth[0], result[0] )
  assert_in_delta( truth[1], result[1], tol )

  result = ge.discreteSegmentAndRatio(1.5)
  truth = [ 1, 0.5 ]
  assert_equal( truth[0], result[0] )
  assert_in_delta( truth[1], result[1], tol )

  result = ge.discreteSegmentAndRatio(2.0)
  truth = [ 1, 1.0 ]
  assert_equal( truth[0], result[0] )
  assert_in_delta( truth[1], result[1], tol )
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

 def test_find_next_s_value_for_length_for_one_half_length_edge_out_of_range
  ge = GridEdger.new(one_edge_grid,1)
  start = -1.0
  assert_nil ge.lengthToS(start,1.0)
  start = 5.0
  assert_nil ge.lengthToS(start,1.0)
 end

 def test_find_next_s_value_for_length_for_one_half_length_edge
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-12

  start = 0.0
  length = 2.0
  assert_in_delta( 1.0, ge.lengthToS(start,length), tol )

  start = 0.0
  length = 1.0
  assert_in_delta( 0.5, ge.lengthToS(start,length), tol )

  start = 0.5
  length = 0.5
  assert_in_delta( 0.75, ge.lengthToS(start,length), tol )
 end

 def test_find_next_s_value_for_length_for_end_point
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-12
  start = 0.9
  length = 1.0
  assert_in_delta( 1.0, ge.lengthToS(start,length), tol )
 end

 def test_find_next_s_value_for_length_for_two_edge
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-8

  start = 0.0
  length = 1.0
  target = 0.69874993561325
  assert_in_delta( target, ge.lengthToS(start,length), tol )

  start = 0.0
  length = 3.0
  target = 1.25
  assert_in_delta( target, ge.lengthToS(start,length), tol )

  start = 0.75
  length = 1.0
  target = 1.11560137232357
  assert_in_delta( target, ge.lengthToS(start,length), tol )
 end

 def test_discretize_single_edge_with_two_nodes
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-12

  length = 2.0
  assert_equal ge, ge.discretize(length)
  assert_equal 2, ge.nodes

  assert_in_delta( 0.0, ge.nodeS(0), tol )
  assert_in_delta( 1.0, ge.nodeS(1), tol )
 end

 def test_discretize_single_edge_with_three_nodes
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-12

  length = 1.0
  assert_equal ge, ge.discretize(length)
  assert_equal 3, ge.nodes

  assert_in_delta( 0.0, ge.nodeS(0), tol )
  assert_in_delta( 0.5, ge.nodeS(1), tol )
  assert_in_delta( 1.0, ge.nodeS(2), tol )
 end

 def test_discretize_single_edge_with_2_plus_nodes
  ge = GridEdger.new(one_edge_grid,1)
  tol = 1.0e-12

  length = 1.5
  assert_equal ge, ge.discretize(length)
  assert_equal 3, ge.nodes

  assert_in_delta( 0.00, ge.nodeS(0), tol )
  assert_in_delta( 0.75, ge.nodeS(1), tol )
  assert_in_delta( 1.00, ge.nodeS(2), tol )
 end

 def test_discretize_double_edge_with_a_few_nodes
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-12

  length = 2.0
  assert_equal ge, ge.discretize(length)
  assert_equal 5, ge.nodes

  assert_in_delta( 0.0000000000000, ge.nodeS(0), tol )
  assert_in_delta( 1.0795481390678, ge.nodeS(1), tol )
  assert_in_delta( 1.4524358428491, ge.nodeS(2), tol )
  assert_in_delta( 1.7509707149822, ge.nodeS(3), tol )
  assert_in_delta( 2.0000000000000, ge.nodeS(4), tol )
 end

 def test_discretize_double_edge_by_iterating_length
  ge = GridEdger.new(two_edge_grid,1)
  tol = 1.0e-12

  assert_equal ge, ge.discretizeEvenly
 end

 def test_discretize_single_edge_into_two_edges_via_insert
  grid = one_edge_grid
  ge = GridEdger.new(grid,1)

  assert_nil ge.insert
  length = 1.0
  ge.discretize(length)

  assert_equal ge, ge.insert

  assert_equal 6, grid.nnode
 end

end
