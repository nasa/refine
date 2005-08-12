#!/usr/bin/env ruby
#
# Mobility test for the ring c lib
#
# $Id$

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Ring').build

require 'test/unit'
require 'Ring/Ring'

class TestRing < Test::Unit::TestCase

 EMPTY = -1

 def set_up
  @ring = Ring.new
 end
 def setup ; set_up ; end

 def test_create_initializes_sizes_to_zero
  assert_equal 0, @ring.segments
  assert_equal 0, @ring.triangles
 end

 # segments

 def test_adding_segment_returns_ring
  assert_equal @ring, @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
 end

 def test_adding_new_segment_increments_total
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_equal 1, @ring.segments
 end

 def test_adding_existing_segment_returns_nil
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_nil @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
 end

 def test_adding_existing_segment_does_not_increment_total
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_equal 1, @ring.segments
 end

 def test_adding_mirror_segment_returns_ring
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_equal @ring, @ring.addSegment(4,2,[4.1,4.2],[2.1,2.2])
 end

 def test_adding_mirror_segment_removes_it_without_add
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addSegment(4,2,[4.1,4.2],[2.1,2.2])
  assert_equal 0, @ring.segments
 end

 # retrieving segment

 def test_retreve_segment_out_of_range_empty
  assert_nil @ring.segment(-1)
  assert_nil @ring.segment(0)
 end

 def test_retreve_segment_out_of_range_one
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_nil @ring.segment(1)
 end

 def test_retreve_segment_0
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  segment = @ring.segment(0)
  assert_equal 2, segment[0]
  assert_equal 4, segment[1]
  tol = 1.0e-15
  assert_in_delta 2.1, segment[2][0], tol
  assert_in_delta 2.2, segment[2][1], tol
  assert_in_delta 4.1, segment[3][0], tol
  assert_in_delta 4.2, segment[3][1], tol
 end

 # look for node in segments

 def test_do_not_find_node_that_does_not_exist_in_segements
  assert_nil @ring.segmentsContainNode(2)
 end

 def test_find_node_in_segement_returns_nodal_uv
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  tol = 1.0e-15
  assert_in_delta 2.1, @ring.segmentsContainNode(2)[0], tol
  assert_in_delta 2.2, @ring.segmentsContainNode(2)[1], tol

  assert_in_delta 4.1, @ring.segmentsContainNode(4)[0], tol
  assert_in_delta 4.2, @ring.segmentsContainNode(4)[1], tol
 end

 # evaluate canidate segment intersection with ring

 def add_right_ring_segments(n0,n1,n2)
  @ring.addSegment(0,1,n0,n1)
  @ring.addSegment(1,2,n1,n2)
  @ring.addSegment(2,0,n2,n0)
 end

 def test_ring_surrounds_a_completing_segment
  n0 = [0.0,0.0] ; n1 = [1.0,0.0] ; n2 = [0.0,1.0]
  add_right_ring_segments(n0,n1,n2)
  assert @ring.surroundsSegment(1,0,n1,n0)
  assert @ring.surroundsSegment(2,1,n2,n1)
  assert @ring.surroundsSegment(0,2,n0,n2)  
 end

 def test_ring_does_not_surrounds_a_duplicate_segment
  n0 = [0.0,0.0] ; n1 = [1.0,0.0] ; n2 = [0.0,1.0]
  add_right_ring_segments(n0,n1,n2)
  assert !@ring.surroundsSegment(0,1,n0,n1)
  assert !@ring.surroundsSegment(1,2,n1,n2)
  assert !@ring.surroundsSegment(2,0,n2,n0)  
 end

 # triangles

 def test_adding_triangle_returns_nil_if_base_is_not_a_segment
  assert_nil @ring.addTriangle(2,4,5,[5.1,5.2])
 end

 def test_adding_triangle_does_not_increment_total_if_base_is_not_a_segment
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_equal 0, @ring.triangles
 end

 def test_adding_triangle_returns_ring_if_base_is_a_segment
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_equal @ring, @ring.addTriangle(2,4,5,[5.1,5.2])
 end

 def test_adding_triangle_increments_total_if_base_is_a_segment
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_equal 1, @ring.triangles
 end

 # triangles add segments

 def test_adding_triangle_adds_three_sides_as_segments_creating_two
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_equal 2, @ring.segments
 end

 def test_adding_triangle_adds_three_sides_as_segments_creating_segment_1_2
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addSegment(4,5,[4.1,4.2],[5.1,5.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_equal 1, @ring.segments
 end

 def test_adding_triangle_adds_three_sides_as_segments_creating_segment_2_0
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addSegment(5,2,[5.1,5.2],[2.1,2.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_equal 1, @ring.segments
 end

 def test_adding_triangle_adds_three_sides_as_segments_creating_none
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addSegment(4,5,[4.1,4.2],[5.1,5.2])
  @ring.addSegment(5,2,[5.1,5.2],[2.1,2.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_equal 0, @ring.segments
 end

 def test_adding_triangle_adds_segments_with_correct_nodes_and_uvs
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])

  tol = 1.0e-15

  segment = @ring.segment(0)
  assert_equal 5, segment[0]
  assert_equal 4, segment[1]

  assert_in_delta 5.1, segment[2][0], tol
  assert_in_delta 5.2, segment[2][1], tol
  assert_in_delta 4.1, segment[3][0], tol
  assert_in_delta 4.2, segment[3][1], tol

  segment = @ring.segment(1)
  assert_equal 2, segment[0]
  assert_equal 5, segment[1]
  tol = 1.0e-15
  assert_in_delta 2.1, segment[2][0], tol
  assert_in_delta 2.2, segment[2][1], tol
  assert_in_delta 5.1, segment[3][0], tol
  assert_in_delta 5.2, segment[3][1], tol
 end

 # retrieving triangles

 def test_retreve_triangle_out_of_range_empty
  assert_nil @ring.triangle(-1)
  assert_nil @ring.triangle(0)
 end

 def test_retreive_triangle
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  @ring.addTriangle(2,4,5,[5.1,5.2])
  assert_nil @ring.triangle(1)

  triangle = @ring.triangle(0)
  assert_equal 2, triangle[0]
  assert_equal 4, triangle[1]
  assert_equal 5, triangle[2]
  tol = 1.0e-15
  assert_in_delta 2.1, triangle[3][0], tol
  assert_in_delta 2.2, triangle[3][1], tol
  assert_in_delta 4.1, triangle[4][0], tol
  assert_in_delta 4.2, triangle[4][1], tol
  assert_in_delta 5.1, triangle[5][0], tol
  assert_in_delta 5.2, triangle[5][1], tol
 end

end
