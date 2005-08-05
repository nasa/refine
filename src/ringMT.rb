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

 # retrieving triangles


end
