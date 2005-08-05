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

 def test_adding_segment_returns_ring
  assert_equal @ring, @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
 end

 def test_adding_segment_increments_total
  @ring.addSegment(2,4,[2.1,2.2],[4.1,4.2])
  assert_equal 1, @ring.segments
 end

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

 puts 'test_adding_triangle_also_adds_three_sements'
 puts 'adding segment in opposite direction removes it from ring'

end
