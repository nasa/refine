#!/usr/bin/env ruby
#
# Mobility test for near c lib
#
# $Id$

exit 1 unless system 'ruby makeRubyExtension.rb Near master_header.h'

require 'test/unit'
require 'Near/Near'

class TestNear < Test::Unit::TestCase

 def testStoreIndex
  near = Near.new(5,0,0,0,0)
  assert_equal 5, near.index
 end

 def testNodePointToPointClearance
  point1 = Near.new(5,0,0,0,0)
  point2 = Near.new(5,1,0,0,0)
  assert_equal 1, point1.clearance(point2)

  point2 = Near.new(5,2,0,0,0)
  assert_equal 2, point1.clearance(point2)

  point2 = Near.new(5,0,1,0,0)
  assert_equal 1, point1.clearance(point2)

  point2 = Near.new(5,0,0,1,0)
  assert_equal 1, point1.clearance(point2)
 end

 def testNodeClearanceWithSphereRadius
  sphere1 = Near.new(5,0,0,0,1)
  sphere2 = Near.new(5,1,0,0,1)
  assert_equal 0, sphere1.clearance(sphere2)

  sphere2 = Near.new(5,3,0,0,1)
  assert_equal 1, sphere1.clearance(sphere2)
 end
 
 def testInitializeWithNoChildren
  near = Near.new(5,0,0,0,0)
  assert_equal( -1, near.rightIndex)
  assert_equal( -1, near.leftIndex)
 end
 
 def testPopulateLeftChildFirst
  near  = Near.new(5,0,0,0,0)
  child = Near.new(6,1,0,0,0)
  assert_equal near, near.insert(child)
  assert_equal(-1, near.rightIndex)
  assert_equal 6, near.leftIndex
 end
 
 def testPopulateRightChildSecond
  near  = Near.new(5,0,0,0,0)
  child6 = Near.new(6,1,0,0,0)
  child7 = Near.new(7,1,0,0,0)
  assert_equal near, near.insert(child6)
  assert_equal near, near.insert(child7)
  assert_equal 7, near.rightIndex
  assert_equal 6, near.leftIndex
 end
 
 def testPopulateLeftCloser
  near  = Near.new(5,0,0,0,0)
  child6 = Near.new(6,1,0,0,0)
  child7 = Near.new(7,2,0,0,0)
  child8 = Near.new(8,1,0,0,0)
  assert_equal near, near.insert(child6)
  assert_equal near, near.insert(child7)
  assert_equal near, near.insert(child8)
  assert_equal 7, near.rightIndex
  assert_equal 6, near.leftIndex
  assert_equal 8, child6.leftIndex
 end
 
 def testPopulateRightCloser
  near  = Near.new(5,0,0,0,0)
  child6 = Near.new(6,1,0,0,0)
  child7 = Near.new(7,2,0,0,0)
  child8 = Near.new(8,2,0,0,0)
  assert_equal near, near.insert(child6)
  assert_equal near, near.insert(child7)
  assert_equal near, near.insert(child8)
  assert_equal 7, near.rightIndex
  assert_equal 6, near.leftIndex
  assert_equal 8, child7.leftIndex
 end
 
 def testInitializeChildDistance
  near  = Near.new(5,0,0,0,0)
  assert_equal 0, near.farChild
  assert_equal 0, near.rightDistance
  assert_equal 0, near.leftDistance
 end

end
