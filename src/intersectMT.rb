#!/usr/bin/env ruby
#
# Mobility test for intersect c lib
#
# $Id$

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Intersect').build

require 'test/unit'
require 'Intersect/Intersect'

class TestIntersect < Test::Unit::TestCase

 def set_up
  @intersect = Intersect.new
 end
 def setup ; set_up ; end

 def testPlaneAndNodeInAndOut
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  assert_equal true, @intersect.triangleNode(v0,v1,v2,[0.3,0.3,0])
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[0,-0.5,0])
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[1,1,0])
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[-0.5,0,0])
 end

 def testPlaneAndNodeNick
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[0,0,0])
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[1,0,0])
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[0,1,0])
  assert_equal false, @intersect.triangleNode(v0,v1,v2,[0.5,0.5,0])
 end

 def testPlaneAndSegmentThroughMiddle
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  n0 = [0.3,0.3,1]
  n1 = [0.3,0.3,1]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [0.3,0.3,-1]
  assert_equal true, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [0.3,0.3,0]
  assert_equal true, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n1 = [0.3,0.3,0]
  assert_equal true, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [0.3,0.3,-0]
  assert_equal true, @intersect.triangleSegment(v0,v1,v2,n0,n1)
 end

 def testPlaneAndSegmentThroughNick
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  n0 = [0.5,0.5,10]
  n1 = [0.5,0.5,-1]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [0,0,18]
  n1 = [0,0,-3]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
 end

 def testPlaneAndSegmentThroughMiss
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  n0 = [1,1,10]
  n1 = [1,1,-1]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [0,-1,18]
  n1 = [0,-1,-3]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
 end

 def testPlaneAndSegmentAdjecent
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  n0 = [0,0,0]
  n1 = [1,0,0]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [1,0,0]
  n1 = [0,1,0]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [0,1,0]
  n1 = [0,0,0]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
 end

 def testPlaneAndSegmentCoplanarPerce
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  n0 = [0.3,0.3,0]
  n1 = [1,1,0]
  assert_equal true, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [1,1,0]
  n1 = [0.3,0.3,0]
  assert_equal true, @intersect.triangleSegment(v0,v1,v2,n0,n1)
  n0 = [1,1,0]
  n1 = [2,2,0]
  assert_equal false, @intersect.triangleSegment(v0,v1,v2,n0,n1)
 end

end
