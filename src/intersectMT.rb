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

 def testNodeSide
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  nn = [0,0,-1]
  n0 = [0,0,0]
  np = [0,0,1]
  assert_equal( -1, @intersect.side(v0,v1,v2,nn))
  assert_equal(  0, @intersect.side(v0,v1,v2,n0))
  assert_equal(  1, @intersect.side(v0,v1,v2,np))
 end

 def testPlaneAndSegemntThroughMiddle
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  n0 = [0.3,0.3,1]
  n1 = [0,3,0.3,1]
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

end
