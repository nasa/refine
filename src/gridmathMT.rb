#!/usr/bin/env ruby
#
# Mobility test for gridmath c lib
#
# $Id$

exit 1 unless system 'ruby makeRubyExtension.rb GridMath master_header.h'

require 'test/unit'
require 'GridMath/GridMath'

class TestGridMath < Test::Unit::TestCase

 def set_up
  @gm = GridMath.new
 end
 def setup ; set_up ; end

 def testSubtractVector
  v1 = [1,2,3]
  v2 = [7,0,6]
  assert_equal [-6,2,-3], @gm.subtractVector(v1,v2)
  v1 = [13,6,23]
  v2 = [2,10,16]
  assert_equal [11,-4,7], @gm.subtractVector(v1,v2)
 end

 def testDotProduct
  v1 = [1,2,3]
  v2 = [1,2,3]
  assert_equal 14, @gm.dotProduct(v1,v2)
  v1 = [2,0,5]
  v2 = [1,20,3]
  assert_equal 17, @gm.dotProduct(v1,v2)
 end

 def testCrossProduct
  v1 = [1,0,0]
  v2 = [0,1,0]
  assert_equal [0,0,1],  @gm.crossProduct(v1,v2)
  assert_equal [0,0,-1], @gm.crossProduct(v2,v1)
  assert_equal [0,0,0], @gm.crossProduct(v1,v1)
  assert_equal [0,0,0], @gm.crossProduct(v2,v2)
 end

 def testVectorLength
  assert_equal 5, @gm.vectorLength([3,4,0])
  assert_equal 0, @gm.vectorLength([0,0,0])
  assert_equal 1, @gm.vectorLength([1,0,0])
  assert_equal 10, @gm.vectorLength([0,10,0])
 end

 def testNormalizeVector
  assert_equal [1,0,0], @gm.vectorNormalize([5,0,0])
  assert_equal [0,1,0], @gm.vectorNormalize([0,6,0])
  assert_equal [0,0,-1], @gm.vectorNormalize([0,0,-3])
 end

end
