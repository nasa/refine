#!/usr/bin/env ruby
#
# Mobility test for gridmath c lib
#
# $Id$

exit 1 unless system 'ruby makeRubyExtension.rb GridMath master_header.h'

require 'test/unit'
require 'GridMath/GridMath'

class TestGridMath < Test::Unit::TestCase

 def testSubtractVector
  gm = GridMath.new
  v1 = [1,2,3]
  v2 = [7,0,6]
  assert_equal [-6,2,-3], gm.subtractVector(v1,v2)
  v1 = [13,6,23]
  v2 = [2,10,16]
  assert_equal [11,-4,7], gm.subtractVector(v1,v2)
 end

end
