#!/usr/bin/env ruby
#
# Mobility test for intersect c lib
#
# $Id$

exit 1 unless system 'ruby makeRubyExtension.rb Intersect master_header.h'

require 'test/unit'
require 'Intersect/Intersect'

class TestIntersect < Test::Unit::TestCase

 def testNodeSide
  intersect = Intersect.new
  v0 = [0,0,0]
  v1 = [1,0,0]
  v2 = [0,1,0]
  nn = [0,0,-1]
  n0 = [0,0,0]
  np = [0,0,1]
  assert_equal( -1, intersect.side(v0,v1,v2,nn))
  assert_equal(  0, intersect.side(v0,v1,v2,n0))
  assert_equal(  1, intersect.side(v0,v1,v2,np))
 end

end
