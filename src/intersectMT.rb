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
  tri0 = [0,0,0]
  tri1 = [1,0,0]
  tri2 = [0,1,0]
  noden = [0,0,-1]
  node0 = [0,0,0]
  nodep = [0,0,1]
  assert_equal( -1, intersect.side(tri0,tri1,tri2,noden))
  assert_equal(  0, intersect.side(tri0,tri1,tri2,node0))
  assert_equal(  1, intersect.side(tri0,tri1,tri2,nodep))
 end

end
