#!/usr/bin/env ruby
#
# Mobility test for octree c lib
#
# $Id$

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Octree').build

require 'test/unit'
require 'Octree/Octree'

class TestOctree < Test::Unit::TestCase

 EMPTY = -1

 def testInitBoundingBox
  octree = Octree.new(-1,1,4,5,8,9)
  assert_equal [-1,1,4,5,8,9], octree.boundingBox
 end

 def testAddDataToEmptyOctant
  octree = Octree.new(0,1,0,1,0,1)
  assert_equal 1, octree.nOctant
  identity = [1,0,0,1,0,1]
  node = [0.5,0.5,0.5]
  data = [1,2,3,4,5,6]  
 end

end
