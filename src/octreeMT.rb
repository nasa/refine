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

 def testStoreBoundingBox
  octree = Octree.new(-1,1,4,5,8,9)
  assert_equal [-1,1,4,5,8,9], octree.boundingBox
 end

end
