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

 UnitBoundingBox = [0,1, 0,1, 0,1]

 def testInitBoundingBox
  boundingBox = [-1,1,4,5,8,9]
  octree = Octree.new(boundingBox)
  assert_equal boundingBox, octree.boundingBox
  unitOctree = Octree.new(UnitBoundingBox)
  assert_equal UnitBoundingBox, unitOctree.boundingBox
 end

 def testInitOneOctant
  octree = Octree.new(UnitBoundingBox)
  assert_equal 1, octree.nOctant
 end

 def testQwertyEmptyOctant
  octree = Octree.new(UnitBoundingBox)
  identity = [1,0,0, 0,1,0, 0,0,1]
  center = [0.5,0.5,0.5]
  assert_equal identity, octree.query(center)
 end

 def testQwertyConstantOctant
  octree = Octree.new(UnitBoundingBox)
  data   = [1,2,3, 4,5,6, 7,8,9]
  center = [0.5,0.5,0.5]
  assert_equal octree, octree.insert(center,data)
  assert_equal data, octree.query(center)
 end

 def testTwoOctant2
  octree = Octree.new(UnitBoundingBox)
  rightData = [0,0,0, 0,0,0, 0,0,0]
  rightNode = [0.7, 0.5, 0.5]
  assert_equal octree, octree.insert(rightNode, rightData)
  leftData = [1,1,1, 1,1,1, 1,1,1]
  leftNode = [0.2, 0.5, 0.5]
  assert_equal octree, octree.insert(leftNode, leftData)

  assert_equal leftData, octree.query(leftNode)
  assert_equal rightData, octree.query(rightNode)
 end

end
