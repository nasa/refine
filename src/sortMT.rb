#!/usr/bin/env ruby
#
# Mobility test for sort c lib
#
# $Id$

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Sort').build

require 'test/unit'
require 'Sort/Sort'

class TestSort < Test::Unit::TestCase

 EMPTY = -1

 def testZeroLengthHeapSort
  assert_equal [], Sort.Heap([])
 end

 def testOneLengthHeapSort
  assert_equal [0], Sort.Heap([7])
 end

 def testTwoLengthHeapSort
  assert_equal [0,1], Sort.Heap([5,7])
 end

end
