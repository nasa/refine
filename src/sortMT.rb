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

 def testLengthZeroHeapSort
  assert_equal [], Sort.Heap([])
 end

 def testLengthOneHeapSort
  assert_equal [0], Sort.Heap([7])
 end

 def testLengthTwoHeapSort
  assert_equal [0,1], Sort.Heap([5,7])
 end

 def testLengthZeroBinarySearch
  assert_equal EMPTY, Sort.Search([],-1)
  assert_equal EMPTY, Sort.Search([],0)
  assert_equal EMPTY, Sort.Search([],1)
 end

end
