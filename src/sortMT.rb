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

 def testLengthTwoHeapSortOrdered
  assert_equal [0,1], Sort.Heap([5,7])
 end

 def testLengthTwoHeapSortReversed
  assert_equal [1,0], Sort.Heap([67,2])
 end

 def testLengthThreeHeapSort012
  assert_equal [0,1,2], Sort.Heap([5,7,9])
 end

 def testLengthThreeHeapSort021
  assert_equal [0,2,1], Sort.Heap([5,10,9])
 end

 def testLengthThreeHeapSort102
  assert_equal [1,0,2], Sort.Heap([2,1,3])
 end

 def testLengthThreeHeapSort120
  assert_equal [1,2,0], Sort.Heap([17,5,15])
 end

 def testLengthThreeHeapSort201
  assert_equal [2,0,1], Sort.Heap([234,465,2])
 end

 def testLengthThreeHeapSort210
  assert_equal [2,1,0], Sort.Heap([1234,465,2])
 end

 def XtestLengthZeroBinarySearch
  assert_equal EMPTY, Sort.Search([],-1)
  assert_equal EMPTY, Sort.Search([],0)
  assert_equal EMPTY, Sort.Search([],1)
 end

end
