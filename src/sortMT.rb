#!/usr/bin/env ruby
#
# Mobility test for sort c lib
#
# $Id$

Dir.chdir ENV['srcdir'] if ENV['srcdir']
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

 def testLengthThreeHeapSortRepeated220
  assert_equal [2,1,0], Sort.Heap([1234,1234,2])
 end

 def testLotsOfReversedNumbers
  n=5000
  input = []
  n.times { |index| input.push n-index }
  output = Sort.Heap input
  (n-1).times { |index| assert(input[output[index]]<=input[output[index+1]])}
 end

 def testLotsOfDuplicateRandomNumbers1
  srand 12098245
  n=5000
  input = []
  n.times { |index| input.push rand(n/2) }
  output = Sort.Heap input
  (n-1).times { |index| assert(input[output[index]]<=input[output[index+1]])}
 end

 def testLotsOfDuplicateRandomNumbers2
  srand 24589879
  n=5000
  input = []
  n.times { |index| input.push rand(n/2) }
  output = Sort.Heap input
  (n-1).times { |index| assert(input[output[index]]<=input[output[index+1]])}
 end

 def testLotsOfDuplicateRandomNumbers3
  srand 339040234
  n=5000
  input = []
  n.times { |index| input.push rand(n/2) }
  output = Sort.Heap input
  (n-1).times { |index| assert(input[output[index]]<=input[output[index+1]])}
 end

 def testLengthZeroBinarySearch
  assert_equal EMPTY, Sort.Search([],-1)
  assert_equal EMPTY, Sort.Search([],0)
  assert_equal EMPTY, Sort.Search([],1)
 end

 def testLengthOneBinarySearch
  assert_equal EMPTY, Sort.Search([7],3)
  assert_equal EMPTY, Sort.Search([7],10)
  assert_equal 0,     Sort.Search([5],5)
 end

 def testBinarySearch2
  sorted = [1,10]
  assert_equal 0,     Sort.Search(sorted,1)
  assert_equal EMPTY, Sort.Search(sorted,5)
  assert_equal 1,     Sort.Search(sorted,10)
 end

 def testBinarySearch3
  sorted = [1,5,10]
  assert_equal 0,     Sort.Search(sorted,1)
  assert_equal EMPTY, Sort.Search(sorted,2)
  assert_equal 1,     Sort.Search(sorted,5)
  assert_equal EMPTY, Sort.Search(sorted,7)
  assert_equal 2,     Sort.Search(sorted,10)
 end


 def testBinarySearch5
  sorted = [1,2,3,4,5]
  assert_equal 0, Sort.Search(sorted,1)
  assert_equal 1, Sort.Search(sorted,2)
  assert_equal 2, Sort.Search(sorted,3)
  assert_equal 3, Sort.Search(sorted,4)
  assert_equal 4, Sort.Search(sorted,5)
 end

 def testBinarySearch6
  sorted = [1,2,3,4,5,6]
  assert_equal 0, Sort.Search(sorted,1)
  assert_equal 1, Sort.Search(sorted,2)
  assert_equal 2, Sort.Search(sorted,3)
  assert_equal 3, Sort.Search(sorted,4)
  assert_equal 4, Sort.Search(sorted,5)
  assert_equal 5, Sort.Search(sorted,6)
 end

end
