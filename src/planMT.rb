#!/usr/bin/env ruby
#
# Mobility test for the plan c lib
#
# $Id$

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Plan').build

require 'test/unit'
require 'Plan/Plan'

class TestQueue < Test::Unit::TestCase

 EMPTY = -1

# create plan with intial size and growth stride
# a plan has item
#   planAddItemWithPriority
#   planItems
# each item has a priority
#   planItemPriorty
# the items priority is sorted to give its ranking 
#   planItemWithThisRanking

 def test_create_initializes_sizes
  max_size = 273
  chunk    = 1527
  plan = Plan.new(max_size,chunk)
  assert_equal 0,        plan.size
  assert_equal max_size, plan.max_size
  assert_equal chunk,    plan.chunk_size
 end

 def test_create_initializes_minimum_sizes
  plan = Plan.new(0,0)
  assert_equal 0, plan.size
  assert_equal 1, plan.max_size
  assert_equal 1, plan.chunk_size
 end

 def test_create_initializes_with_garbage_input
  plan = Plan.new(-53,-1)
  assert_equal 0, plan.size
  assert_equal 1, plan.max_size
  assert_equal 1, plan.chunk_size
 end

 def test_extra_items_resize_plan_by_chunk
  assert_equal 0, 1
 end

 def test_rankings_derived_from_priorities
  plan = Plan.new(6,6)
  assert_equal EMPTY, plan.itemWithThisRanking(0)
  assert_nil plan.deriveRankingsFromPriorities
  assert_equal EMPTY, plan.itemWithThisRanking(0)
  plan.addItemWithPriority(0,2.0)
  plan.addItemWithPriority(1,1.0)
  plan.addItemWithPriority(2,0.0)
  plan.addItemWithPriority(3,3.0)
  plan.addItemWithPriority(4,4.0)
  plan.addItemWithPriority(5,5.0)
  assert_equal EMPTY, plan.itemWithThisRanking(0)
  assert_equal grid, plan.deriveRankingsFromPriorities
  assert_equal 2, plan.itemWithThisRanking(0)
  assert_equal 1, plan.itemWithThisRanking(1)
  assert_equal 0, plan.itemWithThisRanking(2)
  3.upto(5) do |rank|
   assert_equal rank, plan.itemWithThisRanking(rank)
  end
  plan.addItemWithPriority(5,-5.0)
  assert_equal grid, plan.deriveRankingsFromPriorities
  assert_equal 5, plan.itemWithThisRanking(0)
 end


end
