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


end
