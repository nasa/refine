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

# a plan has item
#   planAddItemWithPriority
#   planItems
# each item has a priority
#   planItemPriorty
# the items priority is sorted to give its ranking 
#   planItemWithThisRanking

end
