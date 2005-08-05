#!/usr/bin/env ruby
#
# Mobility test for the plan c lib
#
# $Id$

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Ring').build

require 'test/unit'
require 'Ring/Ring'

class TestQueue < Test::Unit::TestCase

 EMPTY = -1

 def set_up
  @ring = Ring.new
 end
 def setup ; set_up ; end

 def test_create_initializes_sizes_to_zero
  assert_equal 0,        @ring.segments
  assert_equal 0,        @ring.triangles
 end

end
