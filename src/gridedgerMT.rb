#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for gridedger c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('GridEdger').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMath/GridMath'
require 'GridMetric/GridMetric'
require 'GridEdger/GridEdger'

class Grid
 include GridMetric
end

class TestGridEdger < Test::Unit::TestCase

 EMPTY = (-1)
 
 def test_to_make_test_unit_happy
  assert true
 end

end
