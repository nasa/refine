#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for gridshape c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'Header2Wrapper'
Header2Wrapper.new('GridShape').convert

require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('GridShape').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridShape/GridShape'

class Grid
 include GridShape
end

class TestGridShape < Test::Unit::TestCase

 def testStupidity
  assert true
  assert false
 end

end
