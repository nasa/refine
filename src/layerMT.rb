#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

exit 1 unless system 'ruby makeRubyExtension.rb Layer master_header.h'

require 'test/unit'
require 'Layer/Layer'

class TestLayer < Test::Unit::TestCase

 def testInit
  l = Layer.new
  assert_equal 0, l.nfront
 end

end
