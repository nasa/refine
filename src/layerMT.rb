#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb Layer grid.h master_header.h'

require 'test/unit'
require 'Grid/Grid'
require 'Layer/Layer'

class TestLayer < Test::Unit::TestCase

 def testInit
  grid = Grid.new(2,0,0,0)
  layer = Layer.new(grid)
  assert_equal 2, layer.maxnode
  assert_equal 0, layer.nfront
 end

end
