#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for adj c lib

exit unless system 'ruby extconf.rb Adj'
exit unless system 'make'

require 'test/unit'
require 'Adj'

class TestSampleUnit < Test::Unit::TestCase

 def set_up
  @adj = Adj.new(4,4)
 end

 def testCreateAdj
  assert_equal 4, @adj.nnode
 end

end
