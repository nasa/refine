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
  @bigNode = 10000
 end

 def testCreate
  assert_equal 4, @adj.nnode
 end

 def testRegister
  assert_not_nil @adj.register(0,1)
  assert_nil @adj.register(4,1)
  assert_nil @adj.register(@bigNode,1)
 end

 def testIterator
  assert_equal false, @adj.valid
  assert_equal false, @adj.more
  assert_equal( -1, @adj.current)
  
  assert_nil @adj.first(@bigNode);
  assert_equal @adj, @adj.first(0);
  assert_equal false, @adj.valid;
  
  assert_equal @adj, @adj.register( 2, 299 )
  assert_equal @adj, @adj.first(2)
  assert_equal 299,  @adj.current
  assert_equal true,  @adj.valid
  assert_equal false, @adj.more
  assert_equal @adj, @adj.next
  assert_equal false, @adj.valid
  
  assert_equal @adj, @adj.register( 3, 398 )
  assert_equal @adj, @adj.register( 3, 399 )
  @adj.first(3);
  assert_equal true,  @adj.valid
  assert_equal true,  @adj.more
  
  100.times {@adj.next} # abusive use of next
 end

end
