#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for adj c lib

exit 1 unless system 'ruby makeRubyExtension.rb Adj master_header.h'

require 'test/unit'
require 'Adj/Adj'

class TestAdj < Test::Unit::TestCase

 def set_up
  @adj = Adj.new(4,4)
  @bigNode = 10000
 end

 def testCreate
  assert_equal 4, @adj.nnode
 end

 def testRegister
  assert_nil @adj.register(-1,1)
  assert_not_nil @adj.register(0,1)
  assert_nil @adj.register(4,1)
  assert_nil @adj.register(@bigNode,1)
 end

 def testRegisterMax
  4.times {|i| assert_not_nil @adj.register(0,i)}
  assert_nil @adj.register(1,0)
 end

 def testIterator
  assert_equal false, @adj.valid
  assert_equal false, @adj.more
  assert_equal( -1, @adj.item)
  
  assert_equal @adj,  @adj.first(@bigNode);
  assert_equal false, @adj.valid;
  assert_equal @adj,  @adj.first(-1);
  assert_equal false, @adj.valid;
  
  assert_equal @adj,  @adj.register( 2, 299 )
  assert_equal @adj,  @adj.first(2)
  assert_equal 299,   @adj.item
  assert_equal true,  @adj.valid
  assert_equal false, @adj.more
  assert_equal @adj,  @adj.next
  assert_equal false, @adj.valid
  
  assert_equal @adj,  @adj.register( 3, 398 )
  assert_equal @adj,  @adj.register( 3, 399 )
  assert_equal @adj,  @adj.first(3);
  assert_equal true,  @adj.valid
  assert_equal true,  @adj.more
  assert_equal @adj,  @adj.next
  assert_equal true,  @adj.valid
  assert_equal false, @adj.more
  
  100.times {@adj.next} # abusive use of next
 end

 def testAddAndRemove
  assert_equal false, @adj.exists(1,0)
  assert_nil          @adj.remove(1,0)
  assert_equal @adj,  @adj.register(1,0)
  assert_equal true,  @adj.exists(1,0)
  assert_equal @adj,  @adj.remove(1,0)
  assert_equal false, @adj.exists(1,0)
  assert_nil          @adj.remove(1,0)
 end
 
 def testMultipleExists
  assert_equal false, @adj.exists(1,198)
  assert_equal @adj,  @adj.register(1,198)
  assert_equal @adj,  @adj.register(2,198)
  assert_equal @adj,  @adj.register(1,199)
  
  assert_equal true,  @adj.exists(1,198)
  assert_equal true,  @adj.exists(1,199)
  @adj.remove(1,198)
  assert_equal false, @adj.exists(1,198)
  assert_equal true,  @adj.exists(1,199)
  @adj.register(1,198)
  assert_equal true,  @adj.exists(1,198)
  assert_equal true,  @adj.exists(1,199)
 end

 def testDegree
  assert_equal    0, @adj.degree(0)
  assert_equal @adj, @adj.register( 0, 299 )
  assert_equal    1, @adj.degree(0)
 end

end
