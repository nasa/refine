#!/usr/bin/env ruby
#

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Adj').build

require 'test/unit'
require 'Adj/Adj'

class TestAdj < Test::Unit::TestCase

 def set_up
  @adj = Adj.new(4,6,3)
  @bigNode = 10000
 end
 def setup ; set_up ; end

 def testCreate
  assert_equal 4, @adj.nnode
  assert_equal 6, @adj.nadj
  assert_equal 3, @adj.chunkSize
 end

 def testCreateZeroSize
  localAdj = Adj.new(0,0,0)
  assert_equal 1, localAdj.nnode
  assert_equal 1, localAdj.nadj
  assert_equal 1, localAdj.chunkSize
 end

 def testRegister
  assert_nil @adj.register(-1,1)
  assert_not_nil @adj.register(0,1)
  assert_nil @adj.register(4,1)
  assert_nil @adj.register(@bigNode,1)
 end

 def testRegisterAllocatesNewCunk
  6.times {|i| assert_not_nil @adj.register(0,i)}
  assert_equal 6, @adj.nadj
  assert_not_nil  @adj.register(1,0)
  assert_equal 9, @adj.nadj
 end

 def testReallocKeepsOldRegisteredElements
  degree = 10
  4.times { |node| degree.times { |i| @adj.register(node,10*i+node) } }
  assert_equal degree, @adj.degree(0) 
  assert_equal degree, @adj.degree(1) 
  assert_equal degree, @adj.degree(2) 
  assert_equal degree, @adj.degree(3) 
 end

 def testReallocateNNodesUp
  assert_equal 4, @adj.nnode
  assert_not_nil  @adj.realloc(9)
  assert_equal 9, @adj.nnode
  assert_equal @adj,  @adj.register( 8, 899 )
  assert_not_nil  @adj.realloc(0)
  assert_equal 1, @adj.nnode
 end

 def testReallocateDownWithoutMemoryLeak
  assert_equal @adj, @adj.register( 2, 299 )
  assert_not_nil  @adj.realloc(0)
  assert_equal 1, @adj.nnode
  6.times {|i| assert_not_nil @adj.register(0,i)}
  assert_equal  6, @adj.nadj
 end

 def testIteratorInitialState
  assert_equal false, @adj.valid
  assert_equal false, @adj.more
  assert_equal( -1, @adj.item)
 end
  
 def testIteratorForBigNode
  assert_equal @adj,  @adj.first(@bigNode);
  assert_equal false, @adj.valid;
  assert_equal @adj,  @adj.first(-1);
  assert_equal false, @adj.valid;
 end

 def testIterator1Element
  assert_equal @adj,  @adj.register( 2, 299 )
  assert_equal @adj,  @adj.first(2)
  assert_equal 299,   @adj.item
  assert_equal true,  @adj.valid
  assert_equal false, @adj.more
  assert_equal @adj,  @adj.next
  assert_equal false, @adj.valid
 end

 def testIterator2Elements
  assert_equal @adj,  @adj.register( 3, 398 )
  assert_equal @adj,  @adj.register( 3, 399 )
  assert_equal @adj,  @adj.first(3);
  assert_equal true,  @adj.valid
  assert_equal 399,   @adj.item
  assert_equal true,  @adj.more
  assert_equal @adj,  @adj.next
  assert_equal true,  @adj.valid
  assert_equal 398,   @adj.item
  assert_equal false, @adj.more
 end
 
 def testAbusiveUseOfNext
  100.times {@adj.next}
  assert_equal @adj,  @adj.register( 1, 498 )
  assert_equal @adj,  @adj.register( 1, 499 )
  assert_equal @adj,  @adj.first(3)
  100.times {@adj.next}
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
  assert_equal @adj, @adj.register( 0, 298 )
  assert_equal    2, @adj.degree(0)
 end

end
