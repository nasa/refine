#!/usr/bin/env ruby
#
# Mobility test for queue c lib
#
# $Id$

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Queue').build

require 'test/unit'
require 'Queue/Queue'

class TestQueue < Test::Unit::TestCase

 EMPTY = -1

 def testInitializedWithOneTransaction
  assert_equal 1, Queue.new.transactions
 end

 def testNewTransactionIncrementsTransactions
  q = Queue.new
  assert_equal q, q.newTransaction
  assert_equal 2, q.transactions
  assert_equal q, q.newTransaction
  assert_equal 3, q.transactions
 end

 def testResetSetsTheNumberOfTransactionsToOne
  q = Queue.new.newTransaction.newTransaction
  assert_equal q, q.reset
  assert_equal 1, q.transactions
 end

 def testTransactionInitializedWithZeroRemoved
  q = Queue.new
  assert_equal EMPTY, q.removed(-1)
  assert_equal 0,     q.removed(0)
  assert_equal EMPTY, q.removed(1)
  q.newTransaction
  assert_equal 0,     q.removed(0)
  assert_equal 0,     q.removed(1)
 end

 def testTransactionInitializedWithZeroAdded
  q = Queue.new
  assert_equal EMPTY, q.added(-1)
  assert_equal 0,     q.added(0)
  assert_equal EMPTY, q.added(1)
  q.newTransaction
  assert_equal 0,     q.added(0)
  assert_equal 0,     q.added(1)
 end

 def testRemoveTransaction
  nodes = [0,1,2,3]
  q = Queue.new
  assert_equal q,     q.remove(nodes)
  assert_equal 1,     q.removed(0)
  assert_equal nodes, q.removedNodes(0)
 end

 def testAddTransaction
  cellId = 7
  nodes = [0,1,2,3,cellId]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new
  assert_equal q,     q.add(nodes,xyzs)
  assert_equal 1,     q.added(0)
  assert_equal nodes, q.addedNodes(0)
  assert_equal xyzs,  q.addedXYZs(0)
 end

 def testResetZerosOutPreviouslyRemovedNodes
  nodes = [0,1,2,3]
  q = Queue.new.remove([9,9,9,9]).reset
  assert_equal 0, q.removed(0)
  q.remove(nodes)
  assert_equal nodes, q.removedNodes(0)  
 end

 def testReallocMemoryForLotsOfTransactions
  q = Queue.new
  10000.times { q.newTransaction.remove([0,1,2,3]) }
 end

end
