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

 def testInitializedWithNoTransactions
  assert_equal 0, Queue.new.transactions
 end

 def testNewTransactionIncrementsTransactions
  q = Queue.new
  assert_equal q, q.newTransaction
  assert_equal 1, q.transactions
  assert_equal q, q.newTransaction
  assert_equal 2, q.transactions
 end

 def testResetZerosOutTheNumberOfTransactions
  q = Queue.new.newTransaction.newTransaction
  assert_equal q, q.reset
  assert_equal 0, q.transactions
 end

 def testTransactionInitializedWithZeroNodes
  q = Queue.new
  assert_equal EMPTY, q.transactionNodes(-1)
  assert_equal 0,     q.transactionNodes(0)
  assert_equal EMPTY, q.transactionNodes(1)
  q.newTransaction
  assert_equal 0,     q.transactionNodes(0)
  assert_equal 0,     q.transactionNodes(1)
 end

 def testAddNodeToTransaction
  nodes=[23,65,78]
  q = Queue.new
  assert_equal q, q.addNode(nodes[0])
  q.newTransaction
  assert_equal q, q.addNode(nodes[1])
  assert_equal q, q.addNode(nodes[2])
  assert_equal 1, q.transactionNodes(0)
  assert_equal 2, q.transactionNodes(1)
  q.transactions.times do |t|
   assert_equal nodes[t], q.addedNode(t)
  end
  assert_equal EMPTY, q.addedNode(-1)
  assert_equal EMPTY, q.addedNode(5)
 end

 def XtestResetZerosOutPreviouslyAddedNodes
 end

 def XtestReallocMemoryWhenRunOutOfTransactions
 end

end
