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

 def testTransactionInitializedWithZeroRemovedCells
  q = Queue.new
  assert_equal EMPTY, q.removedCells(-1)
  assert_equal 0,     q.removedCells(0)
  assert_equal EMPTY, q.removedCells(1)
  q.newTransaction
  assert_equal 0,     q.removedCells(0)
  assert_equal 0,     q.removedCells(1)
 end

 def testTransactionInitializedWithZeroAddedCells
  q = Queue.new
  assert_equal EMPTY, q.addedCells(-1)
  assert_equal 0,     q.addedCells(0)
  assert_equal EMPTY, q.addedCells(1)
  q.newTransaction
  assert_equal 0,     q.addedCells(0)
  assert_equal 0,     q.addedCells(1)
 end

 def testRemoveTransaction
  nodes = [0,1,2,3]
  q = Queue.new
  assert_equal q,     q.removeCell(nodes)
  assert_equal 1,     q.removedCells(0)
  assert_equal nodes, q.removedCellNodes(0)
 end

 def testAddTransaction
  cellId = 7
  nodes = [0,1,2,3,cellId]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new
  assert_equal q,     q.addCell(nodes,xyzs)
  assert_equal 1,     q.addedCells(0)
  assert_equal nodes, q.addedCellNodes(0)
  assert_equal xyzs,  q.addedCellXYZs(0)
 end

 def testResetZerosOutPreviouslyRemovedCellNodes
  nodes = [0,1,2,3]
  q = Queue.new.removeCell([9,9,9,9]).reset
  assert_equal 0, q.removedCells(0)
  q.removeCell(nodes)
  assert_equal nodes, q.removedCellNodes(0)  
 end

 def testResetZerosOutPreviouslyRemovedAndAddedCellNodes
  nodesCellId = [0,1,2,3,7]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new.addCell([9,9,9,9,9],[1,2,3, 1,2,3, 1,2,3, 1,2,3]).reset
  assert_equal 0, q.addedCells(0)
  q.addCell(nodesCellId,xyzs)
  assert_equal nodesCellId, q.addedCellNodes(0)  
  assert_equal xyzs, q.addedCellXYZs(0)  
 end

 def testReallocMemoryForLotsOfTransactions
  q = Queue.new
  10000.times { q.newTransaction }
 end

 def testReallocMemoryForLotsOfCellAdds
  nodes = [0,1,2,3,4]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new
  10000.times { q.addCell(nodes,xyzs) }
 end

 def testReallocMemoryForLotsOfCellRemoves
  nodes = [0,1,2,3]
  q = Queue.new
  10000.times { q.removeCell(nodes) }
 end

end
