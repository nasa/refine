#!/usr/bin/env ruby
#
# Mobility test for queue c lib
#
# $Id$

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Queue').build

require 'test/unit'
require 'Queue/Queue'

class TestQueue < Test::Unit::TestCase

 EMPTY = -1

 def testInitializedWithOneTransaction
  assert_equal 1, Queue.new(3).transactions
 end

 def testInitializedWithNodeSize
  assert_equal 3, Queue.new(3).nodeSize
  assert_equal 15, Queue.new(15).nodeSize
  assert_equal 26, Queue.new(26).nodeSize
 end

 def testNewTransactionIncrementsTransactions
  q = Queue.new(3)
  assert_equal q, q.newTransaction
  assert_equal 2, q.transactions
  assert_equal q, q.newTransaction
  assert_equal 3, q.transactions
 end

 def testResetSetsTheNumberOfTransactionsToOne
  q = Queue.new(3).newTransaction.newTransaction
  assert_equal q, q.reset
  assert_equal 1, q.transactions
 end

 def testTransactionInitializedWithZeroRemovedCells
  q = Queue.new(3)
  assert_equal EMPTY, q.removedCells(-1)
  assert_equal 0,     q.removedCells(0)
  assert_equal EMPTY, q.removedCells(1)
  q.newTransaction
  assert_equal 0,     q.removedCells(0)
  assert_equal 0,     q.removedCells(1)
 end

 def testTransactionInitializedWithZeroAddedCells
  q = Queue.new(3)
  assert_equal EMPTY, q.addedCells(-1)
  assert_equal 0,     q.addedCells(0)
  assert_equal EMPTY, q.addedCells(1)
  q.newTransaction
  assert_equal 0,     q.addedCells(0)
  assert_equal 0,     q.addedCells(1)
 end

 def testTransactionInitializedWithZeroRemovedFaces
  q = Queue.new(3)
  assert_equal EMPTY, q.removedFaces(-1)
  assert_equal 0,     q.removedFaces(0)
  assert_equal EMPTY, q.removedFaces(1)
  q.newTransaction
  assert_equal 0,     q.removedFaces(0)
  assert_equal 0,     q.removedFaces(1)
 end

 def testTransactionInitializedWithZeroAddedFaces
  q = Queue.new(3)
  assert_equal EMPTY, q.addedFaces(-1)
  assert_equal 0,     q.addedFaces(0)
  assert_equal EMPTY, q.addedFaces(1)
  q.newTransaction
  assert_equal 0,     q.addedFaces(0)
  assert_equal 0,     q.addedFaces(1)
 end

 def testRemoveCellTransaction
  nodes = [0,1,2,3]
  q = Queue.new(3)
  assert_equal q,     q.removeCell(nodes)
  assert_equal 1,     q.removedCells(0)
  assert_equal nodes, q.removedCellNodes(0)
 end

 def testAddCellTransaction
  cellId = 7
  nodes = [0,1,2,3,cellId,5,5,5,5]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  assert_equal q,     q.addCell(nodes,xyzs)
  assert_equal 1,     q.addedCells(0)
  assert_equal nodes, q.addedCellNodes(0)
  assert_equal xyzs,  q.addedCellXYZs(0)
 end

 def testAddCellTransactionWithBigNodeSize
  cellId = 7
  nodes = [0,1,2,3,cellId,5,5,5,5]
  xyzs = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
          10,11,12,13,14,15,16,17,18,19,
          20,21,22,23,24,25,26,27,28,29,
          30,31,32,33,34,35,36,37,38,39]
  q = Queue.new(10)
  assert_equal q,     q.addCell(nodes,xyzs)
  assert_equal 1,     q.addedCells(0)
  assert_equal nodes, q.addedCellNodes(0)
  assert_equal xyzs,  q.addedCellXYZs(0)
 end

 def testRemoveFaceTransaction
  nodes = [0,1,2]
  q = Queue.new(3)
  assert_equal q,     q.removeFace(nodes)
  assert_equal 1,     q.removedFaces(0)
  assert_equal nodes, q.removedFaceNodes(0)
 end

 def testAddFaceTransaction
  faceId = 9
  nodes = [0,1,2,faceId]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3)
  assert_equal q,     q.addFace(nodes,uvs)
  assert_equal 1,     q.addedFaces(0)
  assert_equal nodes, q.addedFaceNodes(0)
  assert_equal uvs,   q.addedFaceUVs(0)
 end

 def testResetZerosOutPreviouslyRemovedCellNodes
  nodes = [0,1,2,3]
  q = Queue.new(3).removeCell([9,9,9,9]).reset
  assert_equal 0, q.removedCells(0)
  q.removeCell(nodes)
  assert_equal nodes, q.removedCellNodes(0)  
 end

 def testResetZerosOutPreviouslyRemovedAndAddedCellNodes
  nodesCellId = [0,1,2,3,5,6,7,8,9]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3).addCell([1,2,3,4,5,6,7,8,9],[1,2,3, 1,2,3, 1,2,3, 1,2,3]).reset
  assert_equal 0, q.addedCells(0)
  q.addCell(nodesCellId,xyzs)
  assert_equal nodesCellId, q.addedCellNodes(0)  
  assert_equal xyzs, q.addedCellXYZs(0)  
 end

 def testResetZerosOutPreviouslyRemovedFaceNodes
  nodes = [0,1,2]
  q = Queue.new(3).removeFace([9,9,9]).reset
  assert_equal 0, q.removedFaces(0)
  q.removeFace(nodes)
  assert_equal nodes, q.removedFaceNodes(0)  
 end

 def testResetZerosOutPreviouslyRemovedAndAddedFaceNodes
  nodesFaceId = [0,1,2,8]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3).addFace([9,9,9,9],[1,2, 1,2, 1,2]).reset
  assert_equal 0, q.addedFaces(0)
  q.addFace(nodesFaceId,uvs)
  assert_equal nodesFaceId, q.addedFaceNodes(0)  
  assert_equal uvs, q.addedFaceUVs(0)  
 end

 def testReallocMemoryForLotsOfTransactions
  q = Queue.new(3)
  10000.times { q.newTransaction }
 end

 def testReallocMemoryForLotsOfCellAdds
  nodes = [0,1,2,3,4,5,6,7,8]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  10000.times { q.addCell(nodes,xyzs) }
 end

 def testReallocMemoryForLotsOfCellRemoves
  nodes = [0,1,2,3]
  q = Queue.new(3)
  10000.times { q.removeCell(nodes) }
 end

 def testReallocMemoryForLotsOfFaceAdds
  nodes = [0,1,2,3]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3)
  10000.times { q.addFace(nodes,uvs) }
 end

 def testReallocMemoryForLotsOfFaceRemoves
  nodes = [0,1,2]
  q = Queue.new(3)
  10000.times { q.removeFace(nodes) }
 end

 def testReallocMemoryForLotsOfTransactionsPlus
  nodes = [0,1,2,3,4,5,6,7,8]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  10000.times do
   q.newTransaction
   q.addCell(nodes,xyzs).removeCell(nodes)
   q.addFace(nodes,xyzs).removeFace(nodes)
  end
 end

 def testExposeTotalNumberOfRemovedCells
  nodes = [0,1,2,3]
  q = Queue.new(3)
  assert_equal 0, q.totalRemovedCells
  q.removeCell(nodes)
  assert_equal 1, q.totalRemovedCells
  q.removeCell(nodes)
  assert_equal 2, q.totalRemovedCells
 end

 def testSerializeTransactions
  q = Queue.new(3)
  assert_equal [3,1,0,0,0,0, 0, 0, 0, 0], q.dump
  q.newTransaction.newTransaction
  assert_equal [3,3,0,0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0], q.dump
 end

 def testSerializeRemoveCells
  nodes = [5,6,7,8]
  q = Queue.new(3).removeCell(nodes)
  assert_equal [3,1,1,0,0,0, 1, 5,6,7,8, 0, 0, 0 ], q.dump
 end

 def testSerializeAddCells3
  nodes = [5,6,7,8,9,10,11,12,13]
  xyz = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
  q = Queue.new(3).addCell(nodes,xyz)
  assert_equal [3,1,0,1,0,0, 0, 1, 5,6,7,8,9,10,11,12,13, 0, 0, 
    0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2  ], q.dump
 end

 def testSerializeAddCells1
  nodes = [5,6,7,8,9,10,11,12,13]
  xyz = [0.1,0.2,0.3,0.4]
  q = Queue.new(1).addCell(nodes,xyz)
  assert_equal [1,1,0,1,0,0, 0, 1, 5,6,7,8,9,10,11,12,13, 0, 0, 
    0.1,0.2,0.3,0.4 ], q.dump
 end

 def testSerializeRemoveFaces
  nodes = [5,6,7]
  q = Queue.new(3).removeFace(nodes)
  assert_equal [3,1,0,0,1,0, 0, 0, 1, 5,6,7, 0 ], q.dump
 end

 def testSerializeAddFaces
  nodes = [5,6,7,8]
  uv = [0.1,0.2,0.3,0.4,0.5,0.6]
  q = Queue.new(3).addFace(nodes,uv)
  assert_equal [3,1,0,0,0,1, 0, 0, 0, 1, 5,6,7,8,
    0.1,0.2,0.3,0.4,0.5,0.6 ], q.dump
 end

 def testLoadSerializedQueuewithZeroTransactions
  q = Queue.new(3).newTransaction.newTransaction
  assert_equal 3, q.transactions
  assert_equal q, q.load([3,1,0,0,0,0, 0, 0, 0, 0], [])
  assert_equal 1, q.transactions
 end

 def testLoadFailsWhenNodeSizeChanges
  q = Queue.new(3)
  assert_nil q.load([5,1,0,0,0,0, 0, 0, 0, 0], [])
 end

 def testLoadSerializedQueuewithThreeTransactions
  q = Queue.new(3)
  assert_equal 1, q.transactions
  assert_equal q, q.load([3,3,0,0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0], [])
  assert_equal 3, q.transactions
 end

 def testLoadSerializedQueuewithRemovedCell
  q = Queue.new(3)
  assert_equal q, q.load([3,1,1,0,0,0, 1, 5,6,7,8, 0, 0, 0 ], [])
  assert_equal 1, q.removedCells(0)
  assert_equal [5,6,7,8], q.removedCellNodes(0)
 end

 def testLoadSerializedQueuewithAddedCell
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,1,0,0, 0, 1, 5,6,7,8,9,10,11,12,13, 0, 0 ], 
			 [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2])
  assert_equal 1, q.addedCells(0)
  assert_equal [5,6,7,8,9,10,11,12,13], q.addedCellNodes(0)
  assert_equal [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2], q.addedCellXYZs(0)
 end

 def testLoadSerializedQueuewithRemovedFace
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,0,1,0, 0, 0, 1, 5,6,7, 0 ], [])
  assert_equal 1, q.removedFaces(0)
  assert_equal [5,6,7], q.removedFaceNodes(0)
 end

 def testLoadSerializedQueuewithAddedFace
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,0,0,1, 0, 0, 0, 1, 5,6,7,8 ], 
			 [0.1,0.2,0.3,0.4,0.5,0.6])
  assert_equal 1, q.addedFaces(0)
  assert_equal [5,6,7,8], q.addedFaceNodes(0)
  assert_equal [0.1,0.2,0.3,0.4,0.5,0.6], q.addedFaceUVs(0)
 end

 def testSerialIsStatbleAfterARemalloc
  nodes = [5,6,7,8,9,10,11,12,13]
  xyz = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
  q = Queue.new(3)
  n=1000
  n.times { 
   q.newTransaction
   q.removeCell(nodes).addCell(nodes,xyz)
   q.removeCell(nodes).addCell(nodes,xyz)
   q.removeFace(nodes).addFace(nodes,xyz)
   q.removeFace(nodes).addFace(nodes,xyz)
  } 
  i = q.dumpInt
  f = q.dumpFloat
  l = Queue.new(3).load(i,f)
  assert_equal n+1, l.transactions
  assert_equal i, l.dumpInt
  assert_equal f, l.dumpFloat
 end

end
