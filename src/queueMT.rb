#!/usr/bin/env ruby
#
# Mobility test for queue c lib
#


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

 def testTransactionInitializedWithZeroRemovedEdges
  q = Queue.new(3)
  assert_equal EMPTY, q.removedEdges(-1)
  assert_equal 0,     q.removedEdges(0)
  assert_equal EMPTY, q.removedEdges(1)
  q.newTransaction
  assert_equal 0,     q.removedEdges(0)
  assert_equal 0,     q.removedEdges(1)
 end

 def testTransactionInitializedWithZeroAddedEdges
  q = Queue.new(3)
  assert_equal EMPTY, q.addedEdges(-1)
  assert_equal 0,     q.addedEdges(0)
  assert_equal EMPTY, q.addedEdges(1)
  q.newTransaction
  assert_equal 0,     q.addedEdges(0)
  assert_equal 0,     q.addedEdges(1)
 end

 def testRemoveCellTransaction
  nodes = [0,1,2,3]
  nodeParts = [10,11,12,13]
  q = Queue.new(3)
  assert_equal q,         q.removeCell(nodes,nodeParts)
  assert_equal 1,         q.removedCells(0)
  assert_equal nodes,     q.removedCellNodes(0)
  assert_equal nodeParts, q.removedCellNodeParts(0)
 end

 def testAddCellTransaction
  nodes = [0,1,2,3]
  nodeParts = [5,5,5,5]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  assert_equal q,         q.addCell(nodes,nodeParts,xyzs)
  assert_equal 1,         q.addedCells(0)
  assert_equal nodes,     q.addedCellNodes(0)
  assert_equal nodeParts, q.addedCellNodeParts(0)
  assert_equal xyzs,      q.addedCellXYZs(0)
 end

 def testAddCellTransactionWithBigNodeSize
  nodes = [0,1,2,3]
  nodeParts = [5,5,5,5]
  xyzs = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
          10,11,12,13,14,15,16,17,18,19,
          20,21,22,23,24,25,26,27,28,29,
          30,31,32,33,34,35,36,37,38,39]
  q = Queue.new(10)
  assert_equal q,     q.addCell(nodes,nodeParts,xyzs)
  assert_equal xyzs,  q.addedCellXYZs(0)
 end

 def testRemoveFaceTransaction
  nodes = [0,1,2]
  nodeParts = [10,11,12]
  q = Queue.new(3)
  assert_equal q,         q.removeFace(nodes,nodeParts)
  assert_equal 1,         q.removedFaces(0)
  assert_equal nodes,     q.removedFaceNodes(0)
  assert_equal nodeParts, q.removedFaceNodeParts(0)
 end

 def testAddFaceTransaction
  faceId = 9
  nodes = [0,1,2]
  nodeParts = [40,41,42]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3)
  assert_equal q,         q.addFace(nodes,faceId,nodeParts,uvs)
  assert_equal 1,         q.addedFaces(0)
  assert_equal nodes,     q.addedFaceNodes(0)
  assert_equal faceId,    q.addedFaceId(0)
  assert_equal nodeParts, q.addedFaceNodeParts(0)
  assert_equal uvs,       q.addedFaceUVs(0)
 end

 def testRemoveEdgeTransaction
  nodes = [0,1]
  nodeParts = [10,11]
  q = Queue.new(3)
  assert_equal q,         q.removeEdge(nodes,nodeParts)
  assert_equal 1,         q.removedEdges(0)
  assert_equal nodes,     q.removedEdgeNodes(0)
  assert_equal nodeParts, q.removedEdgeNodeParts(0)
 end

 def testAddEdgeTransaction
  edgeId = 9
  nodes = [0,1]
  nodeParts = [40,41]
  t = [ 10,11 ]
  q = Queue.new(3)
  assert_equal q,         q.addEdge(nodes,edgeId,nodeParts,t)
  assert_equal 1,         q.addedEdges(0)
  assert_equal nodes,     q.addedEdgeNodes(0)
  assert_equal edgeId,    q.addedEdgeId(0)
  assert_equal nodeParts, q.addedEdgeNodeParts(0)
  assert_equal t,         q.addedEdgeTs(0)
 end

 def testResetZerosOutPreviouslyRemovedCellNodes
  nodes = [0,1,2,3]
  nodeParts = [10,11,12,13]
  q = Queue.new(3).removeCell([9,9,9,9],[8,8,8,8]).reset
  assert_equal 0, q.removedCells(0)
  q.removeCell(nodes,nodeParts)
  assert_equal nodes,     q.removedCellNodes(0)  
  assert_equal nodeParts, q.removedCellNodeParts(0)  
 end

 def testResetZerosOutPreviouslyRemovedAndAddedCellNodes
  nodes = [0,1,2,3]
  nodeParts = [6,7,8,9]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3).addCell([1,2,3,4],[6,7,8,9],[1,2,3, 1,2,3, 1,2,3, 1,2,3]).reset
  assert_equal 0,         q.addedCells(0)
  q.addCell(nodes, nodeParts, xyzs)
  assert_equal nodes,     q.addedCellNodes(0)  
  assert_equal nodeParts, q.addedCellNodeParts(0)  
  assert_equal xyzs,      q.addedCellXYZs(0)  
 end

 def testResetZerosOutPreviouslyRemovedFaceNodes
  nodes = [0,1,2]
  nodeParts = [20,21,22]
  q = Queue.new(3).removeFace([9,9,9],[8,8,8]).reset
  assert_equal 0, q.removedFaces(0)
  q.removeFace(nodes,nodeParts)
  assert_equal nodes,     q.removedFaceNodes(0)  
  assert_equal nodeParts, q.removedFaceNodeParts(0)  
 end

 def testResetZerosOutPreviouslyRemovedAndAddedFaceNodes
  faceId = 8
  nodes = [0,1,2]
  nodeParts = [50,51,52]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3).addFace([9,9,9],60,[7,7,7],[1,2, 1,2, 1,2]).reset
  assert_equal 0, q.addedFaces(0)
  q.addFace(nodes,faceId,nodeParts,uvs)
  assert_equal nodes,     q.addedFaceNodes(0)
  assert_equal faceId,    q.addedFaceId(0)
  assert_equal nodeParts, q.addedFaceNodeParts(0)
  assert_equal uvs,       q.addedFaceUVs(0)
  end

 def testResetZerosOutPreviouslyRemovedEdgeNodes
  nodes = [0,1]
  nodeParts = [20,21]
  q = Queue.new(3).removeEdge([9,9],[8,8]).reset
  assert_equal 0, q.removedEdges(0)
  q.removeEdge(nodes,nodeParts)
  assert_equal nodes,     q.removedEdgeNodes(0)  
  assert_equal nodeParts, q.removedEdgeNodeParts(0)  
 end

 def testResetZerosOutPreviouslyRemovedAndAddedEdgeNodes
  edgeId = 8
  nodes = [0,1]
  nodeParts = [50,51]
  t = [ 30,31 ]
  q = Queue.new(3).addEdge([9,9],60,[7,7],[1,2]).reset
  assert_equal 0, q.addedEdges(0)
  q.addEdge(nodes,edgeId,nodeParts,t)
  assert_equal nodes,     q.addedEdgeNodes(0)
  assert_equal edgeId,    q.addedEdgeId(0)
  assert_equal nodeParts, q.addedEdgeNodeParts(0)
  assert_equal t,         q.addedEdgeTs(0)
 end

 def testResetCurrentTransactionFaces
  faceId = 8
  nodes = [0,1,2]
  nodeParts = [50,51,52]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3).addFace(nodes,faceId,nodeParts,uvs)
  q.newTransaction
  q.addFace([9,9,9],60,[7,7,7],[1,2, 1,2, 1,2])
  q.resetCurrentTransaction
  assert_equal 2, q.transactions
  assert_equal 1, q.addedFaces(0)
  assert_equal 0, q.addedFaces(1)
 end

 def testReallocMemoryForLotsOfTransactions
  q = Queue.new(3)
  10000.times { q.newTransaction }
 end

 def testReallocMemoryForLotsOfCellAdds
  nodes = [0,1,2,3]
  nodeParts = [5,6,7,8]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  10000.times { q.addCell(nodes, nodeParts, xyzs) }
 end

 def testReallocMemoryForLotsOfCellRemoves
  nodes = [0,1,2,3]
  nodeParts = [10,11,12,13]
  q = Queue.new(3)
  10000.times { q.removeCell(nodes,nodeParts) }
 end

 def testReallocMemoryForLotsOfFaceAdds
  faceId = 7
  nodes = [0,1,2]
  nodeParts = [50,51,52]
  uvs = [ 10,11, 20,21, 30,31 ]
  q = Queue.new(3)
  10000.times { q.addFace(nodes,faceId,nodeParts,uvs) }
 end

 def testReallocMemoryForLotsOfFaceRemoves
  nodes = [0,1,2]
  nodeParts = [50,51,52]
  q = Queue.new(3)
  10000.times { q.removeFace(nodes,nodeParts) }
 end

 def testReallocMemoryForLotsOfEdgeAdds
  edgeId = 7
  nodes = [0,1]
  nodeParts = [51,52]
  t = [ 20,21 ]
  q = Queue.new(3)
  10000.times { q.addEdge(nodes,edgeId,nodeParts,t) }
 end

 def testReallocMemoryForLotsOfEdgeRemoves
  nodes = [0,1]
  nodeParts = [50,51]
  q = Queue.new(3)
  10000.times { q.removeEdge(nodes,nodeParts) }
 end

 def testReallocMemoryForLotsOfTransactionsPlus
  nodes = [0,1,2,3]
  id = 4
  nodeParts = [5,6,7,8]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  10000.times do
   q.newTransaction
   q.addCell(nodes,nodeParts,xyzs).removeCell(nodes,nodeParts)
   q.addFace(nodes,id,nodeParts,xyzs).removeFace(nodes,nodeParts)
   q.addEdge(nodes,id,nodeParts,xyzs).removeEdge(nodes,nodeParts)
  end
 end

 def testExposeTotalNumberOfRemovedCells
  nodes = [0,1,2,3]
  nodeParts = [10,11,12,13]
  q = Queue.new(3)
  assert_equal 0, q.totalRemovedCells
  q.removeCell(nodes,nodeParts)
  assert_equal 1, q.totalRemovedCells
  q.removeCell(nodes,nodeParts)
  assert_equal 2, q.totalRemovedCells
 end

 def testSerializeTransactions
  q = Queue.new(3)
  assert_equal [3,1,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0], q.dump
  q.newTransaction.newTransaction
  assert_equal [3,3,0,0,0,0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0], q.dump
 end

 def testSerializeRemoveCells
  nodes = [5,6,7,8]
  nodeParts = [15,16,17,18]
  q = Queue.new(3).removeCell(nodes,nodeParts)
  assert_equal [3,1,1,0,0,0,0,0, 1, 5,6,7,8,  15,16,17,18, 0, 0, 0, 0, 0 ], q.dump
 end

 def testSerializeAddCells3
  nodes = [5,6,7,8]
  nodeParts = [10,11,12,13]
  xyz = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
  q = Queue.new(3).addCell(nodes,nodeParts,xyz)
  assert_equal [3,1,0,1,0,0,0,0, 0, 1, 5,6,7,8, 10,11,12,13, 0, 0, 0, 0, 
    0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2  ], q.dump
 end

 def testSerializeAddCells1
  nodes = [5,6,7,8]
  nodeParts = [10,11,12,13]
  xyz = [0.1,0.2,0.3,0.4]
  q = Queue.new(1).addCell(nodes,nodeParts,xyz)
  assert_equal [1,1,0,1,0,0,0,0, 0, 1, 5,6,7,8, 10,11,12,13, 0, 0, 0, 0, 
    0.1,0.2,0.3,0.4 ], q.dump
 end

 def testSerializeRemoveFaces
  nodes = [5,6,7]
  nodeParts = [15,16,17]
  q = Queue.new(3).removeFace(nodes,nodeParts)
  assert_equal [3,1,0,0,1,0,0,0, 0, 0, 1, 5,6,7, 15,16,17, 0, 0, 0 ], q.dump
 end

 def testSerializeAddFaces
  faceId = 99
  nodes = [5,6,7]
  nodeParts = [25,26,27]
  uv = [0.1,0.2,0.3,0.4,0.5,0.6]
  q = Queue.new(3).addFace(nodes,faceId,nodeParts,uv)
  assert_equal [3,1,0,0,0,1,0,0, 0, 0, 0, 1, 5,6,7, 99, 25,26,27, 0, 0,
    0.1,0.2,0.3,0.4,0.5,0.6 ], q.dump
 end

 def testSerializeRemoveEdges
  nodes = [5,6]
  nodeParts = [15,16]
  q = Queue.new(3).removeEdge(nodes,nodeParts)
  assert_equal [3,1,0,0,0,0,1,0, 0, 0, 0, 0, 1, 5,6, 15,16, 0 ], q.dump
 end

 def testSerializeAddEdges
  edgeId = 99
  nodes = [5,6]
  nodeParts = [25,26]
  uv = [0.1,0.2]
  q = Queue.new(3).addEdge(nodes,edgeId,nodeParts,uv)
  assert_equal [3,1,0,0,0,0,0,1, 0, 0, 0, 0, 0, 1, 5,6, 99, 25,26, 
    0.1,0.2 ], q.dump
 end

 def testLoadSerializedQueuewithZeroTransactions
  q = Queue.new(3).newTransaction.newTransaction
  assert_equal 3, q.transactions
  assert_equal q, q.load([3,1,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0], [])
  assert_equal 1, q.transactions
 end

 def testLoadFailsWhenNodeSizeChanges
  q = Queue.new(3)
  assert_nil q.load([5,1,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0], [])
 end

 def testLoadSerializedQueuewithThreeTransactions
  q = Queue.new(3)
  assert_equal 1, q.transactions
  assert_equal q, q.load([3,3,0,0,0,0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0], [])
  assert_equal 3, q.transactions
 end

 def testLoadSerializedQueuewithRemovedCell
  q = Queue.new(3)
  assert_equal q, q.load([3,1,1,0,0,0,0,0, 1, 5,6,7,8, 15,16,17,18, 0, 0, 0, 0, 0 ], [])
  assert_equal 1, q.removedCells(0)
  assert_equal [5,6,7,8], q.removedCellNodes(0)
  assert_equal [15,16,17,18], q.removedCellNodeParts(0)
 end

 def testLoadSerializedQueuewithAddedCell
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,1,0,0,0,0, 0, 1, 5,6,7,8, 10,11,12,13, 0, 0, 0, 0 ], 
			 [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2])
  assert_equal 1, q.addedCells(0)
  assert_equal [5,6,7,8], q.addedCellNodes(0)
  assert_equal [10,11,12,13], q.addedCellNodeParts(0)
  assert_equal [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2], q.addedCellXYZs(0)
 end

 def testLoadSerializedQueuewithRemovedFace
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,0,1,0,0,0, 0, 0, 1, 5,6,7, 35,36,37, 0, 0, 0 ], [])
  assert_equal 1, q.removedFaces(0)
  assert_equal [5,6,7], q.removedFaceNodes(0)
  assert_equal [35,36,37], q.removedFaceNodeParts(0)
 end

 def testLoadSerializedQueuewithAddedFace
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,0,0,1,0,0, 0, 0, 0, 1, 5,6,7, 88, 45,46,47, 0, 0 ], 
			 [0.1,0.2,0.3,0.4,0.5,0.6])
  assert_equal 1, q.addedFaces(0)
  assert_equal [5,6,7], q.addedFaceNodes(0)
  assert_equal 88,      q.addedFaceId(0)
  assert_equal [45,46,47], q.addedFaceNodeParts(0)
  assert_equal [0.1,0.2,0.3,0.4,0.5,0.6], q.addedFaceUVs(0)
 end

 def testLoadSerializedQueuewithRemovedEdge
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,0,0,0,1,0, 0, 0, 0, 0, 1, 5,6, 35,36, 0 ], [])
  assert_equal 1, q.removedEdges(0)
  assert_equal [5,6], q.removedEdgeNodes(0)
  assert_equal [35,36], q.removedEdgeNodeParts(0)
 end

 def testLoadSerializedQueuewithAddedEdge
  q = Queue.new(3)
  assert_equal q, q.load([3,1,0,0,0,0,0,1, 0, 0, 0, 0, 0, 1, 5,6, 88, 45,46 ], 
			 [1.1,1.2])
  assert_equal 1, q.addedEdges(0)
  assert_equal [5,6], q.addedEdgeNodes(0)
  assert_equal 88,      q.addedEdgeId(0)
  assert_equal [45,46], q.addedEdgeNodeParts(0)
  assert_equal [1.1,1.2], q.addedEdgeTs(0)
 end

 def testSerialIsStableAfterARemalloc
  nodes = [5,6,7,8]
  id = 9
  nodeParts = [10,11,12,13]
  nodes = [5,6,7,8,9,10,11,12,13]
  xyz = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
  q = Queue.new(3)
  n=1000
  n.times { 
   q.newTransaction
   q.removeCell(nodes,nodeParts).addCell(nodes,nodeParts,xyz)
   q.removeCell(nodes,nodeParts).addCell(nodes,nodeParts,xyz)
   q.removeFace(nodes,nodeParts).addFace(nodes,id,nodeParts,xyz)
   q.removeFace(nodes,nodeParts).addFace(nodes,id,nodeParts,xyz)
   q.removeEdge(nodes,nodeParts).addEdge(nodes,id,nodeParts,xyz)
   q.removeEdge(nodes,nodeParts).addEdge(nodes,id,nodeParts,xyz)
  } 
  i = q.dumpInt
  f = q.dumpFloat
  l = Queue.new(3).load(i,f)
  assert_equal n+1, l.transactions
  assert_equal i, l.dumpInt
  assert_equal f, l.dumpFloat
 end

 def testAddRemoveCellAndShift
  nodes = [0,1,2,3]
  nodeParts = [5,5,5,5]
  xyzs = [ 0, 1, 2, 10,11,12, 20,21,22, 30,31,32 ]
  q = Queue.new(3)
  q.addCell(nodes,nodeParts,xyzs)
  q.addCell(nodes,nodeParts,xyzs)
  q.removeCell(nodes,nodeParts)

  assert_equal q, q.globalShiftNode(3,2)
  shiftedNodes = [0,1,2,5]
  assert_equal shiftedNodes, q.addedCellNodes(0)
  assert_equal shiftedNodes, q.addedCellNodes(1)
  assert_equal shiftedNodes, q.removedCellNodes(0)
 end

end
