#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

exit 1 unless system 'ruby makeRubyExtension.rb Grid adj.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridMetric adj.c grid.c gridStruct.h master_header.h'
exit 1 unless system 'ruby makeRubyExtension.rb GridSwap adj.c grid.c gridmetric.h gridStruct.h master_header.h'

require 'test/unit'
require 'Grid/Grid'
require 'GridSwap/GridSwap'
require 'GridMetric/GridMetric'

class Grid
 include GridMetric
 include GridSwap
 def totalVolume
  vol = 0.0
  ncell.times { |cellId| vol += volume(cell(cellId)) }
  vol
 end
end

class TestGridSwap < Test::Unit::TestCase

 def testSwapEdgeNotBeter
  assert_not_nil grid=gemGrid(3, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
  assert_not_nil grid=gemGrid(4, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
  assert_not_nil grid=gemGrid(5, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
  assert_not_nil grid=gemGrid(6, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
  assert_not_nil grid=gemGrid(7, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
  assert_not_nil grid=gemGrid(8, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
  assert_not_nil grid=gemGrid(9, 2.0, 0)
  assert_nil grid.swapEdge(0,1)
 end

 def testSwapEdge3
  assert_not_nil grid=gemGrid(3, 0.1, 0)
  initalVolume = grid.totalVolume
  assert_equal grid, grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 1, grid.cellDegree(0)
  assert_equal 1, grid.cellDegree(1)
  assert_equal 2, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 2, grid.cellDegree(4)
 end

 def testSwapEdge4_0
  assert_not_nil grid=gemGrid(4, 0.1, 0)
  initalVolume = grid.totalVolume
  assert_equal grid, grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 2, grid.cellDegree(0)
  assert_equal 2, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 2, grid.cellDegree(5)
 end

 def testSwapEdge4_1
  assert_not_nil grid=gemGrid(4, 0.1, 1)
  initalVolume = grid.totalVolume
  assert_equal grid, grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 2, grid.cellDegree(0)
  assert_equal 2, grid.cellDegree(1)
  assert_equal 2, grid.cellDegree(2)
  assert_equal 4, grid.cellDegree(3)
  assert_equal 2, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
 end

 def testSwapEdge4_invalid
  assert_not_nil grid=gemGrid(4, 0.1, nil, -0.1)
  grid.swapEdge(0,1)
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
 end


 def testSwapEdge4_gapWithSameFaces
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal 3,    grid.ncell
  assert_equal grid, grid.addFaceUV(0,0.0,10.0,
				    1,1.0,11.0,
				    2,2.0,12.0,11)
  assert_equal grid, grid.addFaceUV(1,1.0,11.0,
				    0,0.0,10.0,
				    5,5.0,15.0,11)
  assert grid.rightHandedBoundary, "original boundary is not right handed"
  assert_equal grid, grid.swapEdge(0,1)
  assert_equal 4,    grid.ncell
  assert_equal 2,    grid.cellDegree(0)
  assert_equal 2,    grid.cellDegree(1)
  # new bc faces
  assert_equal 11,   grid.faceId(0,2,5)
  assert_equal 11,   grid.faceId(1,2,5)  
  assert_nil         grid.faceId(0,1,2)
  assert_nil         grid.faceId(0,1,5) 
  assert grid.rightHandedBoundary, "swapped boundary is not right handed"
  assert_equal [0.0,10.0], grid.nodeUV(0,11)
  assert_equal [1.0,11.0], grid.nodeUV(1,11)
  assert_equal [2.0,12.0], grid.nodeUV(2,11)
  assert_equal [5.0,15.0], grid.nodeUV(5,11)
 end

 def testSwapEdge4_gapWithSameFaceNotBeter
  assert_not_nil     grid=gemGrid(4, 2.0, nil, nil, true)
  assert_equal 3,    grid.ncell
  assert_equal grid, grid.addFace(0,1,2,11)
  assert_equal grid, grid.addFace(0,1,5,11)
  assert_nil         grid.swapEdge(0,1)
  assert_equal 3,    grid.ncell
  assert_equal 3,    grid.cellDegree(0)
  assert_equal 3,    grid.cellDegree(1)
  assert_equal 11,   grid.faceId(0,1,2)
  assert_equal 11,   grid.faceId(0,1,5)  
 end

 def testSwapEdge4_gapWithDifferentFaces
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal 3,    grid.ncell
  assert_equal grid, grid.addFace(0,1,2,2)
  assert_equal grid, grid.addFace(0,1,5,5)
  assert_nil         grid.swapEdge(0,1)
  assert_equal 3,    grid.ncell
  assert_equal 3,    grid.cellDegree(0)
  assert_equal 3,    grid.cellDegree(1)
  assert_equal 2,    grid.faceId(0,1,2)
  assert_equal 5,    grid.faceId(0,1,5)  
 end

 def testSwapEdge4_gapWithSameAndExistingFace0
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal 3,    grid.ncell
  assert_equal grid, grid.addFace(0,1,2,11)
  assert_equal grid, grid.addFace(0,1,5,11)
  assert_equal grid, grid.addFace(0,2,5,20)
  assert_nil         grid.swapEdge(0,1)
  assert_equal 3,    grid.ncell
  assert_equal 3,    grid.cellDegree(0)
  assert_equal 3,    grid.cellDegree(1)
  assert_equal 11,   grid.faceId(0,1,2)
  assert_equal 11,   grid.faceId(0,1,5)  
 end

 def testSwapEdge4_gapWithSameAndExistingFace1
  assert_not_nil     grid=gemGrid(4, nil, nil, nil, true)
  assert_equal 3,    grid.ncell
  assert_equal grid, grid.addFace(0,1,2,11)
  assert_equal grid, grid.addFace(0,1,5,11)
  assert_equal grid, grid.addFace(1,2,5,20)
  assert_nil         grid.swapEdge(0,1)
  assert_equal 3,    grid.ncell
  assert_equal 3,    grid.cellDegree(0)
  assert_equal 3,    grid.cellDegree(1)
  assert_equal 11,   grid.faceId(0,1,2)
  assert_equal 11,   grid.faceId(0,1,5)  
 end

 def testSwapEdge5_0
  assert_not_nil grid=gemGrid(5, 0.1, 0)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 3, grid.cellDegree(0)
  assert_equal 3, grid.cellDegree(1)
  assert_equal 6, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
  assert_equal 2, grid.cellDegree(6)
 end

 def testSwapEdge5_1
  assert_not_nil grid=gemGrid(5, 0.1, 1)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 3, grid.cellDegree(0)
  assert_equal 3, grid.cellDegree(1)
  assert_equal 2, grid.cellDegree(2)
  assert_equal 6, grid.cellDegree(3)
  assert_equal 2, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
  assert_equal 4, grid.cellDegree(6)
 end

 def testSwapEdge5_2
  assert_not_nil grid=gemGrid(5, 0.1, 2)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 3, grid.cellDegree(0)
  assert_equal 3, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 6, grid.cellDegree(4)
  assert_equal 2, grid.cellDegree(5)
  assert_equal 4, grid.cellDegree(6)
 end

 def testSwapEdge5_3
  assert_not_nil grid=gemGrid(5, 0.1, 3)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 3, grid.cellDegree(0)
  assert_equal 3, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 4, grid.cellDegree(3)
  assert_equal 2, grid.cellDegree(4)
  assert_equal 6, grid.cellDegree(5)
  assert_equal 2, grid.cellDegree(6)
 end

 def testSwapEdge5_4
  assert_not_nil grid=gemGrid(5, 0.1, 4)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 3, grid.cellDegree(0)
  assert_equal 3, grid.cellDegree(1)
  assert_equal 2, grid.cellDegree(2)
  assert_equal 4, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 2, grid.cellDegree(5)
  assert_equal 6, grid.cellDegree(6)
 end

 def testSwapEdge5_gapWithSameFaces
  assert_not_nil     grid=gemGrid(5, nil, nil, nil, true)
  assert_equal 4,    grid.ncell
  assert_equal grid, grid.addFace(0,1,2,11)
  assert_equal grid, grid.addFace(1,0,6,11)
  assert grid.rightHandedBoundary, "original boundary is not right handed"
  assert_equal grid, grid.swapEdge(0,1)
  assert_equal 6,    grid.ncell
  assert_equal 3,    grid.cellDegree(0)
  assert_equal 3,    grid.cellDegree(1)
  # new bc faces
  assert_equal 11,   grid.faceId(0,2,6)
  assert_equal 11,   grid.faceId(1,2,6)  
  assert_nil         grid.faceId(0,1,2)
  assert_nil         grid.faceId(0,1,6) 
  assert grid.rightHandedBoundary, "swapped boundary is not right handed"
 end

 def testSwapEdge6OnePoint_0
  assert_not_nil grid=gemGrid(6, 0.1, 0)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
  assert_equal 8, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
  assert_equal 4, grid.cellDegree(6)
  assert_equal 2, grid.cellDegree(7)
 end

 def testSwapEdge6OnePoint_1
  assert_not_nil grid=gemGrid(6, 0.1, 1)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
  assert_equal 2, grid.cellDegree(2)
  assert_equal 8, grid.cellDegree(3)
  assert_equal 2, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
  assert_equal 4, grid.cellDegree(6)
  assert_equal 4, grid.cellDegree(7)
 end

 def testSwapEdge6OnePoint_2
  assert_not_nil grid=gemGrid(6, 0.1, 2)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 8, grid.cellDegree(4)
  assert_equal 2, grid.cellDegree(5)
  assert_equal 4, grid.cellDegree(6)
  assert_equal 4, grid.cellDegree(7)
 end

 def testSwapEdge6OnePoint_3
  assert_not_nil grid=gemGrid(6, 0.1, 3)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 4, grid.cellDegree(3)
  assert_equal 2, grid.cellDegree(4)
  assert_equal 8, grid.cellDegree(5)
  assert_equal 2, grid.cellDegree(6)
  assert_equal 4, grid.cellDegree(7)
 end

 def testSwapEdge6OnePoint_4
  assert_not_nil grid=gemGrid(6, 0.1, 4)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 4, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 2, grid.cellDegree(5)
  assert_equal 8, grid.cellDegree(6)
  assert_equal 2, grid.cellDegree(7)
 end

 def testSwapEdge6OnePoint_5
  assert_not_nil grid=gemGrid(6, 0.1, 5)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 4, grid.cellDegree(0)
  assert_equal 4, grid.cellDegree(1)
  assert_equal 2, grid.cellDegree(2)
  assert_equal 4, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
  assert_equal 2, grid.cellDegree(6)
  assert_equal 8, grid.cellDegree(7)
 end

 def testSwapEdge7OnePoint_0
  assert_not_nil grid=gemGrid(7, 0.1, 0)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal  5, grid.cellDegree(0)
  assert_equal  5, grid.cellDegree(1)
  assert_equal 10, grid.cellDegree(2)
 end

 def testSwapEdge7OnePoint_6
  assert_not_nil grid=gemGrid(7, 0.1, 6)
  initalVolume = grid.totalVolume
  grid.swapEdge(0,1)
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal  5, grid.cellDegree(0)
  assert_equal  5, grid.cellDegree(1)
  assert_equal 10, grid.cellDegree(8)
 end

 def testSwap4
  assert_not_nil grid=gemGrid(4, 0.1, 0)
  initalVolume = grid.totalVolume
  assert_equal grid, grid.swap
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 2, grid.cellDegree(0)
  assert_equal 2, grid.cellDegree(1)
  assert_equal 4, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 2, grid.cellDegree(5)
 end

 def testSwap5
  assert_not_nil grid=gemGrid(5, 0.1, 0)
  initalVolume = grid.totalVolume
  grid.swap
  assert grid.minVolume>0.0, "negative volume cell "+grid.minVolume.to_s
  assert_in_delta initalVolume, grid.totalVolume, 1.0e-15
  assert_equal 3, grid.cellDegree(0)
  assert_equal 3, grid.cellDegree(1)
  assert_equal 6, grid.cellDegree(2)
  assert_equal 2, grid.cellDegree(3)
  assert_equal 4, grid.cellDegree(4)
  assert_equal 4, grid.cellDegree(5)
  assert_equal 2, grid.cellDegree(6)
 end

 def gemGrid(nequ=4, a=nil, dent=nil, x0 = nil, gap = nil)
  a  = a  || 0.1
  x0 = x0 || 1.0
  grid = Grid.new(nequ+2+1,14,14,14)
  n = Array.new
  n.push grid.addNode(  x0,0.0,0.0)
  n.push grid.addNode(-1.0,0.0,0.0)
  nequ.times do |i| 
   angle = 2.0*Math::PI*(i-1)/(nequ)
   s = if (dent==i) then 0.9 else 1.0 end
   n.push grid.addNode(0.0,s*a*Math.sin(angle),s*a*Math.cos(angle)) 
  end
  n.push 2
  ngem = (gap)?(nequ-1):(nequ)
  ngem.times do |i|
   grid.addCell(n[0],n[1],n[i+2],n[i+3])
  end
  grid  
 end

end
