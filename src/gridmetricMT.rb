#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('GridMetric').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMetric/GridMetric'

class Grid
 include GridMetric
end

class TestGridMetric < Test::Unit::TestCase

 def rightTet
  grid = Grid.new(4,1,2,0)
  grid.addCell( 
	       grid.addNode(0.0,0.0,0.0), 
	       grid.addNode(1.0,0.0,0.0), 
	       grid.addNode(0.0,1.0,0.0), 
	       grid.addNode(0.0,0.0,1.0) )
  grid
 end

 def isoTet
  grid = Grid.new(4,1,0,0)
  grid.addNode( 0.000, 0.000, 0.000 )
  grid.addNode( 1.000, 0.000, 0000 )
  grid.addNode( 0.500, 0.866, 0.000 )
  grid.addNode( 0.500, 0.289, 0.823 ) 
  grid.addCell(0,1,2,3)
  grid
 end

 def testIsotropicTet
  assert_in_delta 1.000, isoTet.minAR, 1.0e-4
 end

 def testEdgeLength
  assert_not_nil grid = Grid.new(2,0,0,0)
  assert_equal 0, grid.addNode(0.0,0.0,0.0)
  assert_equal 1, grid.addNode(0.0,0.0,2.0)
  assert_in_delta 2.0, grid.edgeLength(0,1), 1.0e-15
 end
 
 def testMatrixSpaceMap
  assert_not_nil grid = Grid.new(5,1,0,0)
  assert_equal 0, grid.addNode(-1.0,-1.0,-1.0)
  grid.addCell( grid.addNode(0.0,0.0,0.0), 
	       grid.addNode(1.0,0.0,0.0), 
	       grid.addNode(0.0,2.0,0.0), 
	       grid.addNode(0.0,0.0,5.0) )
  1.upto(4) do |n| 
   assert_equal grid, grid.setMap(n, 1.00, 0.00, 0.00,
	 		                   0.25, 0.00,
			                         0.04)
  end
  assert_equal grid, grid.removeNode(0)
  assert_equal grid, grid.sortNodeGridEx
  assert_in_delta 1.0, grid.edgeRatio(0,1), 1.0e-15
  assert_in_delta 1.0, grid.edgeRatio(0,2), 1.0e-15
  assert_in_delta 1.0, grid.edgeRatio(0,3), 1.0e-15
  assert_in_delta 1.0, grid.edgeLength(0,1), 1.0e-15
  assert_in_delta 2.0, grid.edgeLength(0,2), 1.0e-15
  assert_in_delta 5.0, grid.edgeLength(0,3), 1.0e-15
 end

 def testCopySpacing
  assert_not_nil grid = Grid.new(2,0,0,0)
  assert_equal 0, grid.nnode
  assert_nil      grid.copySpacing(0,1)
  assert_equal 0, grid.addNode(0,0,0)
  assert_nil      grid.copySpacing(0,1)
  assert_equal 1, grid.addNode(1,0,0)
  assert_equal grid, grid.setMap(0, 1, 0, 0, 1, 0, 1)
  assert_equal grid, grid.setMap(1, 2, 0, 0, 2, 0, 2)
  assert_equal 1, grid.spacing(0)
  assert_equal Math::sqrt(0.5), grid.spacing(1)
  assert_equal grid, grid.copySpacing(0,1)
  assert_equal 1, grid.spacing(1)  
 end

 def testAverageSpacing
  assert_not_nil grid = Grid.new(3,0,0,0)
  grid.addNode(0,0,0)
  grid.addNode(0,0,0)
  grid.addNode(0,0,0)
  assert_equal grid, grid.setMap(0, 1, 2, 3, 4, 5, 6)
  assert_equal grid, grid.setMap(1,11,12,13,14,15,16)
  assert_equal grid, grid.setMap(2,20,20,20,20,20,20)
  assert_equal [20,20,20,20,20,20], grid.map(2)
  assert_equal grid, grid.setMapMatrixToAverageOfNodes(2,0,1)
  assert_equal [6,7,8,9,10,11], grid.map(2)
 end

 def testMatrixSpaceMapRotate
  assert_not_nil grid = Grid.new(3,1,0,0)
  assert_equal 0, grid.addNode( 0.0,0.0,0.0)
  assert_equal 1, grid.addNode( 2.0,2.0,0.0)
  assert_equal 2, grid.addNode(-1.0,1.0,0.0)
  sr2  = Math::sqrt(2.0)
  0.upto(2) do |n| 
   assert_equal grid, grid.setMap(n, 0.3125, -0.1875,  0.0000,
	 		                      0.3125,  0.0000,
		                               1.0000)
  end
  assert_equal grid, grid.sortNodeGridEx
  assert_in_delta 2*sr2, grid.edgeLength(0,1), 1.0e-15
  assert_in_delta sr2,   grid.edgeLength(0,2), 1.0e-15
  assert_in_delta 1.0, grid.edgeRatio(0,1), 1.0e-15
  assert_in_delta 1.0, grid.edgeRatio(0,2), 1.0e-15
 end

 def testTriDiagAlreadyDiag
  grid = Grid.new(0,0,0,0)
  m = [1,0,0,
         2,0,
           3]
  d = grid.triDiag(m)
  assert_equal 1, d[0]
  assert_equal 2, d[1]
  assert_equal 3, d[2]
  e = grid.triOffDiag(m)
  assert_equal 0, e[0]
  assert_equal 0, e[1]
  assert_equal 0, e[2]
  q = grid.triDiagTransform(m)
  assert_equal 1, q[0]
  assert_equal 0, q[1]
  assert_equal 0, q[2]

  assert_equal 0, q[3]
  assert_equal 1, q[4]
  assert_equal 0, q[5]

  assert_equal 0, q[6]
  assert_equal 0, q[7]
  assert_equal 1, q[8]
 end

 def testTriDiagAlreadyTriDiag
  grid = Grid.new(0,0,0,0)
  m = [1,4,0,
         2,5,
           3]
  d = grid.triDiag(m)
  assert_equal 1, d[0]
  assert_equal 2, d[1]
  assert_equal 3, d[2]
  e = grid.triOffDiag(m)
  assert_equal 4, e[0]
  assert_equal 5, e[1]
  assert_equal 0, e[2]
  q = grid.triDiagTransform(m)
  assert_equal 1, q[0]
  assert_equal 0, q[1]
  assert_equal 0, q[2]

  assert_equal 0, q[3]
  assert_equal 1, q[4]
  assert_equal 0, q[5]

  assert_equal 0, q[6]
  assert_equal 0, q[7]
  assert_equal 1, q[8]
 end

 def testTriDiagFull
  grid = Grid.new(0,0,0,0)
  m = [1,3,4,
         2,5,
           6]
  d = grid.triDiag(m)
  assert_equal 1, d[0]
  assert_equal 2+0.8*(2*0.6*5+0.8*(6-2)), d[1]
  assert_in_delta 6-0.8*(2*0.6*5+0.8*(6-2)), d[2], 1.0e-15
  e = grid.triOffDiag(m)
  assert_equal 5, e[0]
  assert_in_delta 5-0.6*(2*0.6*5+0.8*(6-2)), e[1], 1.0e-15
  assert_equal 0, e[2]
  q = grid.triDiagTransform(m)
  assert_equal 1, q[0]
  assert_equal 0, q[1]
  assert_equal 0, q[2]

  assert_equal 0, q[3]
  assert_equal 0.6, q[4]
  assert_equal 0.8, q[5]

  assert_equal 0, q[6]
  assert_equal 0.8, q[7]
  assert_equal(-0.6, q[8])
 end

 def testEigTriDiagAlreadyDiag111
  grid = Grid.new(0,0,0,0)
  d = [1,1,1]
  e = [0,0,0]
  q = [1,0,0, 0,1,0, 0,0,1]
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 1, eig[0]
  assert_equal 1, eig[1]
  assert_equal 1, eig[2]
  vect = grid.vectTriDiag(d,e,q)
  assert_equal 1, vect[0]
  assert_equal 0, vect[1]
  assert_equal 0, vect[2]

  assert_equal 0, vect[3]
  assert_equal 1, vect[4]
  assert_equal 0, vect[5]

  assert_equal 0, vect[6]
  assert_equal 0, vect[7]
  assert_equal 1, vect[8]
 end

 def testEigTriDiagAlreadyDiag132
  grid = Grid.new(0,0,0,0)
  d = [1,3,2]
  e = [0,0,0]
  q = [1,0,0, 0,1,0, 0,0,1]
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 3, eig[0]
  assert_equal 2, eig[1]
  assert_equal 1, eig[2]
  vect = grid.vectTriDiag(d,e,q)
  assert_equal 0, vect[0]
  assert_equal 1, vect[1]
  assert_equal 0, vect[2]

  assert_equal 0, vect[3]
  assert_equal 0, vect[4]
  assert_equal 1, vect[5]

  assert_equal 1, vect[6]
  assert_equal 0, vect[7]
  assert_equal 0, vect[8]
 end

 def XtestEigTriDiag2x2plus1
  grid = Grid.new(0,0,0,0)
  m = [0.5,0.5,0,
           3.5,0,
               1]
  d = grid.triDiag(m)
  e = grid.triOffDiag(m)
  q = grid.triDiagTransform(m)
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 3, eig[0]
  assert_equal 2, eig[1]
  assert_equal 1, eig[2]
 end

 def XtestEigTriDiagFull211212
  grid = Grid.new(0,0,0,0)
  m = [2,1,1,
         2,1,
           2]
  d = grid.triDiag(m)
  e = grid.triOffDiag(m)
  q = grid.triDiagTransform(m)
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 4, eig[0]
  assert_in_delta 1, eig[1], 1.0e-15
  assert_in_delta 1, eig[2], 1.0e-15
  vect = grid.vectTriDiag(d,e,q)
  invsqrt2 = 0.707106781186547
  assert_in_delta(0.57735027,vect[0],1.0e-7)
  assert_in_delta(0.57735027,vect[3],1.0e-7)
  assert_in_delta(0.57735027,vect[6],1.0e-7)
  assert_in_delta(0,vect[1],1.0e-7)
  assert_in_delta(invsqrt2,vect[4],1.0e-7)
  assert_in_delta(-invsqrt2,vect[7],1.0e-7)
  assert_in_delta( 0.81649658,vect[2],1.0e-7)
  assert_in_delta(-0.40824829,vect[5],1.0e-7)
  assert_in_delta(-0.40824829,vect[7],1.0e-7)
 end

 def XtestEigTriDiagAlreadyTridiagFull13n471
  grid = Grid.new(0,0,0,0)
  m = [13,-4,0,
           7,0,
             1]
  d = grid.triDiag(m)
  e = grid.triOffDiag(m)
  q = grid.triDiagTransform(m)
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 15, eig[0]
  assert_in_delta 5, eig[1], 1.0e-15
  assert_in_delta 1, eig[2], 1.0e-15
 end

 def XtestEigTriDiagAlreadyTridiagFull1554
  grid = Grid.new(0,0,0,0)
  m = [4, 0, 0,
         13,-7,
             7]
  d = grid.triDiag(m)
  e = grid.triOffDiag(m)
  q = grid.triDiagTransform(m)
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 15, eig[0]
  assert_in_delta 5, eig[1], 1.0e-15
  assert_in_delta 4, eig[2], 1.0e-15
 end

 def XtestEigTriDiagAlreadyTridiagSplit1554
  grid = Grid.new(0,0,0,0)
  m = [13,0,-4,
          4, 0,
             7]
  d = grid.triDiag(m)
  e = grid.triOffDiag(m)
  q = grid.triDiagTransform(m)
  eig = grid.eigTriDiag(d,e,q)
  assert_equal 15, eig[0]
  assert_in_delta 5, eig[1], 1.0e-15
  assert_in_delta 4, eig[2], 1.0e-15
 end

 def testMatrixEigenValues
  assert_not_nil grid = Grid.new(1,0,0,0)

  assert_not_nil eigsys = grid.eigenSystem( [ 1.0, 0.0, 0.0, 
                                                   2.0, 0.0, 
                                                        3.0 ])
  eig = eigsys[0]
  assert_in_delta 3.0, eig[0], 1.0e-15
  assert_in_delta 2.0, eig[1], 1.0e-15
  assert_in_delta 1.0, eig[2], 1.0e-15

  assert_not_nil eigsys = grid.eigenSystem( [ 0.3125, -0.1875,  0.0000,
                                                       0.3125,  0.0000,
                                                                1.0000 ] )
  eig = eigsys[0]
  assert_in_delta 1.000, eig[0], 1.0e-15
  assert_in_delta 0.500, eig[1], 1.0e-15
  assert_in_delta 0.125, eig[2], 1.0e-15

  assert_not_nil eigsys = grid.eigenSystem( [ 1.0, 0.0, 0.0, 
                                                   1.0, 0.0, 
                                                        1.0 ])
  eig = eigsys[0]
  assert_in_delta 1.0, eig[0], 1.0e-15
  assert_in_delta 1.0, eig[1], 1.0e-15
  assert_in_delta 1.0, eig[2], 1.0e-15

  assert_not_nil eigsys = grid.eigenSystem( [ 2.0, 0.0, 0.0, 
                                                   1.0, 0.0, 
                                                        1.0 ])
  eig = eigsys[0]
  assert_in_delta 2.0, eig[0], 1.0e-15
  assert_in_delta 1.0, eig[1], 1.0e-15
  assert_in_delta 1.0, eig[2], 1.0e-15

  assert_not_nil eigsys = grid.eigenSystem( 
[ 5669.182266666660325,    0.000000000000379,    0.000000000000497,
                        5669.182266666660325,    0.000000000000436,
                                              5669.182266666660325])
  ans = 5669.18226
  eig = eigsys[0]
  assert_in_delta ans, eig[0], 1.0e-5
  assert_in_delta ans, eig[1], 1.0e-5
  assert_in_delta ans, eig[2], 1.0e-5

 end

 def testMatrixEigenSystem
  assert_not_nil grid = Grid.new(1,0,0,0)

  assert_not_nil eigsys = grid.eigenSystem( [ 1.0, 0.0, 0.0, 
                                                   2.0, 0.0, 
                                                        3.0 ])
  assert_not_nil eig = eigsys[0]
  assert_not_nil v1 = eigsys[1]
  assert_not_nil v2 = eigsys[2]
  assert_not_nil v3 = eigsys[3]

  assert_in_delta 3.0, eig[0], 1.0e-15
  assert_in_delta 2.0, eig[1], 1.0e-15
  assert_in_delta 1.0, eig[2], 1.0e-15

  assert_in_delta 0.0, v1[0], 1.0e-15
  assert_in_delta 0.0, v1[1], 1.0e-15
  assert_in_delta 1.0, v1[2], 1.0e-15
  assert_in_delta 0.0, v2[0], 1.0e-15
  assert_in_delta 1.0, v2[1], 1.0e-15
  assert_in_delta 0.0, v2[2], 1.0e-15
  assert_in_delta 1.0, v3[0], 1.0e-15
  assert_in_delta 0.0, v3[1], 1.0e-15
  assert_in_delta 0.0, v3[2], 1.0e-15

  assert_not_nil eigsys = grid.eigenSystem( [ 0.3125, -0.1875,  0.0000,
                                                       0.3125,  0.0000,
                                                                1.0000 ] )
  assert_not_nil eig = eigsys[0]
  assert_not_nil v1 = eigsys[1]
  assert_not_nil v2 = eigsys[2]
  assert_not_nil v3 = eigsys[3]

  assert_in_delta 1.000, eig[0], 1.0e-15
  assert_in_delta 0.500, eig[1], 1.0e-15
  assert_in_delta 0.125, eig[2], 1.0e-15

  os2 = 1.0/Math::sqrt(2.0)

  assert_in_delta 0.0, v1[0], 1.0e-15
  assert_in_delta 0.0, v1[1], 1.0e-15
  assert_in_delta( -1.0, v1[2], 1.0e-15)
  assert_in_delta( -os2, v2[0], 1.0e-15)
  assert_in_delta os2, v2[1], 1.0e-15
  assert_in_delta 0.0, v2[2], 1.0e-15
  assert_in_delta( -os2, v3[0], 1.0e-15)
  assert_in_delta( -os2, v3[1], 1.0e-15)
  assert_in_delta 0.0, v3[2], 1.0e-15
 end

 def testConvertMetricToJacobian
  assert_not_nil grid = Grid.new(1,0,0,0)
  assert_not_nil jacob = grid.convertMetricToJacobian( [ 1.0, 0.0, 0.0, 
                                                              4.0, 0.0, 
                                                                   9.0 ] )
  assert_in_delta 0.0, jacob[0], 1.0e-15
  assert_in_delta 0.0, jacob[1], 1.0e-15
  assert_in_delta 3.0, jacob[2], 1.0e-15
  assert_in_delta 0.0, jacob[3], 1.0e-15
  assert_in_delta 2.0, jacob[4], 1.0e-15
  assert_in_delta 0.0, jacob[5], 1.0e-15
  assert_in_delta(-1.0, jacob[6], 1.0e-15)
  assert_in_delta 0.0, jacob[7], 1.0e-15
  assert_in_delta 0.0, jacob[8], 1.0e-15

  assert_not_nil jacob = grid.convertMetricToJacobian([0.3125, -0.1875,  0.0,
                                                                0.3125,  0.0,
                                                                         1.0])
  assert_in_delta 0.0, jacob[0], 1.0e-15
  assert_in_delta 0.0, jacob[1], 1.0e-15
  assert_in_delta(-1.0, jacob[2], 1.0e-15)
  assert_in_delta(-0.5, jacob[3], 1.0e-15)
  assert_in_delta 0.5, jacob[4], 1.0e-15
  assert_in_delta 0.0, jacob[5], 1.0e-15
  assert_in_delta 0.25, jacob[6], 1.0e-15
  assert_in_delta 0.25, jacob[7], 1.0e-15
  assert_in_delta 0.0, jacob[8], 1.0e-15

 end

 def testFindLargestRatioEdge
  assert_not_nil grid = isoTet.resetSpacing
  grid.scaleSpacing(0,0.50)
  grid.scaleSpacing(1,0.95)

  assert_equal 1, grid.largestRatioEdge(0)
  assert_equal 0, grid.largestRatioEdge(1)
  assert_equal 0, grid.largestRatioEdge(2)
  assert_equal 0, grid.largestRatioEdge(3)
 end

 def testFindSmallestRatioEdge
  assert_not_nil grid = isoTet.resetSpacing
  grid.scaleSpacing(0,2.00)
  grid.scaleSpacing(1,1.05)

  assert_equal 1, grid.smallestRatioEdge(0)
  assert_equal 0, grid.smallestRatioEdge(1)
  assert_equal 0, grid.smallestRatioEdge(2)
  assert_equal 0, grid.smallestRatioEdge(3)
 end

# need this?
 def testAverageEdgeLength
  assert_not_nil grid = isoTet
  4.times do |i| 
   assert_in_delta 1.0, grid.averageEdgeLength(i), 3.0e-3, "node #{i} length"
  end
 end

 def testSpacingFunction
  assert_not_nil grid = isoTet
  4.times do |i| 
   assert_in_delta 1.0, grid.spacing(i), 1.0e-15, "node #{i} spacing"
  end
  assert_equal grid, grid.resetSpacing
  4.times do |i| 
   assert_in_delta grid.averageEdgeLength(i), grid.spacing(i), 1.0e-15, "n#{i}"
  end
  assert_equal grid, grid.scaleSpacing(0,0.5)
  assert_in_delta grid.averageEdgeLength(0)*0.5, grid.spacing(0), 1.0e-15
 end

 def testScaleSpacingFunctionSphere
  assert_not_nil grid = rightTet
  spacing = grid.spacing(0)
  assert_equal grid, grid.scaleSpacingSphere(1.0,0.0,0.0,1.1,0.5)
  delta = 1.0e-15
  assert_in_delta spacing*0.5, grid.spacing(0), delta
  assert_in_delta spacing*0.5, grid.spacing(1), delta
  assert_in_delta spacing,     grid.spacing(2), delta
  assert_in_delta spacing,     grid.spacing(3), delta
 end

 def testVolumeMetrics
  assert_not_nil grid = rightTet
  nodes = [0,1,2,3]
  ar = 0.8399473666
  assert_in_delta 1.0/6.0, grid.volume(nodes), 1.0e-15
  assert_in_delta 1.0/6.0, grid.minVolume, 1.0e-15
  assert_in_delta ar, grid.ar(nodes), 1.0e-10
  assert_in_delta ar, grid.minAR, 1.0e-10
  assert_in_delta ar, grid.nodeAR(0), 1.0e-10
  assert_in_delta ar, grid.nodeAR(1), 1.0e-10
  assert_in_delta ar, grid.nodeAR(2), 1.0e-10
  assert_in_delta ar, grid.nodeAR(3), 1.0e-10
 end

 def testARDerivatives
  assert_not_nil grid = Grid.new(4,1,0,0)

  grid.addCell( 
	       grid.addNode(0.0,0.0,0.0), 
	       grid.addNode(1.0,0.0,0.0), 
	       grid.addNode(0.0,1.0,0.0), 
	       grid.addNode(0.0,0.0,1.0) )
  nodes = [0,1,2,3]
  ar    = 0.8399473666
  tol   = 1.0e-10
  deriv = -0.3733099407

  ans = grid.cellARDerivative(nodes)
  assert_in_delta ar,    ans[0], tol
  assert_in_delta deriv, ans[1], tol
  assert_in_delta deriv, ans[2], tol
  assert_in_delta deriv, ans[3], tol  

  ans = grid.nodeARDerivative(0)
  assert_in_delta ar,    ans[0], tol
  assert_in_delta deriv, ans[1], tol
  assert_in_delta deriv, ans[2], tol
  assert_in_delta deriv, ans[3], tol  

  deriv = 0.1866549704

  ans = grid.nodeARDerivative(1)
  assert_in_delta ar,    ans[0], tol
  assert_in_delta 0.0,   ans[1], tol
  assert_in_delta deriv, ans[2], tol
  assert_in_delta deriv, ans[3], tol  

  ans = grid.nodeARDerivative(2)
  assert_in_delta ar,    ans[0], tol
  assert_in_delta deriv, ans[1], tol
  assert_in_delta 0.0,   ans[2], tol
  assert_in_delta deriv, ans[3], tol  

  ans = grid.nodeARDerivative(3)
  assert_in_delta ar,    ans[0], tol
  assert_in_delta deriv, ans[1], tol
  assert_in_delta deriv, ans[2], tol
  assert_in_delta 0.0,   ans[3], tol  

  assert_nil             grid.storeARDerivative(10)
  assert_equal grid,     grid.storeARDerivative(0)
 end

 def testStoreARDerivative
  assert_not_nil grid = Grid.new(5,2,0,0)
  grid.addCell( 
	       grid.addNode(0.0,0.0,0.0), 
	       grid.addNode(1.0,0.0,0.0), 
	       grid.addNode(0.0,1.0,0.0), 
	       grid.addNode(0.0,0.0,1.0) )
  grid.addCell( 
	       1, 
	       2, 
	       3, 
	       grid.addNode(0.7,0.7,0.7) )
  assert_equal 0,        grid.storedARDegree
  assert_nil             grid.storeARDerivative(10)
  assert_equal 1,        grid.cellDegree(0)
  assert_equal grid,     grid.storeARDerivative(0)
  assert_equal 1,        grid.storedARDegree
  assert_equal 2,        grid.cellDegree(1)
  assert_equal grid,     grid.storeARDerivative(1)
  assert_equal 2,        grid.storedARDegree
 end

 def testRightHandedFaces
  assert_not_nil grid = Grid.new(4,1,2,0)

  grid.addCell( grid.addNode(0.0,0.0,0.0), grid.addNode(1.0,0.0,0.0), 
	        grid.addNode(0.0,1.0,0.0), grid.addNode(0.0,0.0,1.0) )
  grid.addFace(0,1,2,11)
  assert_equal true,  grid.rightHandedFace(0)
  grid.addFace(0,2,3,11)
  assert_equal true,  grid.rightHandedBoundary
  assert_equal grid,  grid.removeFace(grid.findFace(0,1,2))
  grid.addFace(0,2,1,11)
  assert_equal false, grid.rightHandedFace(0)
  assert_equal false, grid.rightHandedBoundary
 end

 def testFaceMetrics
  assert_not_nil grid = rightTet
  assert_in_delta 0.5, grid.faceArea(0,1,2), 1.0e-15
  assert_in_delta 1.0, grid.faceAR(1,2,3), 1.0e-15
  ar = 0.8284271247
  assert_in_delta ar, grid.faceAR(0,1,2), 1.0e-8
  mr = 0.8660254038
  assert_in_delta 1.0, grid.faceMR(1,2,3), 1.0e-14
  assert_in_delta mr,  grid.faceMR(0,1,2), 1.0e-8

  deriv = -0.433
  ans = grid.faceMRDerivative([0,1,2])
  assert_in_delta mr,    ans[0], 1.0e-8
  assert_in_delta deriv, ans[1], 1.0e-4
  assert_in_delta deriv, ans[2], 1.0e-4
  assert_in_delta 0.0,   ans[3], 1.0e-14

  grid.addFace(0,1,2,10)
  grid.addFace(1,2,3,10)
  assert_in_delta mr,  grid.nodeFaceMR(0), 1.0e-8
  assert_in_delta mr,  grid.nodeFaceMR(1), 1.0e-8
  assert_in_delta 1.0, grid.nodeFaceMR(3), 1.0e-8
  assert_in_delta mr,  grid.minFaceMR, 1.0e-8

  ans = grid.nodeFaceMRDerivative(0)
  assert_in_delta mr,    ans[0], 1.0e-8
  assert_in_delta deriv, ans[1], 1.0e-4
  assert_in_delta deriv, ans[2], 1.0e-4
  assert_in_delta 0.0,   ans[3], 1.0e-15  
  ans = grid.nodeFaceMRDerivative(1)
  assert_in_delta mr,    ans[0], 1.0e-8
  assert_in_delta 0.0, ans[1], 1.0e-4
  assert_in_delta( -deriv, ans[2], 1.0e-4)
  assert_in_delta 0.0,   ans[3], 1.0e-15  
  ans = grid.nodeFaceMRDerivative(3)
  assert_in_delta 1.0,    ans[0], 1.0e-8
  assert_in_delta 0.0, ans[1], 1.0e-4
  assert_in_delta 0.0, ans[2], 1.0e-4
  assert_in_delta 0.0,   ans[3], 1.0e-15  
 end
 
 def testFaceMRDerivative
  x1 =0.061
  y1 =0.045
  z1 =0.077

  x2 =1.561
  y2 =0.145
  z2 =0.277

  x3 =-0.161
  y3 =2.545
  z3 =-0.0177

  assert_not_nil grid = Grid.new(3,0,0,0)
  assert_equal 0, grid.addNode(x1,y1,z1)
  assert_equal 1, grid.addNode(x2,y2,z2)
  assert_equal 2, grid.addNode(x3,y3,z3)

  delta = 1.0e-7
  ans = grid.FaceMRDerivative(x1,y1,z1,x2,y2,z2,x3,y3,z3)

  dx = grid.FaceMRDerivative(x1+delta,y1,z1,x2,y2,z2,x3,y3,z3)
  dx = (dx[0]-ans[0])/delta
  assert_in_delta dx, ans[1], 10.0*delta, "dx"

  dy = grid.FaceMRDerivative(x1,y1+delta,z1,x2,y2,z2,x3,y3,z3)
  dy = (dy[0]-ans[0])/delta
  assert_in_delta dy, ans[2], 10.0*delta, "dy"

  dz = grid.FaceMRDerivative(x1,y1,z1+delta,x2,y2,z2,x3,y3,z3)
  dz = (dz[0]-ans[0])/delta
  assert_in_delta dz, ans[3], 10.0*delta, "dz"

  assert_in_delta grid.faceMR(0,1,2),ans[0],1e-14
 end

 def testCellMR

  nodes = [0,1,2,3]

  grid = rightTet
  node0 = grid.nodeXYZ(0)
  node1 = grid.nodeXYZ(1)
  node2 = grid.nodeXYZ(2)
  node3 = grid.nodeXYZ(3)
  assert_in_delta 0.840, grid.cellMeanRatio(node0,node1,node2,node3), 1.0e-4

  grid = isoTet
  node0 = grid.nodeXYZ(0)
  node1 = grid.nodeXYZ(1)
  node2 = grid.nodeXYZ(2)
  node3 = grid.nodeXYZ(3)
  assert_in_delta 1.0, grid.cellMeanRatio(node0,node1,node2,node3), 1.0e-4

 end

 def testCellMRDerivative

  x0=0.1
  y0=0.2
  z0=0.3
  node0 = [x0,y0,z0]
  node1 = [1.12,-0.1,-0.2]
  node2 = [0.01,1.05,0.12]
  node3 = [0.34,0.22,0.99]
  assert_not_nil grid = Grid.new(4,0,0,0)
  grid.addNode(node0[0],node0[1],node0[2])
  grid.addNode(node1[0],node1[1],node1[2])
  grid.addNode(node2[0],node2[1],node2[2])
  grid.addNode(node3[0],node3[1],node3[2])

  delta = 1.0e-7
  ans = grid.cellMeanRatioDerivative(node0,node1,node2,node3)

  dx = grid.cellMeanRatioDerivative([x0+delta,y0,z0],node1,node2,node3)
  dx = (dx[0]-ans[0])/delta
  assert_in_delta dx, ans[1], 10.0*delta, "dx"

  dy = grid.cellMeanRatioDerivative([x0,y0+delta,z0],node1,node2,node3)
  dy = (dy[0]-ans[0])/delta
  assert_in_delta dy, ans[2], 10.0*delta, "dy"

  dz = grid.cellMeanRatioDerivative([x0,y0,z0+delta],node1,node2,node3)
  dz = (dz[0]-ans[0])/delta
  assert_in_delta dz, ans[3], 10.0*delta, "dz"

  assert_in_delta grid.cellMeanRatio(node0,node1,node2,node3), 
   ans[0], 1.0e-15, "function"
 

 end

end
