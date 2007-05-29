#!/usr/bin/env ruby
#
# Mobility test for gridmath c lib
#
# $Id$

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('GridMath').build

require 'test/unit'
require 'GridMath/GridMath'

class Array
 def mat_show
   puts
  if self.length==12
   printf "%10.5f %10.5f %10.5f %10.5f\n", self[0], self[3], self[6], self[9]
   printf "%10.5f %10.5f %10.5f %10.5f\n", self[1], self[4], self[7], self[10]
   printf "%10.5f %10.5f %10.5f %10.5f\n", self[2], self[5], self[8], self[11]
  else
   printf "%10.5f %10.5f %10.5f\n", self[0], self[2], self[4]
   printf "%10.5f %10.5f %10.5f\n", self[1], self[3], self[5]
  end
 end
end

class TestGridMath < Test::Unit::TestCase

 def set_up
  @gm = GridMath.new
  @sin45 = 0.707106781186547
  @sin30 = 0.5
  @cos30 = 0.866025403784439
 end
 def setup ; set_up ; end

 def testSubtractVector
  v1 = [1,2,3]
  v2 = [7,0,6]
  assert_equal [-6,2,-3], @gm.subtractVector(v1,v2)
  v1 = [13,6,23]
  v2 = [2,10,16]
  assert_equal [11,-4,7], @gm.subtractVector(v1,v2)
 end

 def testDotProduct
  v1 = [1,2,3]
  v2 = [1,2,3]
  assert_equal 14, @gm.dotProduct(v1,v2)
  v1 = [2,0,5]
  v2 = [1,20,3]
  assert_equal 17, @gm.dotProduct(v1,v2)
 end

 def testCrossProduct
  v1 = [1,0,0]
  v2 = [0,1,0]
  assert_equal [0,0,1],  @gm.crossProduct(v1,v2)
  assert_equal [0,0,-1], @gm.crossProduct(v2,v1)
  assert_equal [0,0,0], @gm.crossProduct(v1,v1)
  assert_equal [0,0,0], @gm.crossProduct(v2,v2)
 end

 def testVectorLength
  assert_equal 5, @gm.vectorLength([3,4,0])
  assert_equal 0, @gm.vectorLength([0,0,0])
  assert_equal 1, @gm.vectorLength([1,0,0])
  assert_equal 10, @gm.vectorLength([0,10,0])
 end

 def testNormalizeVector
  assert_equal [1,0,0], @gm.vectorNormalize([5,0,0])
  assert_equal [0,1,0], @gm.vectorNormalize([0,6,0])
  assert_equal [0,0,-1], @gm.vectorNormalize([0,0,-3])
 end

 def testOrthogonalizeVector
  axle = [0,0,1]
  assert_equal [1,0,0], @gm.vectorOrthogonalize([1,0,0],axle)
  assert_equal [0,1,0], @gm.vectorOrthogonalize([0,1,0],axle)
  assert_equal [1,0,0], @gm.vectorOrthogonalize([1,0,1],axle)
  assert_equal [0,1,0], @gm.vectorOrthogonalize([0,1,1],axle)
 end

 def test_BarycentricCoordinate_of_vertices
  tol =1.0e-14
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]
  xyz3 = [0,0,1]

  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,xyz0)
  ans = [1.0,0.0,0.0,0.0]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }

  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,xyz1)
  ans = [0.0,1.0,0.0,0.0]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }

  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,xyz2)
  ans = [0.0,0.0,1.0,0.0]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }

  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,xyz3)
  ans = [0.0,0.0,0.0,1.0]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }
 end

 def test_BarycentricCoordinate_of_center1
  tol = 1.0e-14
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]
  xyz3 = [0,0,1]
  target = [0.25,0.25,0.25]
  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,target)
  ans = [0.25,0.25,0.25,0.25]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }
 end

 def test_BarycentricCoordinate_of_center2
  tol = 1.0e-14
  xyz0 = [0,0,0]
  xyz1 = [4,0,0]
  xyz2 = [0,4,0]
  xyz3 = [0,0,4]
  target = [1.0,1.0,1.0]
  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,target)
  ans = [0.25,0.25,0.25,0.25]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }
 end

 def test_BarycentricCoordinate_of_outside_x
  tol = 1.0e-14
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]
  xyz3 = [0,0,1]
  target = [2.0,0.0,0.0]
  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,target)
  ans = [-1.0,2.0,0.0,0.0]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }
 end

 def test_BarycentricCoordinate_of_oppsite_0
  tol = 1.0e-14
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]
  xyz3 = [0,0,1]
  target = [1.0,1.0,1.0]
  bary = @gm.barycentricCoordinate(xyz0,xyz1,xyz2,xyz3,target)
  ans = [-2.0,1.0,1.0,1.0]
  4.times { |i| assert_in_delta( ans[i], bary[i], tol, 'element '+i.to_s ) }
 end

 def testRotateDirectionEndPoints
  v0 = [1,0,0]
  v1 = [0,1,0]
  axle = [0,0,1]
  tol = 1.0e-12
  result = @gm.rotateDirection(v0,v1,axle,0)
  assert_in_delta v0[0], result[0], tol
  assert_in_delta v0[1], result[1], tol
  assert_in_delta v0[2], result[2], tol
  result = @gm.rotateDirection(v0,v1,axle,1)
  assert_in_delta v1[0], result[0], tol
  assert_in_delta v1[1], result[1], tol
  assert_in_delta v1[2], result[2], tol
  negaxle = [0,0,-1]
  result = @gm.rotateDirection(v1,v0,negaxle,1)
  assert_in_delta v0[0], result[0], tol
  assert_in_delta v0[1], result[1], tol
  assert_in_delta v0[2], result[2], tol
  result = @gm.rotateDirection(v1,v0,negaxle,0)
  assert_in_delta v1[0], result[0], tol
  assert_in_delta v1[1], result[1], tol
  assert_in_delta v1[2], result[2], tol
 end

 def testRotateDirectionMiddleOrthog
  v0 = [1,0,0]
  v1 = [0,1,0]
  axle = [0,0,1]
  result = @gm.rotateDirection(v0,v1,axle,0.5)
  assert_in_delta @sin45, result[0], 1.0e-7
  assert_in_delta @sin45, result[1], 1.0e-7
  assert_in_delta 0,      result[2], 1.0e-7
  result = @gm.rotateDirection(v0,v1,axle,(1.0/3.0))
  assert_in_delta @cos30, result[0], 1.0e-7
  assert_in_delta @sin30, result[1], 1.0e-7
  assert_in_delta 0,      result[2], 1.0e-7
  result = @gm.rotateDirection(v0,v1,axle,(2.0/3.0))
  assert_in_delta @sin30, result[0], 1.0e-7
  assert_in_delta @cos30, result[1], 1.0e-7
  assert_in_delta 0,      result[2], 1.0e-7
 end

 def testRotateDirectionMiddleObtose
  v0 = [1,0,0]
  v1 = [-@sin45,@sin45,0]
  axle = [0,0,1]
  result = @gm.rotateDirection(v0,v1,axle,(2.0/3.0))
  assert_in_delta 0, result[0], 1.0e-7
  assert_in_delta 1, result[1], 1.0e-7
  assert_in_delta 0, result[2], 1.0e-7
  result = @gm.rotateDirection(v0,v1,axle,(1.0/3.0))
  assert_in_delta @sin45, result[0], 1.0e-7
  assert_in_delta @sin45, result[1], 1.0e-7
  assert_in_delta 0,      result[2], 1.0e-7
 end

 def testRotateDirectionMiddleOrthogWithEqualSkew
  v0 = [@sin45,0,@sin45]
  v1 = [0,@sin45,-@sin45]
  axle = [0,0,1]
  result = @gm.rotateDirection(v0,v1,axle,0.0)
  assert_in_delta v0[0], result[0], 1.0e-7
  assert_in_delta v0[1], result[1], 1.0e-7
  assert_in_delta v0[2], result[2], 1.0e-7
  result = @gm.rotateDirection(v0,v1,axle,0.5)
  assert_in_delta @sin45, result[0], 1.0e-7
  assert_in_delta @sin45, result[1], 1.0e-7
  assert_in_delta 0,      result[2], 1.0e-7
  result = @gm.rotateDirection(v0,v1,axle,1.0)
  assert_in_delta v1[0], result[0], 1.0e-7
  assert_in_delta v1[1], result[1], 1.0e-7
  assert_in_delta v1[2], result[2], 1.0e-7
 end

 def testRotateDirectionMiddleOrthogWithUnequalSkew
  v0 = [@cos30,0,-@sin30]
  v1 = [0,@sin30,@cos30]
  axle = [0,0,1]
  result = @gm.rotateDirection(v0,v1,axle,(1.0/3.0))
  assert_in_delta 0.0,   result[2], 1.0e-7
  result = @gm.rotateDirection(v0,v1,axle,(2.0/3.0))
  assert_in_delta @sin30,   result[2], 1.0e-7
 end

 def testTriDiagAlreadyDiag
  m = [1,0,0,
         2,0,
           3]
  d = @gm.triDiag(m)
  assert_equal 1, d[0]
  assert_equal 2, d[1]
  assert_equal 3, d[2]
  e = @gm.triOffDiag(m)
  assert_equal 0, e[0]
  assert_equal 0, e[1]
  assert_equal 0, e[2]
  q = @gm.triDiagTransform0(m)
  assert_equal 1, q[0]
  assert_equal 0, q[1]
  assert_equal 0, q[2]
  q = @gm.triDiagTransform1(m)
  assert_equal 0, q[0]
  assert_equal 1, q[1]
  assert_equal 0, q[2]
  q = @gm.triDiagTransform2(m)
  assert_equal 0, q[0]
  assert_equal 0, q[1]
  assert_equal 1, q[2]
 end

 def testTriDiagAlreadyTriDiag
  m = [1,4,0,
         2,5,
           3]
  d = @gm.triDiag(m)
  assert_equal 1, d[0]
  assert_equal 2, d[1]
  assert_equal 3, d[2]
  e = @gm.triOffDiag(m)
  assert_equal 4, e[0]
  assert_equal 5, e[1]
  assert_equal 0, e[2]
  q = @gm.triDiagTransform0(m)
  assert_equal 1, q[0]
  assert_equal 0, q[1]
  assert_equal 0, q[2]
  q = @gm.triDiagTransform1(m)
  assert_equal 0, q[0]
  assert_equal 1, q[1]
  assert_equal 0, q[2]
  q = @gm.triDiagTransform2(m)
  assert_equal 0, q[0]
  assert_equal 0, q[1]
  assert_equal 1, q[2]
 end

 def testTriDiagFull
  m = [1,3,4,
         2,5,
           6]
  d = @gm.triDiag(m)
  assert_equal 1, d[0]
  assert_equal 2+0.8*(2*0.6*5+0.8*(6-2)), d[1]
  assert_in_delta 6-0.8*(2*0.6*5+0.8*(6-2)), d[2], 1.0e-15
  e = @gm.triOffDiag(m)
  assert_equal 5, e[0]
  assert_in_delta 5-0.6*(2*0.6*5+0.8*(6-2)), e[1], 1.0e-15
  assert_equal 0, e[2]
  q = @gm.triDiagTransform0(m)
  assert_equal 1, q[0]
  assert_equal 0, q[1]
  assert_equal 0, q[2]
  q = @gm.triDiagTransform1(m)
  assert_equal 0,   q[0]
  assert_equal 0.6, q[1]
  assert_equal 0.8, q[2]
  q = @gm.triDiagTransform2(m)
  assert_equal 0,    q[0]
  assert_equal 0.8,  q[1]
  assert_equal(-0.6, q[2])
 end

 def testEigTriDiagAlreadyDiag111
  d = [1,1,1]
  e = [0,0,0]
  q0 = [1,0,0]
  q1 = [0,1,0]
  q2 = [0,0,1]
  eig = @gm.eigTriDiag(d,e)
  assert_equal 1, eig[0]
  assert_equal 1, eig[1]
  assert_equal 1, eig[2]
  v0 = @gm.vectTriDiag0(d,e,q0,q1,q2)
  assert_equal 1, v0[0]
  assert_equal 0, v0[1]
  assert_equal 0, v0[2]
  v1 = @gm.vectTriDiag1(d,e,q0,q1,q2)
  assert_equal 0, v1[0]
  assert_equal 1, v1[1]
  assert_equal 0, v1[2]
  v2 = @gm.vectTriDiag2(d,e,q0,q1,q2)
  assert_equal 0, v2[0]
  assert_equal 0, v2[1]
  assert_equal 1, v2[2]
 end

 def testEigTriDiagAlreadyDiag132
  d = [1,3,2]
  e = [0,0,0]
  q0 = [1,0,0]
  q1 = [0,1,0]
  q2 = [0,0,1]
  eig = @gm.eigTriDiag(d,e)
  assert_equal 3, eig[0]
  assert_equal 2, eig[1]
  assert_equal 1, eig[2]
  v0 = @gm.vectTriDiag0(d,e,q0,q1,q2)
  assert_equal 0, v0[0]
  assert_equal 1, v0[1]
  assert_equal 0, v0[2]
  v1 = @gm.vectTriDiag1(d,e,q0,q1,q2)
  assert_equal 0, v1[0]
  assert_equal 0, v1[1]
  assert_equal 1, v1[2]
  v2 = @gm.vectTriDiag2(d,e,q0,q1,q2)
  assert_equal 1, v2[0]
  assert_equal 0, v2[1]
  assert_equal 0, v2[2]
 end

 def testEigTriDiag2x2plus1
  m = [0.5,0.5,0,
           3.5,0,
               1]
  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta 3.58113883, eig[0], 1.0e-7
  assert_in_delta 1,          eig[1], 1.0e-7
  assert_in_delta 0.41886117, eig[2], 1.0e-7
 end

 def testEigTriDiagFull211212
  tol = 1.0e-7
  m = [2,1,1,
         2,1,
           2]
  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  q0 = @gm.triDiagTransform0(m)
  q1 = @gm.triDiagTransform1(m)
  q2 = @gm.triDiagTransform2(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta 4, eig[0], 1.0e-15
  assert_in_delta 1, eig[1], 1.0e-15
  assert_in_delta 1, eig[2], 1.0e-15
  v0 = @gm.vectTriDiag0(d,e,q0,q1,q2)
  v1 = @gm.vectTriDiag1(d,e,q0,q1,q2)
  v2 = @gm.vectTriDiag2(d,e,q0,q1,q2)
  invsqrt2 = 0.707106781186547
  assert_in_delta( 0.57735027, v0[0],tol)
  assert_in_delta( 0.57735027, v0[1],tol)
  assert_in_delta( 0.57735027, v0[2],tol)
  assert_in_delta( 0.81649658, v2[0],tol)
  assert_in_delta(-0.40824829, v2[1],tol)
  assert_in_delta(-0.40824829, v2[2],tol)
  assert_in_delta( 0.0,        v1[0],tol)
  assert_in_delta( invsqrt2,   v1[1],tol)
  assert_in_delta(-invsqrt2,   v1[2],tol)
 end

 def testEigTriDiag13n471
  m = [13,-4,0,
           7,0,
             1]
  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta 15, eig[0], 1.0e-15
  assert_in_delta  5, eig[1], 1.0e-15
  assert_in_delta  1, eig[2], 1.0e-15
 end

 def testEigTriDiag1554
  m = [4, 0, 0,
         13,-7,
             7]
  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta 17.6157731, eig[0], 1.0e-7
  assert_in_delta 4,          eig[1], 1.0e-7
  assert_in_delta 2.38422690, eig[2], 1.0e-7
 end

 def testEigTriDiagSplit1554
  m = [13,0,-4,
          4, 0,
             7]
  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta 15, eig[0], 1.0e-15
  assert_in_delta  5, eig[1], 1.0e-15
  assert_in_delta  4, eig[2], 1.0e-15
 end

 def testEigTriDiagRandom
  m = [  0.22461, 0.43558, 0.12848,
                  0.40385, 0.65227,
                           0.75951]
  d = @gm.triDiag(m)
  assert_in_delta( 0.22461, d[0], 1.0e-5 )
  assert_in_delta( 0.78631, d[1], 1.0e-5 )
  assert_in_delta( 0.37705, d[2], 1.0e-5 )
  e = @gm.triOffDiag(m)
  q0 = @gm.triDiagTransform0(m)
  q1 = @gm.triDiagTransform1(m)
  q2 = @gm.triDiagTransform2(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta(  1.37923, eig[0], 1.0e-5 )
  assert_in_delta(  0.27956, eig[1], 1.0e-5 )
  assert_in_delta( -0.27082, eig[2], 1.0e-5 )
  v0 = @gm.vectTriDiag0(d,e,q0,q1,q2)
  v1 = @gm.vectTriDiag1(d,e,q0,q1,q2)
  v2 = @gm.vectTriDiag2(d,e,q0,q1,q2)
 end

 def reconstruct(e,v0,v1,v2)
  
 end

 def testEigTriDiagMessyDiag
  m = [ 5669.182266666660325,    0.000000000379356,    0.000000000497356,
                              5669.182266666660325,    0.000000000436356,
                                                    5669.182266666660325]
  ans = 5669.18227
  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  q0 = @gm.triDiagTransform0(m)
  q1 = @gm.triDiagTransform1(m)
  q2 = @gm.triDiagTransform2(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta( ans, eig[0], 1.0e-5 )
  assert_in_delta( ans, eig[1], 1.0e-5 )
  assert_in_delta( ans, eig[2], 1.0e-5 )
  v0 = @gm.vectTriDiag0(d,e,q0,q1,q2)
  v1 = @gm.vectTriDiag1(d,e,q0,q1,q2)
  v2 = @gm.vectTriDiag2(d,e,q0,q1,q2)
 end

 def testEigTriDiag100to1
  m = [ 100,  1.0e-5, 1.0e-6,
             1000000, 1.0e-6,
                         100]

  d = @gm.triDiag(m)
  e = @gm.triOffDiag(m)
  eig = @gm.eigTriDiag(d,e)
  assert_in_delta( 1000000, eig[0], 1.0e-5 )
  assert_in_delta(     100, eig[1], 1.0e-5 )
  assert_in_delta(     100, eig[2], 1.0e-5 )
 end

 def testLU3x3eye
  a = [ 1,0,0, 0,1,0, 0,0,1 ]
  lu = @gm.lu3x3(a)
  9.times do |n|
   assert_in_delta( a[n], lu[n], 1.0e-15, n.to_s+" entry bad" )
  end
 end

 def testLU3x3full
  a = [ 3,0,6, 5,8,2, 2,2,8 ]
  g = [ 3,0,2, 5,8,-1, 2,2,6 ]
  lu = @gm.lu3x3(a)
  9.times do |n|
   assert_in_delta( g[n], lu[n], 1.0e-15, n.to_s+" entry bad" )
  end
 end

 def testBackSolve3x3eye
  lu = [ 1,0,0, 0,1,0, 0,0,1 ]
  b = [1,2,3]
  s = @gm.backsolve3x3(lu,b)
  3.times do |n|
   assert_in_delta( b[n], s[n], 1.0e-15, n.to_s+" entry bad" )
  end
 end

 def testMatrixDeterminate
  tol = 1.0e-14
  assert_in_delta( 0.0, @gm.matrixDeterminate([ 1,0,0, 0,1,0, 0,1,0 ]), tol)

  assert_in_delta( 1.0, @gm.matrixDeterminate([ 1,0,0, 0,1,0, 0,0,1 ]), tol)
  assert_in_delta( 1.0, @gm.matrixDeterminate([ 0,1,0, 0,0,1, 1,0,0 ]), tol)
  assert_in_delta( 1.0, @gm.matrixDeterminate([ 0,0,1, 1,0,0, 0,1,0 ]), tol)

  assert_in_delta(-1.0, @gm.matrixDeterminate([ 0,0,1, 0,1,0, 1,0,0 ]), tol)
  assert_in_delta(-1.0, @gm.matrixDeterminate([ 0,1,0, 1,0,0, 0,0,1 ]), tol)
  assert_in_delta(-1.0, @gm.matrixDeterminate([ 1,0,0, 0,0,1, 0,1,0 ]), tol)

  assert_in_delta(-1.0, @gm.matrixDeterminate([-1,0,0, 0,1,0, 0,0,1 ]), tol)
  assert_in_delta(-9.0, @gm.matrixDeterminate([-1,0,0, 0,9,0, 0,0,1 ]), tol)

  assert_in_delta(-9.0, @gm.matrixDeterminate([-1,-2,-2, -4,-3,-4, -4,-4,-3 ]),
                  tol)
 end

 def testFindBaryCornersForRightTriangle
  tol = 1.0e-15
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]

  xyz = xyz0
  bary = [1,0,0]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz) 
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)

  xyz = xyz1
  bary = [0,1,0]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz) 
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)

  xyz = xyz2
  bary = [0,0,1]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz) 
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)
 end

 def testBaryInterpolateForRightTriangle
  tol = 1.0e-15
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]

  xyz = [1.0/2.0,1.0/2.0,0]
  bary = [0,1.0/2.0,1.0/2.0]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz)
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)

  xyz = [1.0/3.0,1.0/3.0,0]
  bary = [1.0/3.0,1.0/3.0,1.0/3.0]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz) 
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)
 end

 def testBaryInterpolateOutsideRightTriangle
  tol = 1.0e-15
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]

  xyz = [1.0,1.0,0.0]
  bary = [-1.0,1.0,1.0]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz) 
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)

  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [1,1,0]

  xyz = [0.0,1.0,0.0]
  bary = [1.0,-1.0,1.0]
  ans = @gm.barycentricCoordinateTri(xyz0,xyz1,xyz2,xyz) 
  assert_in_delta(bary[0],ans[0],tol)
  assert_in_delta(bary[1],ans[1],tol)
  assert_in_delta(bary[2],ans[2],tol)
 end

 def test_gaussianElimination_for_little_2x3
  a = [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ]
  ans = @gm.gaussianElimination(2,3,a)

  tol = 1.0e-14
  assert_in_delta(1.0,ans[0],tol)
  assert_in_delta(0.0,ans[1],tol)
  assert_in_delta(2.0,ans[2],tol)
  assert_in_delta(1.0,ans[3],tol)
  assert_in_delta(3.0,ans[4],tol)
  assert_in_delta(2.0,ans[5],tol)
 end

 def test_gaussianElimination_for_3_tri_case
  a = [ 0.0, 0.5, -0.5,  -0.5, 0.5, 0.0,  1.0, 1.0, 1.0 ]+
      [ -10.0, -15.5, -5.0 ]
  ans = @gm.gaussianElimination(3,4,a)

  tol = 1.0e-14
  assert_in_delta(0.0,ans[1],tol,'terms not elimated')
  assert_in_delta(0.0,ans[2],tol,'terms not elimated')
  assert_in_delta(0.0,ans[5],tol,'terms not elimated')

  assert_in_delta(1.0,ans[0],tol,'terms eye')
  assert_in_delta(1.0,ans[4],tol,'terms eye')
  assert_in_delta(1.0,ans[8],tol,'terms eye')

  assert_in_delta( 2.0,ans[6],tol);
  assert_in_delta(-2.0,ans[7],tol);

  assert_in_delta(-31.0,ans[9],tol);
  assert_in_delta( 20.0,ans[10],tol);
  assert_in_delta(-10.0-1.0/6.0,ans[11],tol);
 end
 
 def test_gaussianBacksolve_little_2x3
  a = [ 1.0, 0.0, 2.0, 1.0, 3.0, 2.0 ]
  ans = @gm.gaussianBacksolve(2,3,a)

  tol = 1.0e-14

  assert_in_delta(1.0,ans[0],tol,'matrix was molested')
  assert_in_delta(0.0,ans[1],tol,'matrix was molested')
  assert_in_delta(2.0,ans[2],tol,'matrix was molested')
  assert_in_delta(1.0,ans[3],tol,'matrix was molested')

  assert_in_delta(-1.0,ans[4],tol)
  assert_in_delta( 2.0,ans[5],tol)
 end

 def test_gaussianElimination_for_3_tri_case
  a = [ 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, -2.0, 1.0, 
       -31.0, 20.0, -10.0-1.0/6.0 ]
  ans = @gm.gaussianBacksolve(3,4,a)

  tol = 1.0e-14
  assert_in_delta(0.0,ans[1],tol,'matrix was molested')
  assert_in_delta(0.0,ans[2],tol,'matrix was molested')
  assert_in_delta(0.0,ans[5],tol,'matrix was molested')

  assert_in_delta(1.0,ans[0],tol,'matrix was molested')
  assert_in_delta(1.0,ans[4],tol,'matrix was molested')
  assert_in_delta(1.0,ans[8],tol,'matrix was molested')

  assert_in_delta( 2.0,ans[6],tol,'matrix was molested');
  assert_in_delta(-2.0,ans[7],tol,'matrix was molested');

  assert_in_delta(-10.0-1.0/3.0,ans[9],tol);
  assert_in_delta(-1.0/3.0,ans[10],tol);
  assert_in_delta(-10.0-1.0/6.0,ans[11],tol);
 end

 def test_impliedMetric_right_tet
  xyz0 = [0,0,0]
  xyz1 = [1,0,0]
  xyz2 = [0,1,0]
  xyz3 = [0,0,1]
  m = @gm.impliedMetric(xyz0,xyz1,xyz2,xyz3)
  tol = 1.0e-14
  assert_in_delta(1.0,m[0],tol)
  assert_in_delta(0.5,m[1],tol)
  assert_in_delta(0.5,m[2],tol)
  assert_in_delta(1.0,m[3],tol)
  assert_in_delta(0.5,m[4],tol)
  assert_in_delta(1.0,m[5],tol)
 end

 def test_impliedMetric_iso_tet
  xyz0 = [ 0.000, 0.000, 0.000 ]
  xyz1 = [ 1.000, 0.000, 0.000 ]
  xyz2 = [ 0.500, 0.866, 0.000 ]
  xyz3 = [ 0.500, 0.289, 0.817 ]
  m = @gm.impliedMetric(xyz0,xyz1,xyz2,xyz3)
  tol = 0.01
  assert_in_delta(1.0,m[0],tol)
  assert_in_delta(0.0,m[1],tol)
  assert_in_delta(0.0,m[2],tol)
  assert_in_delta(1.0,m[3],tol)
  assert_in_delta(0.0,m[4],tol)
  assert_in_delta(1.0,m[5],tol)
 end

end
