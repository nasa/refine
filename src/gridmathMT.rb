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

class TestGridMath < Test::Unit::TestCase

 def set_up
  @gm = GridMath.new
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

 def testRotateDirectionEndPoints
  v0 = [1,0,0]
  v1 = [0,1,0]
  axle = [0,0,1]
  assert_equal v0, @gm.rotateDirection(v0,v1,axle,0)
  assert_equal v1, @gm.rotateDirection(v0,v1,axle,1)
  negaxle = [0,0,-1]
  assert_equal v1, @gm.rotateDirection(v1,v0,negaxle,0)
  assert_equal v0, @gm.rotateDirection(v1,v0,negaxle,1)
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
  assert_in_delta( 0.57735027, v0[0],1.0e-7)
  assert_in_delta( 0.57735027, v0[1],1.0e-7)
  assert_in_delta( 0.57735027, v0[2],1.0e-7)
  assert_in_delta(-0.81649658, v1[0],1.0e-7)
  assert_in_delta( 0.40824829, v1[1],1.0e-7)
  assert_in_delta( 0.40824829, v1[2],1.0e-7)
  assert_in_delta( 0.0,        v2[0],1.0e-7)
  assert_in_delta( invsqrt2,   v2[1],1.0e-7)
  assert_in_delta(-invsqrt2,   v2[2],1.0e-7)
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

end
