#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for gridshape c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('GridShape').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridShape/GridShape'

class Grid
 include GridShape
end

class TestGridShape < Test::Unit::TestCase

 def setup
  @g = Grid.new(0,0,0,0)
 end
 def set_up; setup; end

 def testForStupidity
  assert true
 end

 def testPlotting
  grid=Grid.new(3,0,1,0)
  grid.addNode(0,0,0)
  grid.addNode(0.5,0.6,0.7)
  grid.addNode(-0.8,0.5,-0.9)
  grid.addFace(0,1,2,1)
  grid.plotMinDeterminateAtSurface
 end

 def testFirstOrderLagrangeJacobian0111
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  where = [0.3, 0.3, 0.3]
  
  jac = @g.shapeJacobian1(n0,n1,n2,n3,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)
 end

 def testFirstOrderLagrangeJacobian1211
  n0 = [1, 0, 0]
  n1 = [2, 0, 0]
  n2 = [1, 1, 0]
  n3 = [1, 0, 1]
  where = [0.3, 0.3, 0.3]
  
  jac = @g.shapeJacobian1(n0,n1,n2,n3,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)
 end

 def testFirstOrderLagrangeJacobian0112
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 2]
  where = [0.3, 0.3, 0.3]
  
  jac = @g.shapeJacobian1(n0,n1,n2,n3,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(2.0,jac[8],tol)
 end

 def testFirstOrderLagrangeJacobian1011inverted
  n0 = [1, 0, 0]
  n1 = [0, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  where = [0.3, 0.3, 0.3]
  
  jac = @g.shapeJacobian1(n0,n1,n2,n3,where)
  tol=1.0e-14
  assert_in_delta(-1.0,jac[0],tol)
  assert_in_delta(-1.0,jac[1],tol)
  assert_in_delta(-1.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)
 end

 def testFirstOrderLagrangeJacobian1011twoFlip
  n0 = [1, 0, 0]
  n1 = [0, 0, 0]
  n2 = [0, 0, 1]
  n3 = [0, 1, 0]
  where = [0.3, 0.3, 0.3]
  
  jac = @g.shapeJacobian1(n0,n1,n2,n3,where)
  tol=1.0e-14
  assert_in_delta(-1.0,jac[0],tol)
  assert_in_delta(-1.0,jac[1],tol)
  assert_in_delta(-1.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(0.0,jac[4],tol)
  assert_in_delta(1.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(1.0,jac[7],tol)
  assert_in_delta(0.0,jac[8],tol)
 end

 def testSecondOrderLagrangeJacobian0111
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  h=0.5
  e01 = [h, 0, 0]
  e02 = [0, h, 0]
  e03 = [0, 0, h]
  e12 = [h, h, 0]
  e13 = [h, 0, h]
  e23 = [0, h, h]
  
  where = [0.3, 0.3, 0.3]
  
  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)

  where = [0.0, 0.0, 0.0]
  
  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)

  where = [0.0, 1.0, 0.0]
  
  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)
 end

 def testSecondOrderLagrangeJacobianIceCreamCone
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  h=0.5
  e01 = [h, 0, 0]
  e02 = [0, h, 0]
  e03 = [0, 0, h]
  e12 = [1, 1, 0]
  e13 = [1, 0, 1]
  e23 = [0, 1, 1]
  
  where = [0.0, 0.0, 0.0]
  
  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)

  where = [0.5, 0.5, 0.0]
  
  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(2.0,jac[0],tol)
  assert_in_delta(1.0,jac[1],tol)
  assert_in_delta(1.0,jac[2],tol)
  assert_in_delta(1.0,jac[3],tol)
  assert_in_delta(2.0,jac[4],tol)
  assert_in_delta(1.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(3.0,jac[8],tol)
 end

 def testSecondOrderLagrangeJacobianPartialIceCreamCone
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  h=0.5
  e01 = [h, 0, 0]
  e02 = [0, h, 0]
  e03 = [0, 0, h]
  e12 = [h, h, 0]
  e13 = [1, 0, 1]
  e23 = [0, 1, 1]
  
  where = [0.5, 0.5, 0.0]

  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(1.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(1.0,jac[2],tol)
  assert_in_delta(0.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(1.0,jac[5],tol)
  assert_in_delta(0.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(3.0,jac[8],tol)
 end

 def testSecondOrderLagrangeJacobianInvertedEdge01
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  h=0.5
  e01 = [1, 1, 1]
  e02 = [0, h, 0]
  e03 = [0, 0, h]
  e12 = [h, h, 0]
  e13 = [h, 0, h]
  e23 = [0, h, h]
  
  where = [0.0, 0.0, 0.0]

  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(3.0,jac[0],tol)
  assert_in_delta(0.0,jac[1],tol)
  assert_in_delta(0.0,jac[2],tol)
  assert_in_delta(4.0,jac[3],tol)
  assert_in_delta(1.0,jac[4],tol)
  assert_in_delta(0.0,jac[5],tol)
  assert_in_delta(4.0,jac[6],tol)
  assert_in_delta(0.0,jac[7],tol)
  assert_in_delta(1.0,jac[8],tol)

  where = [1.0, 0.0, 0.0]

  jac = @g.shapeJacobian2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(-1.0,jac[0],tol)
  assert_in_delta(-2.0,jac[1],tol)
  assert_in_delta(-2.0,jac[2],tol)
  assert_in_delta(-4.0,jac[3],tol)
  assert_in_delta(-3.0,jac[4],tol)
  assert_in_delta(-4.0,jac[5],tol)
  assert_in_delta(-4.0,jac[6],tol)
  assert_in_delta(-4.0,jac[7],tol)
  assert_in_delta(-3.0,jac[8],tol)
 end

 def testSecondOrderLagrangeJacobianDeterminateInvertedEdge01
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  h=0.5
  e01 = [1, 1, 1]
  e02 = [0, h, 0]
  e03 = [0, 0, h]
  e12 = [h, h, 0]
  e13 = [h, 0, h]
  e23 = [0, h, h]
  
  where = [0.0, 0.0, 0.0]

  det = @g.shapeJacobianDet2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta( 3.0,det,tol)

  where = [1.0, 0.0, 0.0]

  det = @g.shapeJacobianDet2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta(-9.0,det,tol)
 end

 def testSecondOrderLagrangeJacobianDeterminateDerivative
  n0 = [0, 0, 0]
  n1 = [1, 0, 0]
  n2 = [0, 1, 0]
  n3 = [0, 0, 1]
  h=0.5
  e01 = [h, 0, 0]
  e02 = [0, h, 0]
  e03 = [0, 0, h]
  e12 = [h, h, 0]
  e13 = [h, 0, h]
  e23 = [0, h, h]
  
  where = [0.0, 0.0, 0.0]

  deriv = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  tol=1.0e-14
  assert_in_delta( 1.0, deriv[0], tol)

  delta = 1.0e-8
  tol=delta

  n0 = [delta, 0, 0]
  p = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  n0 = [-delta, 0, 0]
  m = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  fd = (p[0]-m[0])/delta*0.5

  assert_in_delta( fd, deriv[1], tol)
  assert_in_delta(-3.0, deriv[1], tol)

  n0 = [0, delta, 0]
  p = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  n0 = [0, -delta, 0]
  m = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  fd = (p[0]-m[0])/delta*0.5

  assert_in_delta( fd, deriv[2], tol)
  assert_in_delta(-3.0, deriv[2], tol)

  n0 = [0, 0, delta]
  p = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  n0 = [0, 0, -delta]
  m = @g.shapeJacobianDetDeriv2(n0,n1,n2,n3,e01,e02,e03,e12,e13,e23,where)
  fd = (p[0]-m[0])/delta*0.5

  assert_in_delta( fd, deriv[3], tol)
  assert_in_delta(-3.0, deriv[3], tol)

 end

 def testSecondOrderLagrangeJacobianDeterminateMin
  @g.addNode(0, 0, 0)
  @g.addNode(1, 0, 0)
  @g.addNode(0, 1, 0)
  @g.addNode(0, 0, 1)
  nodes = [0,1,2,3]

  tol=1.0e-14
  assert_in_delta(1.0,@g.minCellJacDet2(nodes),tol)
  @g.setNodeXYZ(3, [0, 0, 2])
  assert_in_delta(2.0,@g.minCellJacDet2(nodes),tol)
 end


end
