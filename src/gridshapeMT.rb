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
  puts jac
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

end
