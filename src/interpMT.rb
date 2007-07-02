#!/usr/bin/env ruby
#
# Mobility test for the interp c lib
#
# $Id: interpMT.rb 2979 2005-08-05 13:07:18Z mikepark $

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Interp').build

require 'test/unit'
require 'Interp/Interp'

class TestInterp < Test::Unit::TestCase

 def set_up
  @i0 = Interp.new(0,1)
 end
 def setup ; set_up ; end

 def test_create
  assert_equal 0, @i0.functionId
  assert_equal 1, @i0.order
 end

 def test_function
  tol = 1.0e-14
  assert_in_delta( 0.0, @i0.function([0.0,0.0,0.0]), tol )
 end

 def test_error
  tol = 1.0e-12
  truth = 121.6
  xyz0 = [0.0,0.0,0.0]
  xyz1 = [1.0,0.0,0.0]
  xyz2 = [0.0,1.0,0.0]
  xyz3 = [0.0,0.0,1.0]
  assert_in_delta( truth, @i0.error(xyz0,xyz1,xyz2,xyz3), tol )
  xyz0 = [1.0,0.0,0.0]
  xyz1 = [2.0,0.0,0.0]
  xyz2 = [1.0,1.0,0.0]
  xyz3 = [1.0,0.0,1.0]
  assert_in_delta( truth, @i0.error(xyz0,xyz1,xyz2,xyz3), tol )
 end

 def test_metric
  tol = 1.0e-14
  res = @i0.metric([0.0,0.0,0.0])
  ans = [16.0, 0.0, 0.0, 64.0, 0.0, 128.0]
  6.times { |i| assert_in_delta( ans[i], res[i], tol, 'element '+i.to_s ) }
 end

end
