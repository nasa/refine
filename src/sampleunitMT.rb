#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for sampleunit c lib

system 'ruby extconf.rb SampleUnit'
system 'make'

require 'test/unit'
require 'SampleUnit'

class TestSampleUnit < Test::Unit::TestCase

  def testSampleUnit
    s = SampleUnit.new
    assert_equal 3, s.sampleUnit(1,2)
  end

end
