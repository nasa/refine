#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for sampleunit c lib

`rm -rf unittest`
`mkdir unittest`
`cp sampleunit.h sampleunit.c sampleunit_ruby.c unittest`
Dir.chdir "unittest"
`ruby ../extconf.rb SampleUnit`
`make`

require 'test/unit'
require 'SampleUnit'

class TestSampleUnit < Test::Unit::TestCase

  def testSampleUnit
    s = SampleUnit.new
    assert_equal 3, s.sampleUnit(1,2)
  end

end
