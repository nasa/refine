#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for sampleunit c lib

exit unless system 'ruby makeRubyExtension.rb SampleUnit'

require 'test/unit'
require 'SampleUnit/SampleUnit'

class TestSampleUnit < Test::Unit::TestCase

  def testSampleUnit
    s = SampleUnit.new
    assert_equal 3, s.sampleUnit(1,2)
  end

end
