#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for sampleunit c lib

exit 1 unless system 'ruby makeRubyExtension.rb SampleUnit master_header.h'

require 'test/unit'
require 'SampleUnit/SampleUnit'

class TestSampleUnit < Test::Unit::TestCase

  def testSampleUnit
    s = SampleUnit.new
    assert_equal 3, s.sampleUnit(1,2)
  end

end
