#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for sampleunit c lib

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'mkmf'

ext = 'SampleUnit'
`mkdir -p SampleUnit`
Dir.chdir ext
$objs = %w[sampleunit.o sampleunit_ruby.o]
create_makefile(ext,'..')
exit 1 unless system "make --quiet --no-print-directory"
Dir.chdir '..'

require 'test/unit'
require 'SampleUnit/SampleUnit'

class TestSampleUnit < Test::Unit::TestCase

  def testSampleUnit
    s = SampleUnit.new
    assert_equal 3, s.sampleUnit(1,2)
  end

end
