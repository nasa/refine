
require 'mkmf'

$objs = %w{sampleunit.o sampleunit_ruby.o}
create_makefile("SampleUnit")
