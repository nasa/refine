#!/usr/bin/env ruby
#

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'extend'

extend_with('ref_cell')

require 'test/unit'
require 'ref_cell/ref_cell'

class TestAdj < Test::Unit::TestCase

 def set_up
  @ref_cell = Ref_Cell.new(4)
 end
 def setup ; set_up ; end

end
