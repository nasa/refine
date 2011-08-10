#!/usr/bin/env ruby
#

Dir.chdir ENV['srcdir'] if ENV['srcdir']

require 'extend'

extend_with('ref_cell')

require 'test/unit'
require 'ref_cell/ref_cell'

class Test_Ref_Cell < Test::Unit::TestCase

 def successful?(call)
   assert_equal 0, call
 end

 def set_up
  @ref_cell = Ref_Cell.new(4)
 end
 def setup ; set_up ; end

 def test_starts_with_zero_cells
   assert_equal 0, @ref_cell.n
 end

 def test_add_cell
   assert_equal 0, @ref_cell.add([0,1,2,3])
 end

 def test_add_cell_past_max_reallocs
   max = @ref_cell.max
   (max+1).times { |cell| assert_equal cell, @ref_cell.add([0,1,2,3]) }
   assert( @ref_cell.max > max )
 end

 def test_cell_nodes
   assert_equal nil, @ref_cell.nodes(0)
   assert_equal nil, @ref_cell.nodes(1000000000)
   @ref_cell.add([0,1,2,3])
   assert_equal [0,1,2,3], @ref_cell.nodes(0)
 end

end
