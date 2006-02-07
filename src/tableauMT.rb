#!/usr/bin/env ruby
#
# Mobility test for the tableau c lib
#
# $Id$

valgrinding = false
if (valgrinding)
 GC.disable # for runnig valgridn quietly
else
 Dir.chdir ENV['srcdir'] if ENV['srcdir']
 require 'RubyExtensionBuilder'
 RubyExtensionBuilder.new('Tableau').build
end

require 'test/unit'
require 'Tableau/Tableau'

class TestTableau < Test::Unit::TestCase

 EMPTY = -1
 TOL = 1.0e-15

 def test_create_initializes_problem_size
  tableau = Tableau.new(1,2)
  assert_equal 1, tableau.constraints
  assert_equal 2, tableau.dimension
 end

 def test_initial_basis_is_after_given_dimension
  tableau = Tableau.new(2,3)
  assert_equal [3, 4], tableau.basis
 end

end
