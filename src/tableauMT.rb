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
  tableau = Tableau.new(3,3)
  assert_equal [3, 4, 5], tableau.basis
  tableau = Tableau.new(4,5)
  assert_equal [5, 6, 7, 8], tableau.basis
 end

 def test_find_right_basis_for_three_triangles
  tableau = Tableau.new(3,3)
  const_mat = [ 0.0, -0.5, 1.0,
                0.5,  0.5, 1.0,
               -0.5,  0.0, 1.0]
  assert_equal tableau, tableau.constraintMatrix(const_mat)
  assert_equal tableau, tableau.constraint([ 0.0, 0.0, 1.0])
  assert_equal tableau, tableau.cost([ 0.0, 0.5, 0.0])
  assert_equal tableau, tableau.solve
  assert_equal [0, 1, 2], tableau.basis
 end

 

end
