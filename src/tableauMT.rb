#!/usr/bin/env ruby
#
# Mobility test for the tableau c lib
#


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

 def test_initial_tableau_for_three_right_triangles
  tableau = Tableau.new(3,3)
  const_mat = [ 0.0, -0.5, 1.0,
                0.5,  0.5, 1.0,
               -0.5,  0.0, 1.0]
  assert_equal tableau, tableau.constraintMatrix(const_mat)
  assert_equal tableau, tableau.constraint([ 0.0, 0.0, 1.0])
  assert_equal tableau, tableau.cost([ 0.0, 0.5, 0.0])
  assert_equal tableau, tableau.init
  t = tableau.tableau
  m = 4
  tol = TOL
  j = 0; truth = [  0.0,  0.0,  0.0,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 1; truth = [  0.0,  0.0, -0.5,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 2; truth = [  0.5,  0.5,  0.5,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 3; truth = [  0.0,-0.5,  0.0,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 4; truth = [  0.0,  1.0,  0.0,  0.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 5; truth = [  0.0,  0.0,  1.0,  0.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 6; truth = [  0.0,  0.0,  0.0,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
 end

 def test_initial_tableau_for_three_slanted_triangles
  tableau = Tableau.new(3,3)
  const_mat = [-0.15, -0.70, 1.0,
                0.50,  0.45, 1.0,
               -0.35,  0.25, 1.0]
  assert_equal tableau, tableau.constraintMatrix(const_mat)
  assert_equal tableau, tableau.constraint([ 0.0, 0.0, 1.0])
  assert_equal tableau, tableau.cost([ -0.095, 0.505, 0.155])
  assert_equal tableau, tableau.init
  t = tableau.tableau
  m = 4
  tol = 1.0e-5
  j = 0; truth = [ 0.0,  0.0,  0.0,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 1; truth = [ -0.095, -0.15, -0.7,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 2; truth = [  0.505,  0.5,  0.45,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 3; truth = [ 0.155,  -0.35, 0.25,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 4; truth = [  0.0,  1.0,  0.0,  0.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 5; truth = [  0.0,  0.0,  1.0,  0.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
  j = 6; truth = [  0.0,  0.0,  0.0,  1.0 ]
  m.times { |i| assert_in_delta truth[i], t[j][i], tol }
 end

 def test_largest_pivot_for_three_skew_triangles
  tableau = Tableau.new(3,3)
  const_mat = [-0.15, -0.70, 1.0,
                0.50,  0.45, 1.0,
               -0.35,  0.25, 1.0]
  assert_equal tableau, tableau.constraintMatrix(const_mat)
  assert_equal tableau, tableau.constraint([ 0.0, 0.0, 1.0])
  assert_equal tableau, tableau.cost([ -0.095, 0.505, 0.155])
  assert_equal tableau, tableau.init
  assert_equal [3, 1], tableau.largestPivot
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
  assert_equal [0, 1, 2], tableau.basis.sort
 end

end
