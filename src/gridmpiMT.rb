#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for grid c lib

require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('GridMPI').build

require 'test/unit'
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMPI/GridMPI'

class Grid
 include GridMPI
end

class TestGridMPI < Test::Unit::TestCase

 def rightTet
  grid = Grid.new(5,2,0,0)
  grid.addCell( 
	       grid.addNode(0,0,0), 
	       grid.addNode(1,0,0), 
	       grid.addNode(0,1,0), 
	       grid.addNode(0,0,1) )
  grid.identityGlobal
  grid.setAllLocal
  grid
 end 

 def testSplitEdgeAcrossProc
  p1 = rightTet
  p2 = rightTet
  p2.setGhost(0)
  p2.setGhost(1)
  p2.setGhost(2)
  p1.setGhost(3)
 end

end
