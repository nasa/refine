#!/usr/bin/env ruby
#
# $Id$
#
# Mobility test for layer c lib

#GC.disable # turn off GC if layer does not mark grid

Dir.chdir ENV['srcdir'] if ENV['srcdir']
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Layer').build

require 'test/unit'
require 'Adj/Adj'
require 'Near/Near'
require 'Intersect/Intersect'
require 'Line/Line'
require 'Grid/Grid'
require 'GridMath/GridMath'
require 'GridMetric/GridMetric'
require 'GridCAD/GridCAD'
require 'Layer/Layer'

class Grid
 include GridMetric
 include GridCAD
end

class TestLayer < Test::Unit::TestCase

 EMPTY = (-1)

 def testInit
  assert_not_nil  grid = Grid.new(2,0,0,0)
  assert_not_nil  layer = Layer.new(grid)
  assert_equal 2, layer.maxnode
  assert_equal 0, layer.ntriangle
  assert_equal 0, layer.maxtriangle
  assert_equal 0, layer.nnormal
  assert_equal 0, layer.maxnormal
 end

 def testInitGC
  assert_not_nil  grid = Grid.new(2,0,0,0)
  assert_not_nil  layer = Layer.new(grid)
  assert_equal 2, layer.maxnode
  assert_not_nil  grid = String.new("hello")
  assert_nil      GC.start
  assert_equal 2, layer.maxnode
 end

 def testInsertTriangleIntoLayer
  layer = Layer.new(Grid.new(4,0,2,0))
  assert_nil               layer.triangle(0)
  assert_nil               layer.triangleNormals(0)
  assert_equal layer,      layer.addTriangle(0,1,2)
  assert_equal 1,          layer.ntriangle
  assert                   layer.maxtriangle>=1
  assert_nil               layer.triangle(-1)
  assert_equal [0,1,2],    layer.triangle(0)
  assert_equal [0,1,2],    layer.triangleNormals(0)
  assert_equal 0,          layer.constrainedSide(0,0)
  assert_equal 0,          layer.constrainedSide(0,1)
  assert_equal 0,          layer.constrainedSide(0,2)
  assert_equal 0,          layer.parentGeomEdge(0,0)
  assert_equal 0,          layer.parentGeomEdge(0,1)
  assert_equal 0,          layer.parentGeomEdge(0,2)
  assert_nil               layer.triangle(1)
 end

 def testPopulateAdvancingFront
  assert_not_nil        grid = Grid.new(4,0,2,0)
  0.upto(3) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  grid.addFace(0,1,3,2)
  assert_equal 2,       grid.nface
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.populateAdvancingFront([1,2])
  assert_equal 2,       layer.ntriangle
  assert_equal [0,1,2], layer.triangle(0)
  assert_equal [0,1,3], layer.triangle(1)
 end

 def testRememberTriangleParentGeomFaces
  layer = Layer.new(Grid.new(0,0,0,0))
  assert_equal false,   layer.parentGeomFace(5)
  assert_equal false,   layer.parentGeomFace(12)
  assert_equal layer,   layer.addParentGeomFace(12)
  assert_equal false,   layer.parentGeomFace(5)
  assert_equal true,    layer.parentGeomFace(12)
  assert_equal layer,   layer.addParentGeomFace(5)
  assert_equal true,    layer.parentGeomFace(5)
  assert_equal true,    layer.parentGeomFace(12)
  assert_nil            layer.addParentGeomFace(5)
 end

 def testAddNormalDynamicAllocates
  layer = Layer.new(Grid.new(235,0,0,0))
  assert_equal(-1,         layer.addNormal(-1) )
  assert_equal(-1,         layer.addNormal(layer.maxnode) )
  assert_equal 0,          layer.nnormal
  assert_equal 0,          layer.addNormal(33)
  assert_equal 1,          layer.addNormal(234)
  assert_equal 2,          layer.nnormal
  assert                   layer.maxnormal>=layer.nnormal
 end

 def testWeGetUniqueNormals
  layer = Layer.new(Grid.new(42,0,0,0))
  assert_equal 0, layer.uniqueNormalId(40)
  assert_equal 0, layer.uniqueNormalId(40)
  assert_equal 1, layer.uniqueNormalId(24)
  assert_equal 0, layer.uniqueNormalId(40)
 end

 def testUniqueNormalSetsANormalsRootToGlobalNode
  layer = Layer.new(Grid.new(42,0,0,0))
  assert_equal 0,  layer.uniqueNormalId(40)
  assert_equal 1,  layer.uniqueNormalId(24)
  assert_equal 40, layer.normalRoot(0)  
  assert_equal 24, layer.normalRoot(1)  
 end

 def testMakeNormals
  assert_not_nil        grid = Grid.new(5,0,2,0)
  0.upto(4) {grid.addNode(0,0,0)}
  grid.addFace(1,2,3,1)
  grid.addFace(1,2,4,2)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.populateAdvancingFront([1,2])
  assert_equal 4,       layer.nnormal
  assert_equal [0,1,2], layer.triangleNormals(0)
  assert_equal [0,1,3], layer.triangleNormals(1)
  assert_equal(-1,       layer.normalRoot(-1))
  assert_equal(-1,       layer.normalRoot(100))
  assert_equal 1,       layer.normalRoot(0)
  assert_equal 2,       layer.normalRoot(1)
  assert_equal 3,       layer.normalRoot(2)
  assert_equal 4,       layer.normalRoot(3)
 end

 def testConstrainNormals
  assert_not_nil        grid = Grid.new(6,0,3,0)
  0.upto(3) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  grid.addFace(0,1,3,2)
  grid.addFace(0,3,4,773)
  assert_not_nil        layer = Layer.new(grid)
  assert_nil            layer.constrainNormal(3)
  assert_equal 0,       layer.constrained(0)
  assert_equal layer,   layer.populateAdvancingFront([1,2])
  assert_equal 2,       layer.ntriangle
  assert_equal 4,       layer.nnormal
  assert_equal 0,       layer.constrained(0)
  assert_equal 0,       layer.constrained(1)
  assert_equal 0,       layer.constrained(2)
  assert_equal 0,       layer.constrained(3)
  assert_nil            layer.constrainNormal(0)
  assert_equal layer,   layer.constrainNormal(773)
  assert_equal 773,     layer.constrained(0)
  assert_equal 0,       layer.constrained(1)
  assert_equal 0,       layer.constrained(2)
  assert_equal 773,     layer.constrained(3)

  assert_equal 0,       layer.nConstrainedSides(0)
  assert_equal 1,       layer.nConstrainedSides(773)
 end

 def testRememberConstrainingFacesAndEdges
  assert_not_nil        grid = Grid.new(3,0,1,0)
  0.upto(2) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.populateAdvancingFront([1])
  assert_nil            layer.constrainNormal(0)
  assert_equal layer,   layer.constrainNormal(2)
  assert_equal layer,   layer.constrainNormal(-1)
  assert_equal false,   layer.constrainingGeometry(-2)
  assert_equal true,    layer.constrainingGeometry(-1)
  assert_equal false,   layer.constrainingGeometry(0)
  assert_equal false,   layer.constrainingGeometry(1)
  assert_equal true,    layer.constrainingGeometry(2)
  assert_equal false,   layer.constrainingGeometry(3)
 end

 def testConstrainNormalForEdge
  assert_not_nil        grid = Grid.new(10,10,10,10)
  0.upto(3) {grid.addNode(0,0,0)}
  grid.addFace(1,2,3,1)
  grid.addFace(1,2,0,2)
  assert_equal grid,    grid.setNGeomNode(2)
  assert_equal grid,    grid.setNGeomEdge(1)
  assert_equal grid,    grid.addGeomEdge(1,0,2)
  grid.addEdge(0,1,1,0.0,1.0)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.populateAdvancingFront([1])
  assert_equal 3,       layer.nnormal
  assert_equal layer,   layer.constrainNormal(-1)
  assert_equal layer,   layer.constrainNormal(2)
  assert_equal(-1,      layer.constrained(0))
  assert_equal 2,       layer.constrained(1)
  assert_equal 0,       layer.constrained(2)
 end

 def testConstrainTriangleSide
  assert_not_nil        grid = Grid.new(6,0,3,0)
  0.upto(3) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal 0,       layer.constrainedSide(0,0)
  assert_equal layer,   layer.populateAdvancingFront([1])
  assert_equal 0,       layer.constrainedSide(0,0)
  assert_equal 0,       layer.constrainedSide(0,1)
  assert_equal 0,       layer.constrainedSide(0,2)
  assert_equal 0,       layer.constrainedSide(0,3)
  assert_equal 0,       layer.constrainedSide(1,0)
  assert_nil   layer.constrainTriangleSide(-1,0,0)
  assert_nil   layer.constrainTriangleSide(layer.nnormal,0,0)
  assert_nil   layer.constrainTriangleSide(0,-1,0)
  assert_nil   layer.constrainTriangleSide(0,layer.nnormal,0)
  assert_equal 0,       layer.constrainedSide(0,0)
  assert_equal 0,       layer.constrainedSide(0,1)
  assert_equal 0,       layer.constrainedSide(0,2)
  assert_equal layer,   layer.constrainTriangleSide(0,1,1)
  assert_equal layer,   layer.constrainTriangleSide(1,2,2)
  assert_equal layer,   layer.constrainTriangleSide(0,2,3)
  assert_equal 1,       layer.constrainedSide(0,0)
  assert_equal 2,       layer.constrainedSide(0,1)
  assert_equal 3,       layer.constrainedSide(0,2)
 end

 def testConstrainTriangleSideWithBCFace
  assert_not_nil        grid = Grid.new(6,0,3,0)
  0.upto(3) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  grid.addFace(1,2,3,2)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal layer,   layer.populateAdvancingFront([1])
  assert_equal 0,       layer.constrainedSide(0,0)
  assert_equal 0,       layer.constrainedSide(0,1)
  assert_equal 0,       layer.constrainedSide(0,2)
  assert_equal layer,   layer.constrainNormal(2)
  assert_equal 0,       layer.constrainedSide(0,0)
  assert_equal 2,       layer.constrainedSide(0,1)
  assert_equal 0,       layer.constrainedSide(0,2)
 end

 def testRememberTriangleEdgeParent
  assert_not_nil        grid = Grid.new(6,0,3,0)
  0.upto(2) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal 0,       layer.parentGeomEdge(0,0)
  assert_equal layer,   layer.populateAdvancingFront([1])
  assert_equal 0,       layer.parentGeomEdge(0,0)
  assert_equal 0,       layer.parentGeomEdge(0,1)
  assert_equal 0,       layer.parentGeomEdge(0,2)
  assert_equal 0,       layer.parentGeomEdge(0,3)
  assert_equal 0,       layer.parentGeomEdge(1,0)
  assert_nil            layer.setParentGeomEdge(-1,0,0)
  assert_nil            layer.setParentGeomEdge(layer.nnormal,0,0)
  assert_nil            layer.setParentGeomEdge(0,-1,0)
  assert_nil            layer.setParentGeomEdge(0,layer.nnormal,0)
  assert_equal 0,       layer.parentGeomEdge(0,0)
  assert_equal 0,       layer.parentGeomEdge(0,1)
  assert_equal 0,       layer.parentGeomEdge(0,2)
  assert_equal 0,       layer.nParentGeomEdgeSegments(0)
  assert_equal 0,       layer.nParentGeomEdgeSegments(1)
  assert_equal layer,   layer.setParentGeomEdge(0,1,1)
  assert_equal layer,   layer.setParentGeomEdge(1,2,2)
  assert_equal layer,   layer.setParentGeomEdge(0,2,3)
  assert_equal 1,       layer.parentGeomEdge(0,0)
  assert_equal 2,       layer.parentGeomEdge(0,1)
  assert_equal 3,       layer.parentGeomEdge(0,2)
  assert_equal 0,       layer.nParentGeomEdgeSegments(0)
  assert_equal 1,       layer.nParentGeomEdgeSegments(1)
 end

 def testFindParentalEdgesForTriangle
  assert_not_nil        grid = Grid.new(10,10,10,10)
  0.upto(2) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  grid.addEdge(0,1,1,0,1)
  grid.addEdge(0,4,1,0,1)
  assert_not_nil        layer = Layer.new(grid)
  assert_nil            layer.findParentGeomEdges
  assert_equal layer,   layer.populateAdvancingFront([1])
  assert_equal 0,       layer.parentGeomEdge(0,0)
  assert_equal 0,       layer.parentGeomEdge(0,1)
  assert_equal 0,       layer.parentGeomEdge(0,2)
  assert_equal 0,       layer.nParentGeomEdgeSegments(1)
  assert_equal layer,   layer.findParentGeomEdges
  assert_equal 1,       layer.parentGeomEdge(0,0)
  assert_equal 0,       layer.parentGeomEdge(0,1)
  assert_equal 0,       layer.parentGeomEdge(0,2)
  assert_equal 1,       layer.nParentGeomEdgeSegments(1)
 end

 def testNormalTriangleNeighbors
  assert_not_nil        grid = Grid.new(4,0,3,0)
  0.upto(3) {grid.addNode(0,0,0)}
  grid.addFace(0,1,2,1)
  grid.addFace(0,1,3,2)
  grid.addFace(1,2,3,3)
  assert_not_nil        layer = Layer.new(grid)
  assert_equal 0,       layer.normalDeg(0)
  assert_nil            layer.normalTriangles(0)
  assert_equal layer,   layer.populateAdvancingFront([1,2,3])
  assert_equal [0,1,2], layer.triangleNormals(0)
  assert_equal [0,1,3], layer.triangleNormals(1)
  assert_equal [1,2,3], layer.triangleNormals(2)
  assert_equal 2,       layer.normalDeg(0)
  assert_equal 3,       layer.normalDeg(1)
  assert_equal [1,0],   layer.normalTriangles(0)
  assert_equal [2,1,0], layer.normalTriangles(1)
  assert_equal [2,0],   layer.normalTriangles(2)
  assert_equal [2,1],   layer.normalTriangles(3)
 end

 def testNormalDirection1Face
  assert_not_nil          grid = Grid.new(3,0,1,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_nil              layer.triangleDirection(-1)
  assert_nil              layer.triangleDirection(0)
  assert_nil              layer.normalDirection(0)
  assert_equal layer,     layer.populateAdvancingFront([1])
  direction = [0.0,0.0,1.0]
  assert_equal direction, layer.triangleDirection(0)
  assert_equal direction, layer.normalDirection(0)
  assert_equal direction, layer.normalDirection(1)
  assert_equal direction, layer.normalDirection(2)
 end

 def testNormalDirection2Face
  assert_not_nil          grid = Grid.new(4,0,2,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addFace(0,1,2,1)
  grid.addFace(0,3,1,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal [0.0,0.0,1.0], layer.triangleDirection(0)
  assert_equal [0.0,1.0,0.0], layer.triangleDirection(1)
  halfSqrt2 = 0.5 * Math::sqrt(2)
  direction = [0.0, halfSqrt2, halfSqrt2]
  tol = 1.0e-14
  assert_in_delta direction[0], layer.normalDirection(0)[0], tol
  assert_in_delta direction[1], layer.normalDirection(0)[1], tol
  assert_in_delta direction[2], layer.normalDirection(0)[2], tol
  assert_in_delta direction[0], layer.normalDirection(1)[0], tol
  assert_in_delta direction[1], layer.normalDirection(1)[1], tol
  assert_in_delta direction[2], layer.normalDirection(1)[2], tol
  assert_equal [0.0,0.0,1.0], layer.normalDirection(2)
  assert_equal [0.0,1.0,0.0], layer.normalDirection(3)
 end

 def testNormalDirectionVisible
  assert_not_nil          grid = Grid.new(5,0,3,0)
  assert_equal 0,         grid.addNode( 0, 0, 0)
  assert_equal 1,         grid.addNode( 1, 0, 0)
  assert_equal 2,         grid.addNode(-1, 0, 0)
  assert_equal 3,         grid.addNode( 0, 1,-0.1)
  assert_equal 4,         grid.addNode( 0, 1, 0.1)
  grid.addFace(0,3,1,1)
  grid.addFace(0,1,4,1)
  grid.addFace(0,4,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  norm = layer.normalDirection(0)
  assert_in_delta(  0.000, norm[0], 1e-3)
  assert_in_delta( -0.287, norm[1], 1e-3)
  assert_in_delta(  0.957, norm[2], 1e-3)
  assert_equal layer,     layer.visibleNormals
  norm = layer.normalDirection(0)
  assert_in_delta(  0.000, norm[0], 1e-3)
  assert_in_delta( -1.000, norm[1], 1e-4)
  assert_in_delta(  0.000, norm[2], 1e-2)
 end

 def testFreezeNormalDirectionFrozenState
  grid = Grid.new(3,0,1,0)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal false,     layer.normalDirectionFrozen(-1)
  assert_equal false,     layer.normalDirectionFrozen(0)
  assert_equal false,     layer.normalDirectionFrozen(1)
  assert_equal false,     layer.normalDirectionFrozen(5)

  assert_nil              layer.normalDirectionFreeze(-1)
  assert_equal layer,     layer.normalDirectionFreeze(1)
  assert_nil              layer.normalDirectionFreeze(5)

  assert_equal false,     layer.normalDirectionFrozen(0)
  assert_equal true,      layer.normalDirectionFrozen(1)
 end

 def testNormalDirectionFeasible
  assert_not_nil          grid = Grid.new(5,0,3,0)
  assert_equal 0,         grid.addNode( 0, 0, 0)
  assert_equal 1,         grid.addNode( 1, 0, 0)
  assert_equal 2,         grid.addNode(-1, 0, 0)
  assert_equal 3,         grid.addNode( 0, 1,-0.1)
  assert_equal 4,         grid.addNode( 0, 1, 0.1)
  grid.addFace(0,3,1,1)
  grid.addFace(0,1,4,1)
  grid.addFace(0,4,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  norm = layer.normalDirection(0)
  assert_in_delta(  0.000, norm[0], 1e-3)
  assert_in_delta( -0.287, norm[1], 1e-3)
  assert_in_delta(  0.957, norm[2], 1e-3)
  assert_equal layer,     layer.feasibleNormals
  norm = layer.normalDirection(0)
  assert_in_delta(  0.000, norm[0], 1e-8)
  assert_in_delta( -1.000, norm[1], 1e-2)
  assert_in_delta(  0.000, norm[2], 1e-1)
 end

 def testProjectNormalToEdge
  assert_not_nil          grid = Grid.new(9,9,9,9)
  assert_equal 0,         grid.addNode( 0, 0, 0)
  assert_equal 1,         grid.addNode( 0, 1, 0)
  assert_equal 2,         grid.addNode( 1, 0, 1)
  assert_equal 3,         grid.addNode( 1, 0, 0)
  grid.addFace(0,1,2,22)
  grid.addEdge(0,3,1,0,1)
  assert_equal grid,      grid.setNGeomEdge(1)
  assert_equal grid,      grid.addGeomEdge(1,0,3)
  assert_not_nil          layer = Layer.new(grid).populateAdvancingFront([22])
  assert_equal layer,     layer.constrainNormal(-1)
  norm = layer.normalDirection(0)
  assert_in_delta(  0.707, norm[0], 1e-3)
  assert_in_delta(  0.000, norm[1], 1e-3)
  assert_in_delta( -0.707, norm[2], 1e-3)
  assert_equal layer, layer.projectNormalsToConstraints
  norm = layer.normalDirection(0)
  assert_in_delta( 1.000, norm[0], 1e-3)
  assert_in_delta( 0.000, norm[1], 1e-3)
  assert_in_delta( 0.000, norm[2], 1e-3)
 end

 def testProjectNormalToFace
  assert_not_nil          grid = Grid.new(9,9,9,9)
  assert_equal 0,         grid.addNode( 0, 0, 0)
  assert_equal 1,         grid.addNode( 0, 1, 0)
  assert_equal 2,         grid.addNode( 1, 0, 1)
  assert_equal 3,         grid.addNode( 1, 0, 0)
  grid.addFace(0,1,2,22)
  grid.addFace(0,3,1,1)
  assert_not_nil          layer = Layer.new(grid).populateAdvancingFront([22])
  assert_equal layer,     layer.constrainNormal(1)
  norm = layer.normalDirection(0)
  assert_in_delta(  0.707, norm[0], 1e-3)
  assert_in_delta(  0.000, norm[1], 1e-3)
  assert_in_delta( -0.707, norm[2], 1e-3)
  assert_equal layer, layer.projectNormalsToConstraints
  norm = layer.normalDirection(0)
  assert_in_delta( 1.000, norm[0], 1e-3)
  assert_in_delta( 0.000, norm[1], 1e-3)
  assert_in_delta( 0.000, norm[2], 1e-3)
 end

 def testAdvanceLayerIntoVolume
  assert_not_nil          grid = Grid.new(7,4,1,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,1,2,1)
  assert_equal true,      grid.rightHandedFace(0)
  assert_equal true,      grid.rightHandedBoundary
  assert       0<         grid.minVolume, "negative volumes"
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal [0,1,2,3], grid.cell(0)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 7,         grid.nnode
  assert_equal [0,0,0.1], grid.nodeXYZ(4)
  assert_equal [1,0,0.1], grid.nodeXYZ(5)
  assert_equal [0,1,0.1], grid.nodeXYZ(6)
  assert_equal layer,     layer.wiggle(0.1)
  assert_equal [0,0,0.2], grid.nodeXYZ(4)
  assert_equal [1,0,0.2], grid.nodeXYZ(5)
  assert_equal [0,1,0.2], grid.nodeXYZ(6)
  assert_equal [4,5,6,3], grid.cell(0)
  assert_equal 4,         grid.ncell
  assert_equal true,      grid.rightHandedFace(0)
  assert_equal true,      grid.rightHandedBoundary
  assert       0<         grid.minVolume, "negative volumes"
 end

 def testAdvanceLayerIntoVolumeWithTwoFrontFacesOfVolumeCell
  assert_not_nil          grid = Grid.new(10,10,10,10)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(1,1,0.1)
  grid.addCell(0,1,2,3)
  grid.addFace(0,1,2,1)
  grid.addFace(1,3,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal [0,1,2,3], grid.cell(0)
  assert_equal false,     layer.cellInLayer(0)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [4,5,6,7], grid.cell(0)
  assert_equal false,     layer.cellInLayer(0)
  assert_equal true,      layer.cellInLayer(1)
  assert_equal true,      layer.cellInLayer(2)
  assert_equal true,      layer.cellInLayer(3)
 end

 def testAdvanceLayerIntoVolumeWithTwoFrontEdgeOfConstrainingFace
  assert_not_nil          grid = Grid.new(10,10,10,10)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(0,1,2,9)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal true,      layer.faceInLayer(0)
  assert_equal false,     layer.faceInLayer(1)
  assert_equal [0,1,2,1], grid.face(0)
  assert_equal [0,1,2,9], grid.face(1)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [0,1,2,1], grid.face(0)
  assert_equal [0,1,2,9], grid.face(1)
 end

 def testAdvanceLayerOnSymPlane
  assert_not_nil          grid = Grid.new(17,14,14,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,3,1,1)
  grid.addFace(0,1,2,2)
  assert_equal true,      grid.rightHandedBoundary
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal 2,         layer.constrained(0)
  assert_equal 0,         layer.constrained(1)
  assert_equal 2,         layer.constrained(2)
  assert_equal 0,         layer.constrainedSide(0,0)
  assert_equal 0,         layer.constrainedSide(0,1)
  assert_equal 2,         layer.constrainedSide(0,2)
  assert_equal [0,1,2,3], grid.cell(0)
  assert_equal [0,3,1,1], grid.face(0)
  assert_equal [0,1,2,2], grid.face(1)
  assert_equal true,      layer.faceInLayer(0)
  assert_equal false,     layer.faceInLayer(1)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 7,         grid.nnode
  assert_equal [4,6,2,5], grid.cell(0)
  assert_equal 4,         grid.ncell
  assert_equal [0,3,1,1], grid.face(0)
  assert_equal [4,6,2,2], grid.face(1)
  assert_equal [0,1,6,2], grid.face(2)
  assert_equal [0,6,4,2], grid.face(3)
  assert_equal true,      layer.faceInLayer(0)
  assert_equal false,     layer.faceInLayer(1)
  assert_equal true,      layer.faceInLayer(2)
  assert_equal true,      layer.faceInLayer(3)
  assert_equal 4,         grid.nface
  assert_equal true,      grid.rightHandedBoundary
 end

 def testAdvanceLayerOnSymPlane0
  assert_not_nil          grid = Grid.new(17,14,14,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,1,2,1)
  grid.addFace(0,3,1,2)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal 2,         layer.constrainedSide(0,0)
  assert_equal 0,         layer.constrainedSide(0,1)
  assert_equal 0,         layer.constrainedSide(0,2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [1,0,5,2], grid.face(2)
  assert_equal [5,0,4,2], grid.face(3)
 end

 def testAdvanceLayerOnSymPlane1
  assert_not_nil          grid = Grid.new(17,14,14,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,1,2,1)
  grid.addFace(1,3,2,2)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal 0,         layer.constrainedSide(0,0)
  assert_equal 2,         layer.constrainedSide(0,1)
  assert_equal 0,         layer.constrainedSide(0,2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [2,1,6,2], grid.face(2)
  assert_equal [6,1,5,2], grid.face(3)
 end

 def testAdvanceLayerOnSymPlane2
  assert_not_nil          grid = Grid.new(17,14,14,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,1,2,1)
  grid.addFace(2,3,0,2)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal 0,         layer.constrainedSide(0,0)
  assert_equal 0,         layer.constrainedSide(0,1)
  assert_equal 2,         layer.constrainedSide(0,2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [0,2,6,2], grid.face(2)
  assert_equal [0,6,4,2], grid.face(3)
 end

 def testAdvanceLayerOnSymPlane0TwoFace
  assert_not_nil          grid = Grid.new(17,14,14,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(1,1,0)
  assert_equal 4,         grid.addNode(0,0,1)
  grid.addFace(1,3,2,1)
  grid.addFace(0,1,2,1)
  grid.addFace(0,4,1,2)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal 0,         layer.constrainedSide(0,0)
  assert_equal 0,         layer.constrainedSide(0,1)
  assert_equal 0,         layer.constrainedSide(0,2)
  assert_equal 2,         layer.constrainedSide(1,0)
  assert_equal 0,         layer.constrainedSide(1,1)
  assert_equal 0,         layer.constrainedSide(1,2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [1,0,5,2], grid.face(3)
  assert_equal [5,0,8,2], grid.face(4)
 end

 def testAdvanceLayerTwiceOnSymPlane
  assert_not_nil          grid = Grid.new(15,17,16,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,3,1,1)
  grid.addFace(0,1,2,2)
  assert_equal true,      grid.rightHandedBoundary
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal [0,1,2,3], grid.cell(0)
  assert_equal [0,3,1,1], grid.face(0)
  assert_equal [0,1,2,2], grid.face(1)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 7,         grid.nnode
  assert_equal [4,6,2,5], grid.cell(0)
  assert_equal [0,5,6,4], grid.cell(1)
  assert_equal [1,0,5,6], grid.cell(2)
  assert_equal [1,0,3,5], grid.cell(3)
  assert_equal 4,         grid.ncell
  assert_equal [0,3,1,1], grid.face(0)
  assert_equal [4,6,2,2], grid.face(1)
  assert_equal [0,1,6,2], grid.face(2)
  assert_equal [0,6,4,2], grid.face(3)
  assert_equal 4,         grid.nface

  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal [7,9,2,8], grid.cell(0)
  assert_equal [0,5,6,4], grid.cell(1)
  assert_equal [1,0,5,6], grid.cell(2)
  assert_equal [1,0,3,5], grid.cell(3)
  assert_equal [4,8,9,7], grid.cell(4)
  assert_equal [4,5,9,8], grid.cell(5)
  assert_equal [6,4,5,9], grid.cell(6)
 
  assert_equal [0,3,1,1], grid.face(0)
  assert_equal [7,9,2,2], grid.face(1)
  assert_equal [0,1,6,2], grid.face(2)
  assert_equal [0,6,4,2], grid.face(3)
  assert_equal [4,6,9,2], grid.face(4)
  assert_equal [4,9,7,2], grid.face(5)

  assert_equal true,      grid.rightHandedBoundary
 end

 def testAdvanceLayerOnEdge
  assert_not_nil          grid = Grid.new(17,14,20,14)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  grid.addFace(0,2,3,1)
  grid.addFace(0,1,2,2)
  grid.addFace(0,3,1,3)
  assert_equal grid,      grid.setNGeomNode(2)
  assert_equal grid,      grid.setNGeomEdge(1)
  assert_equal grid,      grid.addGeomEdge(1,0,2)
  grid.addEdge(0,1,1,0.0,1.0)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal layer,     layer.constrainNormal(-1)
  assert_equal layer,     layer.constrainNormal(2)
  assert_equal layer,     layer.constrainNormal(3)
  assert_equal(-1,        layer.constrained(0))
  assert_equal 2,         layer.constrained(1)
  assert_equal 3,         layer.constrained(2)
  assert_equal [0,2,3,1], grid.face(0)
  assert_equal [0,1,2,2], grid.face(1)
  assert_equal [0,3,1,3], grid.face(2)
  assert_equal [0,1,1],   grid.edge(0)
  assert_equal false,     layer.edgeInLayer(0)
  assert_equal 0,         layer.nEdgeInLayer(-1)
  assert_equal 0,         layer.nEdgeInLayer(1)
  assert_equal 0,         layer.nEdgeInLayer(2)
  assert_equal 0,         layer.edgeEndPoint(1,0)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 7,         grid.nnode
  assert_equal 4,         grid.ncell
  assert_equal 7,         grid.nface
  assert_equal 2,         grid.nedge
  assert_equal [4,1,1],   grid.edge(0)
  assert_equal [0,4,1],   grid.edge(1)
  assert_equal false,     layer.edgeInLayer(0)
  assert_equal true,      layer.edgeInLayer(1)
  assert_equal 1,         layer.nEdgeInLayer(1)
  assert_equal 4,         layer.edgeEndPoint(1,0)
  assert_equal 0.0,       grid.nodeT(0,1)
  assert_equal 1.0,       grid.nodeT(1,1)
  assert_equal 0.1,       grid.nodeT(4,1)
  assert_equal [0,2,3,1], grid.face(0)
  assert_equal [4,1,5,2], grid.face(1)
  assert_equal [4,6,1,3], grid.face(2)
 end

 def testNormalTermination
  assert_not_nil          grid = Grid.new(7,4,1,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,1,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_equal false,     layer.normalTerminated(0)
  assert_equal false,     layer.normalTerminated(1)
  assert_equal false,     layer.normalTerminated(2)
  assert_equal true,      layer.anyActiveNormals
  assert_equal layer,     layer.terminateNormal(0)
  assert_equal layer,     layer.terminateNormal(1)
  assert_equal true,      layer.normalTerminated(0)
  assert_equal true,      layer.normalTerminated(1)
  assert_equal false,     layer.normalTerminated(2)
  assert_equal true,      layer.anyActiveNormals
  assert_equal layer,     layer.terminateNormal(2)
  assert_equal false,     layer.anyActiveNormals
 end
 
 def volumeGrid
  grid = Grid.new(7,4,1,0)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  grid.addFace(0,1,2,1)
  grid
 end

 def volumeLayer(grid)
  layer = Layer.new(grid)
  layer.populateAdvancingFront([1])
 end

 def testTerminateAll
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(0)
  assert_equal layer,     layer.terminateNormal(1)
  assert_equal layer,     layer.terminateNormal(2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 4,         grid.nnode
  assert_equal 1,         grid.ncell
 end

 def testTerminate01
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(0)
  assert_equal layer,     layer.terminateNormal(1)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 5,         grid.nnode
  assert_equal 2,         grid.ncell
  assert_equal [2,0,1,4], grid.cell(1)
 end

 def testTerminate12
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(1)
  assert_equal layer,     layer.terminateNormal(2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 5,         grid.nnode
  assert_equal 2,         grid.ncell
  assert_equal [0,1,2,4], grid.cell(1)
 end

 def testTerminate20
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(2)
  assert_equal layer,     layer.terminateNormal(0)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 5,         grid.nnode
  assert_equal 2,         grid.ncell
  assert_equal [0,1,2,4], grid.cell(1)
 end

 def testTerminate0
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(0)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 6,         grid.nnode
  assert_equal 3,         grid.ncell
  assert_equal [0,1,5,4], grid.cell(1)
  assert_equal [2,0,1,5], grid.cell(2)
 end

 def testTerminate1
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(1)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 6,         grid.nnode
  assert_equal 3,         grid.ncell
  assert_equal [0,1,5,4], grid.cell(1)
  assert_equal [2,0,1,5], grid.cell(2)
 end

 def testTerminate2
  assert_not_nil grid  = volumeGrid
  assert_not_nil layer = volumeLayer(grid)
  assert_equal layer,     layer.terminateNormal(2)
  assert_equal layer,     layer.advanceConstantHeight(0.1)
  assert_equal 6,         grid.nnode
  assert_equal 3,         grid.ncell
  assert_equal [0,5,2,4], grid.cell(1)
  assert_equal [0,1,2,5], grid.cell(2)
 end

 def testAdvanceLayerIntoVolumeWithVariableHeight
  assert_not_nil          grid = Grid.new(7,4,1,0)
  assert_equal 0,         grid.addNode(0,0,0)
  assert_equal 1,         grid.addNode(1,0,0)
  assert_equal 2,         grid.addNode(0,1,0)
  assert_equal 3,         grid.addNode(0,0,1)
  grid.addCell(0,1,2,3)
  assert_equal 1,         grid.ncell
  grid.addFace(0,1,2,1)
  assert_not_nil          layer = Layer.new(grid)
  assert_equal layer,     layer.populateAdvancingFront([1])
  assert_nil              layer.setNormalHeight(-1,0.0)
  assert_nil              layer.setNormalHeight(3,0.0)
  assert_equal layer,     layer.setNormalHeight(0,0.0)
  assert_equal layer,     layer.setNormalHeight(1,0.1)
  assert_equal layer,     layer.setNormalHeight(2,0.2)
  assert_equal true,      layer.normalTerminated(0)
  assert_equal false,     layer.normalTerminated(1)
  assert_equal false,     layer.normalTerminated(2)
  assert_equal layer,     layer.advance
  assert_equal 6,         grid.nnode
  assert_equal [0,0,0],   grid.nodeXYZ(0)
  assert_equal [1,0,0.1], grid.nodeXYZ(4)
  assert_equal [0,1,0.2], grid.nodeXYZ(5)
 end

 def testInitialNormalHeight
  grid = Grid.new(10,10,10,10)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_nil        layer.getNormalHeight(-1)
  assert_nil        layer.getNormalHeight(4)
  assert_equal 1.0, layer.getNormalHeight(0)
 end

 def testSetAllNormalHeights
  grid = Grid.new(10,10,10,10)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1])
  h = 0.546
  assert_equal layer, layer.setHeightOfAllNormals(h)
  assert_equal h, layer.getNormalHeight(0)
  assert_equal h, layer.getNormalHeight(1)
  assert_equal h, layer.getNormalHeight(2)
 end

 def testSetConstantNormalHeight
  grid = Grid.new(10,10,10,10)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal layer, layer.assignPolynomialNormalHeight(0.1,0.0,0.0,
                                                        [0.5,0,0],[1,0,0])
  assert_equal 1.0, layer.getNormalHeight(0)
  assert_equal 0.1, layer.getNormalHeight(1)
  assert_equal 1.0, layer.getNormalHeight(2)
 end


 def testSettingLinearNormalHeightDistribution
  grid = Grid.new(10,10,10,10)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal layer, layer.assignPolynomialNormalHeight(0.1,1.0,1.0,
                                                         [0.0,0,0],[0,1,0])
  assert_equal 0.1, layer.getNormalHeight(0)
  assert_equal 0.1, layer.getNormalHeight(1)
  assert_equal 1.1, layer.getNormalHeight(2)
 end

 def testSetSquareRootNormalHeightDistribution
  grid = Grid.new(10,10,10,10)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,4,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal layer, layer.assignPolynomialNormalHeight(0.2,1.0,0.5,
                                                         [0.0,0,0],[0,1,0])
  assert_equal 0.2, layer.getNormalHeight(0)
  assert_equal 0.2, layer.getNormalHeight(1)
  assert_equal 2.2, layer.getNormalHeight(2)
 end

 def testMixedElementModeToggleSwitch
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal true,          layer.tetrahedraOnly
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal false,         layer.tetrahedraOnly
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal true,          layer.tetrahedraOnly
 end

 def testAllowedToTerminateMixedElementNormals
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  layer = Layer.new(grid).populateAdvancingFront([1]).toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(0)
  assert_equal layer,         layer.terminateNormal(1)
  assert_equal layer,         layer.terminateNormal(2)
  assert_equal true,          layer.normalTerminated(0)
  assert_equal true,          layer.normalTerminated(1)
  assert_equal true,          layer.normalTerminated(2)
 end

 
 def testAdvanceLayerIntoVolumeWithaPrism
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal 3,             layer.nnormal
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 6,             grid.nnode
  assert_equal [0,0,0.1],     grid.nodeXYZ(3)
  assert_equal [1,0,0.1],     grid.nodeXYZ(4)
  assert_equal [0,1,0.1],     grid.nodeXYZ(5)
  assert_equal [0,1,2,3,4,5], grid.prism(0)
 end

 def testAdvanceLayerOnSymPlaneWithaQuadFace
  assert_not_nil               grid = Grid.new(17,14,14,0)
  assert_equal 0,              grid.addNode(0,0,0)
  assert_equal 1,              grid.addNode(1,0,0)
  assert_equal 2,              grid.addNode(0,1,0)
  assert_equal 3,              grid.addNode(0,0,1)
  grid.addFace(0,1,2,4000)
  grid.addFace(0,3,1,6000)
  assert_not_nil               layer = Layer.new(grid)
  assert_equal layer,          layer.populateAdvancingFront([4000])
  assert_equal layer,          layer.constrainNormal(6000)
  assert_equal layer,          layer.toggleMixedElementMode
  assert_equal 0,              grid.nquad
  assert_equal layer,          layer.advanceConstantHeight(0.1)
  assert_equal [1,0,4,5,6000], grid.quad(0)
  assert_equal 1,              grid.nquad
 end

 def testAdvanceLayerOnSymPlaneWithaTreminatedQuadFace0
  assert_not_nil               grid = Grid.new(17,14,14,0)
  assert_equal 0,              grid.addNode(0,0,0)
  assert_equal 1,              grid.addNode(1,0,0)
  assert_equal 2,              grid.addNode(0,1,0)
  assert_equal 3,              grid.addNode(0,0,1)
  grid.addFace(0,1,2,4000)
  grid.addFace(0,3,1,6000)
  assert_not_nil               layer = Layer.new(grid)
  assert_equal layer,          layer.populateAdvancingFront([4000])
  assert_equal layer,          layer.constrainNormal(6000)
  assert_equal layer,          layer.toggleMixedElementMode
  assert_equal 2,              grid.nface
  assert_equal layer,          layer.terminateNormal(0)
  assert_equal layer,          layer.advanceConstantHeight(0.1)
  assert_equal [1,0,4,6000],   grid.face(2)
  assert_equal 3,              grid.nface
 end

 def testAdvanceLayerOnSymPlaneWithaTreminatedQuadFace1
  assert_not_nil               grid = Grid.new(17,14,14,0)
  assert_equal 0,              grid.addNode(0,0,0)
  assert_equal 1,              grid.addNode(1,0,0)
  assert_equal 2,              grid.addNode(0,1,0)
  assert_equal 3,              grid.addNode(0,0,1)
  grid.addFace(0,1,2,4000)
  grid.addFace(0,3,1,6000)
  assert_not_nil               layer = Layer.new(grid)
  assert_equal layer,          layer.populateAdvancingFront([4000])
  assert_equal layer,          layer.constrainNormal(6000)
  assert_equal layer,          layer.toggleMixedElementMode
  assert_equal 2,              grid.nface
  assert_equal layer,          layer.terminateNormal(1)
  assert_equal layer,          layer.advanceConstantHeight(0.1)
  assert_equal [1,0,4,6000],   grid.face(2)
  assert_equal 3,              grid.nface
 end

 def testAdvanceLayerIntoVolumeWithaTerminatedPrismOrTet0
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(1)
  assert_equal layer,         layer.terminateNormal(2)
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 4,             grid.nnode
  assert_equal [0,0,0.1],     grid.nodeXYZ(3)
  assert_equal 1,             grid.ncell
  assert_equal [0,1,2,3],     grid.cell(0)
 end

 def testAdvanceLayerIntoVolumeWithaTerminatedPrismOrTet1
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(0)
  assert_equal layer,         layer.terminateNormal(2)
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 4,             grid.nnode
  assert_equal [1,0,0.1],     grid.nodeXYZ(3)
  assert_equal 1,             grid.ncell
  assert_equal [0,1,2,3],     grid.cell(0)
 end

 def testAdvanceLayerIntoVolumeWithaTerminatedPrismOrTet2
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(0)
  assert_equal layer,         layer.terminateNormal(1)
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 4,             grid.nnode
  assert_equal [0,1,0.1],     grid.nodeXYZ(3)
  assert_equal 1,             grid.ncell
  assert_equal [2,0,1,3],     grid.cell(0)
 end

 def testAdvanceLayerIntoVolumeWithaTerminatedPrismOrPyramid0
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(0)
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 5,             grid.nnode
  assert_equal 0,             grid.ncell
  assert_equal 0,             grid.nprism
  assert_equal 1,             grid.npyramid
  assert_equal [1,2,0,3,4],   grid.pyramid(0)
 end

 def testAdvanceLayerIntoVolumeWithaTerminatedPrismOrPyramid1
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(1)
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 5,             grid.nnode
  assert_equal 0,             grid.ncell
  assert_equal 0,             grid.nprism
  assert_equal 1,             grid.npyramid
  assert_equal [2,0,1,4,3],   grid.pyramid(0)
 end

 def testAdvanceLayerIntoVolumeWithaTerminatedPrismOrPyramid2
  assert_not_nil              grid = Grid.new(10,10,10,10)
  assert_equal 0,             grid.addNode(0,0,0)
  assert_equal 1,             grid.addNode(1,0,0)
  assert_equal 2,             grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  assert_not_nil              layer = Layer.new(grid)
  assert_equal layer,         layer.populateAdvancingFront([1])
  assert_equal layer,         layer.toggleMixedElementMode
  assert_equal layer,         layer.terminateNormal(2)
  assert_equal layer,         layer.advanceConstantHeight(0.1)
  assert_equal 5,             grid.nnode
  assert_equal 0,             grid.ncell
  assert_equal 0,             grid.nprism
  assert_equal 1,             grid.npyramid
  assert_equal [0,1,2,3,4],   grid.pyramid(0)
 end

 def flatTwoFaceGrid
  # y 2---3
  # ^ |0\1|
  # | 0---1 -> x
  grid = Grid.new(100,100,100,100)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addNode(1,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(2,1,3,1)
  grid
 end

 def testFindPreviousTriangleAroundNormal
  grid  = flatTwoFaceGrid
  layer = Layer.new(grid).populateAdvancingFront([1])

  assert_equal(-1, layer.previousTriangle(4,0) )
  assert_equal(-1, layer.previousTriangle(0,4) )

  assert_equal(-1, layer.previousTriangle(0,0) )
  assert_equal  1, layer.previousTriangle(1,0)
  assert_equal(-1, layer.previousTriangle(2,0) )

  assert_equal(-1, layer.previousTriangle(1,1) )
  assert_equal  0, layer.previousTriangle(2,1)
  assert_equal(-1, layer.previousTriangle(3,1) )
 end

 def testFindNextTriangleAroundNormal
  grid  = flatTwoFaceGrid
  layer = Layer.new(grid).populateAdvancingFront([1])

  assert_equal(-1, layer.nextTriangle(4,0) )
  assert_equal(-1, layer.nextTriangle(0,4) )

  assert_equal(-1, layer.nextTriangle(0,0) )
  assert_equal(-1, layer.nextTriangle(1,0) )
  assert_equal  1, layer.nextTriangle(2,0)

  assert_equal  0, layer.nextTriangle(1,1)
  assert_equal(-1, layer.nextTriangle(2,1) )
  assert_equal(-1, layer.nextTriangle(3,1) )
 end

 def testEdgeAngleGross
  grid  = flatTwoFaceGrid
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert            0> layer.edgeAngle(0,0)
  assert            0> layer.edgeAngle(3,0)
  assert            0> layer.edgeAngle(0,3)

  grid.setNodeXYZ(3,[0,0,0.012337])
  assert_in_delta   1, layer.edgeAngle(0,1), 1.0e-3
  grid.setNodeXYZ(3,[0,0,0.707])
  assert_in_delta  45, layer.edgeAngle(0,1), 0.005
  grid.setNodeXYZ(3,[0.5,0.5,1])
  assert_in_delta  90, layer.edgeAngle(0,1), 1.0e-10
  assert_in_delta  90, layer.edgeAngle(1,0), 1.0e-10
  grid.setNodeXYZ(3,[1,1,0.707])
  assert_in_delta 135, layer.edgeAngle(0,1), 0.005
  grid.setNodeXYZ(3,[1,1,0])
  assert_in_delta 180, layer.edgeAngle(0,1), 1.0e-10
  grid.setNodeXYZ(3,[1,1,-0.707])
  assert_in_delta 225, layer.edgeAngle(0,1), 0.005
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  assert_in_delta 270, layer.edgeAngle(0,1), 1.0e-10
  assert_in_delta 270, layer.edgeAngle(1,0), 1.0e-10
  grid.setNodeXYZ(3,[0,0,-0.707])
  assert_in_delta 315, layer.edgeAngle(0,1), 0.005
  grid.setNodeXYZ(3,[0,0,-0.012337])
  assert_in_delta 359, layer.edgeAngle(0,1), 1.0e-3
 end

 def testEdgeAngleFine
  grid  = flatTwoFaceGrid
  layer = Layer.new(grid).populateAdvancingFront([1])

  wiggle = 1.0e-10
  tol    = 1.0e-8
  
  grid.setNodeXYZ(3,[0,0,wiggle])
  assert_in_delta   0, layer.edgeAngle(0,1), tol
  grid.setNodeXYZ(3,[0,0,-wiggle])
  assert_in_delta 360, layer.edgeAngle(0,1), tol

  grid.setNodeXYZ(3,[0.5,0.5,1+1.0e-10])
  assert_in_delta  90, layer.edgeAngle(0,1), tol
  grid.setNodeXYZ(3,[0.5,0.5,1-1.0e-10])
  assert_in_delta  90, layer.edgeAngle(0,1), tol

  grid.setNodeXYZ(3,[1,1,wiggle])
  assert_in_delta 180, layer.edgeAngle(0,1), tol
  grid.setNodeXYZ(3,[1,1,-wiggle])
  assert_in_delta 180, layer.edgeAngle(0,1), tol

  grid.setNodeXYZ(3,[0.5,0.5,-1+1.0e-10])
  assert_in_delta 270, layer.edgeAngle(0,1), tol
  grid.setNodeXYZ(3,[0.5,0.5,-1-1.0e-10])
  assert_in_delta 270, layer.edgeAngle(0,1), tol
 end

 def testFindFirstFaceOffSymm
  grid = Grid.new(20,20,10,0)
  grid.addNode(0.5,0.35,0)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)

  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  assert_equal 0, layer.firstTriangleAfterGap(0)

  grid.addFace(0,3,1,10)
  layer = Layer.new(grid).populateAdvancingFront([10])
  assert_equal 2, layer.firstTriangleAfterGap(0)
 end

 def testInsertBlendForConvextFace
  grid  = flatTwoFaceGrid
  # y 2---3
  # ^ |0\1|
  # | 0---1 -> x
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal [0,1,2], layer.triangleNormals(0)
  assert_equal [2,1,3], layer.triangleNormals(1)
  assert_equal [0,1,2], layer.triangle(0)
  assert_equal [2,1,3], layer.triangle(1)
  assert_equal 0,     layer.nblend
  assert_equal 0,     layer.nRequiredBlends(0,-1.0)
  assert_equal 1,     layer.nRequiredBlends(1,-1.0)
  assert_equal 1,     layer.nRequiredBlends(2,-1.0)
  assert_equal 0,     layer.nRequiredBlends(3,-1.0)
  assert_equal 0,     layer.blendDegree(-1)
  assert_equal 0,     layer.blendDegree(0)
  assert_equal layer, layer.blend(-1.0)
  assert_equal 1,     layer.nblend
  # y 25--3 normals
  # ^ |0\1|
  # | 0--41 -> x
  assert_equal [0,1,2], layer.triangle(0)
  assert_equal [2,1,3], layer.triangle(1)
  assert_equal [4,1,2,5], layer.blendNormals(0)
  assert_equal [0,4,2], layer.triangleNormals(0)
  assert_equal [5,1,3], layer.triangleNormals(1)
  assert_equal 0,       layer.blendDegree(0)
  assert_equal 1,       layer.blendDegree(1)
  assert_equal 0,       layer.blendDegree(4)
 end

 def testAdvanceBlendForConvextFace
  grid  = flatTwoFaceGrid
  # y 2---3
  # ^ |0\1|
  # | 0---1 -> x
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(-1.0)
  assert_equal layer, layer.advanceConstantHeight(0.1)
  assert_equal 4, layer.ntriangle
  assert_equal 9, grid.ncell
  assert_equal 0, layer.nblend
  assert_equal 0, layer.blendDegree(1)
 end


 def fourFaceConvex
  #        2
  #      / y \
  #    /4  ^  0\  
  #  /     |     \
  # 4 -z<- 0 -> x 1
  #  \     |     /
  #    \3  |  1/
  #      \ | / 
  #        3

  grid = Grid.new(10,10,10,10)
  grid.addNode(0, 0, 0)
  grid.addNode(1, 0, 0)
  grid.addNode(0, 1, 0)
  grid.addNode(0,-1, 0)
  grid.addNode(0.01, 0,-1)
  grid.addFace(0,1,2,1)
  grid.addFace(0,3,1,1)
  grid.addFace(0,4,3,1)
  grid.addFace(0,2,4,1)
  grid
 end

 def testBlendTerminatedWithCommonNormal01
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(0.5,0.5,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(1,1,0)
  grid.addNode(0,1,0)
  grid.addNode(-0.5,0.5,top)
  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)
  grid.addFace(0,3,4,10)
  grid.addFace(0,4,5,10)
  grid.addFace(0,5,1,10)
  layer = Layer.new(grid).populateAdvancingFront([10])
  assert_equal 6, layer.nnormal
  layer.blend(275.0)
  assert_equal 7, layer.nnormal  
  assert_equal 1, layer.nblend 
  assert_equal [0,0,5,6], layer.blendNormals(0)
  assert_equal 5,  layer.ntriangle
  layer.advanceConstantHeight(0.1)
  assert_equal 6,  layer.ntriangle
  assert_equal 17, grid.ncell
 end

 def testBlendTerminatedWithCommonNormal10
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(-0.5,0.5,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(1,1,0)
  grid.addNode(0,1,0)
  grid.addNode(0.5,0.5,top)
  grid.addFace(5,1,2,10)
  grid.addFace(5,2,3,10)
  grid.addFace(5,3,4,10)
  grid.addFace(5,4,0,10)
  grid.addFace(5,0,1,10)
  layer = Layer.new(grid).populateAdvancingFront([10])
  layer.blend(275.0)
  assert_equal 7, layer.nnormal  
  assert_equal [0,0,5,6], layer.blendNormals(0)
  layer.advanceConstantHeight(0.1)
  assert_equal 6,  layer.ntriangle
  assert_equal 16, grid.ncell
 end

 def testInsertBlendForTwoConvextFaces
  grid = fourFaceConvex
  #        2
  #      / y \
  #    /4  ^  0\  
  #  /     |     \
  # 4 -z<- 0 -> x 1
  #  \     |     /
  #    \3  |  1/
  #      \ | / 
  #        3
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal [0,1,2], layer.triangleNormals(0)
  assert_equal [0,3,1], layer.triangleNormals(1)
  assert_equal [0,4,3], layer.triangleNormals(2)
  assert_equal [0,2,4], layer.triangleNormals(3)
  assert_equal 2,     layer.nRequiredBlends(0,-1.0)
  assert_equal layer, layer.blend(-1.0)
  assert_equal 2,     layer.nblend
  assert_equal [5,1,6], layer.triangleNormals(0)
  assert_equal [5,3,1], layer.triangleNormals(1)
  assert_equal [0,4,7], layer.triangleNormals(2)
  assert_equal [0,2,4], layer.triangleNormals(3)
  #       2 6
  #      / y \
  #    /4  ^  0\  
  #  /     |     \
  # 4 -z<-0 5-> x 1
  #  \     |     /
  #    \3  |  1/
  #      \ | / 
  #       7 3
  assert_equal [0,5,2,6], layer.blendNormals(1)
  assert_equal [5,0,3,7], layer.blendNormals(0)
 end

 def testAdvanceBlendForConvextFace_ConvertBlendToTriangle
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(-1.0)
  assert_equal 1,       layer.nblend
  layer.advance
  assert_equal 0,       layer.nblend
  assert_equal 4,       layer.ntriangle
  assert_equal [4,1,5], layer.triangleNormals(2)
  assert_equal [5,2,4], layer.triangleNormals(3)
 end

 def testCheckAfterBlendForConvextWithConstrainingFace
  grid = Grid.new(100,100,100,100)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,0,1)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(0,2,3,1)
  grid.addFace(1,0,3,11)
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.constrainNormal(11)
  grid.removeFace(2)
  layer.blend(200.0)

  assert_equal [0,1,5], layer.triangleNormals(0)
  assert_equal [4,2,3], layer.triangleNormals(1)

  # Y    3
  # |    |\
  # +--Z 4-2
  # +--Z 0-5
  # |    |/
  # X    1

  assert_equal [0.0,-1.0,0.0], layer.normalTriangleDirection(0,0)
  assert_equal [0.0,-1.0,0.0], layer.normalTriangleDirection(1,0)
  assert_equal [0.0,-1.0,0.0], layer.normalTriangleDirection(5,0)

  assert_equal [-1.0,0.0,0.0], layer.normalTriangleDirection(4,0)
  assert_equal [-1.0,0.0,0.0], layer.normalTriangleDirection(2,0)
  assert_equal [-1.0,0.0,0.0], layer.normalTriangleDirection(3,0)

 end

 def testAdvanceBlendForConvextFace_AddConstrainingFace
  # Y    3
  # |    |\
  # +--Z 0-2
  # |    |/
  # X    1
  grid = Grid.new(100,100,100,100)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,0,1)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(0,2,3,1)
  grid.addFace(1,0,3,11)
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.constrainNormal(11)
  grid.removeFace(2)
  layer.blend(200.0)
  assert_equal 2, grid.nface
  layer.advanceConstantHeight(0.1)
  assert_equal 9, grid.ncell
  assert_equal 7, grid.nface
 end

 def testCheckSubBlendNormalDirectionForConvextWithConstrainingFace
  grid = Grid.new(100,100,100,100)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,0,1)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(0,2,3,1)
  grid.addFace(1,0,3,11)
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.constrainNormal(11)
  grid.removeFace(2)
  layer.blend(200.0)

  layer.subBlend(44.0)
  assert_equal [4,6,2,7], layer.subBlendNormals(0,0)
  assert_equal [6,0,7,5], layer.subBlendNormals(0,1)

  # Y    3
  # |    |\
  # +--Z 4-2
  #      6-7
  # +--Z 0-5
  # |    |/
  # X    1

  tol = 1.0e-3
  norm = layer.normalDirection(6)
  assert_in_delta( -0.707, norm[0], tol)
  assert_in_delta( -0.707, norm[1], tol)
  assert_in_delta(  0.000, norm[2], tol)
  norm = layer.normalDirection(7)
  assert_in_delta( -0.707, norm[0], tol)
  assert_in_delta( -0.707, norm[1], tol)
  assert_in_delta(  0.000, norm[2], tol)
 end

 # deactivated, it will seg fault when blend axle call with new blend normal
 def XtestForCollidingNormalsBetweenBlendAndConcaveNeighbor
  grid = Grid.new(100,100,100,100)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,0,1)
  grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(0,2,3,1)
  grid.addNode(0,1,1)
  grid.addNode(-1,1,0)
  grid.addFace(3,4,5,1)
  grid.addFace(1,0,3,11)
  grid.addFace(1,0,3,11)
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.constrainNormal(11)
  grid.removeFace(3)
  grid.removeFace(4)
  layer.blend(200.0)
  layer.preventBlendNormalDirectionFromPointingAtNeighbors(0.4)

  assert_equal [0,1,7], layer.triangleNormals(0)
  assert_equal [6,2,3], layer.triangleNormals(1)
  assert_equal [3,4,5], layer.triangleNormals(2)

  #      5
  #      |\
  # Y    3-4
  # |    |\
  # +--Z 6-2
  # +--Z 0-7
  # |    |/
  # X    1

  tol = 1.0e-5
  norm = layer.normalDirection(3)
  assert_in_delta( -0.707107, norm[0], tol)
  assert_in_delta( -0.707107, norm[1], tol)
  assert_in_delta(  0.000000, norm[2], tol)
  norm = layer.normalDirection(6)
  assert_in_delta( -0.707107, norm[0], tol, "6-X")
  assert_in_delta( -0.707107, norm[1], tol, "6-Y")
  assert_in_delta(  0.000000, norm[2], tol, "6-Z")
  norm = layer.normalDirection(2)
  assert_in_delta( -0.707107, norm[0], tol, "3-X")
  assert_in_delta( -0.707107, norm[1], tol, "3-Y")
  assert_in_delta(  0.000000, norm[2], tol, "3-Z")

 end
#layer.writeTecplotFrontGeometry
#grid.writeTecplotSurfaceZone

 def XtestForCollidingNormalsBetweenBlendAndConcaveNeighborSharp
  grid = Grid.new(100,100,100,100)
  assert_equal 0, grid.addNode(0,0,0)
  assert_equal 1, grid.addNode(0.1,1,0)
  assert_equal 2, grid.addNode(0,0,1)
  assert_equal 3, grid.addNode(0,1,0)
  grid.addFace(0,1,2,1)
  grid.addFace(0,2,3,1)
  assert_equal 4, grid.addNode(0,1,1)
  assert_equal 5, grid.addNode(-0.707,0.293,0)
  grid.addFace(3,4,5,1)
  grid.addFace(1,0,3,11)
  grid.addFace(1,0,3,11)
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.constrainNormal(11)
  grid.removeFace(3)
  grid.removeFace(4)
  layer.blend(200.0)
  layer.preventBlendNormalDirectionFromPointingAtNeighbors(0.4)

  layer.writeTecplotFrontGeometry

  assert_equal [0,1,7], layer.triangleNormals(0)
  assert_equal [6,2,3], layer.triangleNormals(1)
  assert_equal [3,4,5], layer.triangleNormals(2)

  #      5
  #      |\
  # Y    3-4
  # |    |\
  # +--Z 6-2
  # +--Z 0-7
  # |    |/
  # X    1

 end
#layer.writeTecplotFrontGeometry
#grid.writeTecplotSurfaceZone

 def testBlendTriplePoint
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(0.5,0.35,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)

  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)
  grid.addFace(0,3,1,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  assert_equal 3,     layer.nRequiredBlends(0,270.0)
  layer.blend(270.0)
  assert_equal 10, layer.nnormal  
  assert_equal 3, layer.nblend 
  assert_equal [5,1,8], layer.triangleNormals(0)
  assert_equal [6,2,9], layer.triangleNormals(1)
  assert_equal [0,3,7], layer.triangleNormals(2)
  assert_equal [5,0,1,7], layer.blendNormals(0)
  assert_equal [6,5,2,8], layer.blendNormals(1)
  assert_equal [0,6,3,9], layer.blendNormals(2)

#     /3   9\
#    / |   | \
#    / 0   6 \
#   / /     \ \
#   //   5   \\
#   7   / \   2
#      1---8

  assert_equal 3,         layer.blendDegree(0)
  assert_equal [0, 1, 2], layer.orderedVertexBlends(0)
  assert_equal [0, 5, 6], layer.orderedVertexNormals(0)
  assert_equal layer, layer.advanceConstantHeight(0.1)
  #layer.writeTecplotFrontGeometry
  assert_equal 12, layer.ntriangle  
  assert_equal 21, grid.ncell  
 end

 def testBlendTriplePointDisorder
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(54.9355,37.5639,9.71613)
  grid.addNode(54.985,37.5706,9.71519)
  grid.addNode(54.9348,37.5737,9.71417)
  grid.addNode(54.9606,37.5627,9.66953)

  grid.addFace(0,1,2,10)
  grid.addFace(0,3,1,10)
  grid.addFace(2,1,3,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  layer.blend(240.0)
  layer.subBlend(30.0)
  #assert_equal [0, 5, 6], layer.orderedVertexNormals(1)
  layer.advanceConstantHeight(0.01)
  #layer.writeTecplotFrontGeometry
 end

 def testSubBlendCount
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_in_delta 270.0, layer.edgeAngle(0,1), 1.0e-8
  layer.blend(180.0)
  assert_equal 1, layer.nblend
  assert_equal EMPTY, layer.nSubBlend(-1)
  assert_equal 1, layer.nSubBlend(0)
  assert_equal EMPTY, layer.nSubBlend(1)
  assert_equal layer, layer.subBlend(44.0)
  assert_equal 2, layer.nSubBlend(0)
 end

 def testSubBlendNormals2
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(180.0)
  layer.subBlend(44.0)
  assert_equal [4,1,2,5], layer.blendNormals(0)
  assert_equal 2, layer.nSubBlend(0)
  assert_equal [4,6,2,7], layer.subBlendNormals(0,0)
  assert_equal [6,1,7,5], layer.subBlendNormals(0,1)
 end

 def testSubBlendNormals3
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(180.0)
  layer.subBlend(29.0)
  assert_equal [4,1,2,5], layer.blendNormals(0)
  assert_equal 3, layer.nSubBlend(0)
  assert_equal [4,6,2,8], layer.subBlendNormals(0,0)
  assert_equal [6,7,8,9], layer.subBlendNormals(0,1)
  assert_equal [7,1,9,5], layer.subBlendNormals(0,2)
 end

 def testsubBlendForTwoConvextFaces2
  grid = fourFaceConvex
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal layer, layer.blend(-1.0)
  assert_equal 2,     layer.nblend
  #       2 6
  #      / y \
  #    /4  ^  0\  
  #  /     |     \
  # 4 -z<-0 5-> x 1
  #  \     |     /
  #    \3  |  1/
  #      \ | / 
  #       7 3
  layer.subBlend(44.0)
  assert_equal 2, layer.nSubBlend(0)
  assert_equal 2, layer.nSubBlend(1)
  assert_equal [5,0,3,7], layer.blendNormals(0)
  assert_equal [0,5,2,6], layer.blendNormals(1)
  assert_equal [5,8,3,10], layer.subBlendNormals(0,0)
  assert_equal [8,0,10,7], layer.subBlendNormals(0,1)
  assert_equal [0,8,2,9], layer.subBlendNormals(1,0)
  assert_equal [8,5,9,6], layer.subBlendNormals(1,1)
 end

 def testsubBlendForTwoConvextFaces3
  grid = fourFaceConvex
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal layer, layer.blend(-1.0)
  assert_equal 2,     layer.nblend
  layer.subBlend(29.0)
  assert_equal 3, layer.nSubBlend(0)
  assert_equal 3, layer.nSubBlend(1)

  assert_equal [ 5, 0, 3, 7], layer.blendNormals(0)

  assert_equal [ 5, 9, 3,12], layer.subBlendNormals(0,0)
  assert_equal [ 9, 8,12,13], layer.subBlendNormals(0,1)
  assert_equal [ 8, 0,13, 7], layer.subBlendNormals(0,2)

  assert_equal [ 0, 5, 2, 6], layer.blendNormals(1)

  assert_equal [ 0, 8, 2,10], layer.subBlendNormals(1,0)
  assert_equal [ 8, 9,10,11], layer.subBlendNormals(1,1)
  assert_equal [ 9, 5,11, 6], layer.subBlendNormals(1,2)
 end

 def testSubBlendTriplePoint2
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(0.5,0.35,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)

  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)
  grid.addFace(0,3,1,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  layer.blend(270.0)
  layer.subBlend(44.0)
  assert_equal 16, layer.nnormal  
  assert_equal 2, layer.nSubBlend(0)
  assert_equal 2, layer.nSubBlend(1)
  assert_equal 2, layer.nSubBlend(2)
  assert_equal [ 5,10, 1,13], layer.subBlendNormals(0,0)
  assert_equal [10, 0,13, 7], layer.subBlendNormals(0,1)
  assert_equal [ 6,11, 2,14], layer.subBlendNormals(1,0)
  assert_equal [11, 5,14, 8], layer.subBlendNormals(1,1)
  assert_equal [ 0,12, 3,15], layer.subBlendNormals(2,0)
  assert_equal [12, 6,15, 9], layer.subBlendNormals(2,1)

  assert_equal [0, 1, 2], layer.orderedVertexBlends(0)
  assert_equal [0,10, 5, 11, 6, 12], layer.orderedVertexNormals(0)
 end

 def testSubBlendTriplePoint3
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(0.5,0.35,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)

  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)
  grid.addFace(0,3,1,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  layer.blend(270.0)
  layer.subBlend(29.0)
  assert_equal 22, layer.nnormal  
  assert_equal 3, layer.nSubBlend(0)
  assert_equal 3, layer.nSubBlend(1)
  assert_equal 3, layer.nSubBlend(2)
  assert_equal [ 5,10, 1,16], layer.subBlendNormals(0,0)
  assert_equal [10,11,16,17], layer.subBlendNormals(0,1)
  assert_equal [11, 0,17, 7], layer.subBlendNormals(0,2)

  assert_equal [ 6,12, 2,18], layer.subBlendNormals(1,0)
  assert_equal [12,13,18,19], layer.subBlendNormals(1,1)
  assert_equal [13, 5,19, 8], layer.subBlendNormals(1,2)

  assert_equal [ 0,14, 3,20], layer.subBlendNormals(2,0)
  assert_equal [14,15,20,21], layer.subBlendNormals(2,1)
  assert_equal [15, 6,21, 9], layer.subBlendNormals(2,2)

  assert_equal [0, 11, 10, 5, 13, 12, 6, 15, 14], layer.orderedVertexNormals(0)
 end

 def testsubBlendForTwoConvextFaces2Advance
  grid = fourFaceConvex
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(-1.0)
  layer.subBlend(44.0)
  layer.advanceConstantHeight(0.1)
  assert_equal 12, layer.ntriangle
  assert_equal 24, grid.ncell
 end

 def testsubBlendForTwoConvextFaces3Advance
  grid = fourFaceConvex
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(-1.0)
  layer.subBlend(29.0)
  layer.advanceConstantHeight(0.1)
  assert_equal 16, layer.ntriangle
  assert_equal 30, grid.ncell
 end

 def testVaribleSubBlendForTwoConvextFaces12Advance
  grid = fourFaceConvex
  grid.setNodeXYZ(2,[0.3,1,-0.3])
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(-1.0)
  layer.subBlend(44.0)
  assert_equal [0, 5, 2, 6], layer.blendNormals(1)
  assert_equal [0, 8, 2, 6], layer.subBlendNormals(1,0)
  assert_equal [8, 5, 6, 6], layer.subBlendNormals(1,1)
  layer.advanceConstantHeight(0.1)
  assert_equal 11, layer.ntriangle
 end

 def testSubBlendTriplePoint2Advance
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(0.5,0.35,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)

  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)
  grid.addFace(0,3,1,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  layer.blend(270.0)
  layer.subBlend(44.0)
  assert_equal 16, layer.nnormal
  assert_equal layer, layer.advanceConstantHeight(0.1)
  assert_equal 21, layer.ntriangle
 end

 def testSubBlendTriplePoint3Advance
  grid = Grid.new(20,20,10,0)
  top = 0.8
  grid.addNode(0.5,0.35,top)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0.5,0.7,0)

  grid.addFace(0,1,2,10)
  grid.addFace(0,2,3,10)
  grid.addFace(0,3,1,10)

  layer = Layer.new(grid).populateAdvancingFront([10])
  layer.blend(270.0)
  layer.subBlend(29.0)
  assert_equal 22, layer.nnormal
  assert_equal layer, layer.advanceConstantHeight(0.1)
  assert_equal 30, layer.ntriangle
 end

 def testExtrudeBlend
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0.5,0.5,-1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  layer.blend(-1.0)
  assert_equal 1,     layer.nblend
  assert_equal 2,     layer.ntriangle
  assert_equal 6,     layer.nnormal
  assert_equal 4,     grid.nnode
  assert_equal [4,1,2,5], layer.blendNormals(0)
  assert_equal layer, layer.extrudeBlend(1,1,1)
  assert_equal 1,     layer.nblend
  assert_equal 6,     layer.ntriangle
  assert_equal 10,    layer.nnormal
  assert_equal 6,     grid.nnode
  assert_equal [2,1,1], grid.nodeXYZ(4)
  assert_equal [1,2,1], grid.nodeXYZ(5)

  # y 25--3 normals
  # ^ |0\1|
  # | 0--41 -> x
  assert_equal [0,1,2], layer.triangle(0)
  assert_equal [2,1,3], layer.triangle(1)
  assert_equal [0,4,2], layer.triangleNormals(0)
  assert_equal [5,1,3], layer.triangleNormals(1)

  assert_equal [6,7,8,9], layer.blendNormals(0)

  # y 
  # ^    5
  # | 2   \  wake nodes
  # |  \   4
  # |   1
  assert_equal [1,4,2], layer.triangle(2)
  assert_equal [2,4,5], layer.triangle(3)
  assert_equal [2,5,4], layer.triangle(4)
  assert_equal [2,4,1], layer.triangle(5)

  # y 
  # ^    8
  # | 2   \  wake normals
  # |  \   6
  # |   4
  assert_equal [4,6,2], layer.triangleNormals(2)
  assert_equal [2,6,8], layer.triangleNormals(3)


  # y 
  # ^    9
  # | 5   \  wake normals
  # |  \   7
  # |   1
  assert_equal [5,9,7], layer.triangleNormals(4)
  assert_equal [5,7,1], layer.triangleNormals(5)
 end

 def facingGrid(z)
  # y 2
  # ^ |0\
  # | 0---1 -> x
  grid = Grid.new(100,100,100,100)
  grid.addNode(0,0,0)
  grid.addNode(1,0,0)
  grid.addNode(0,1,0)
  grid.addNode(0,0,z)
  grid.addNode(1,0,z)
  grid.addNode(0,1,z)
  grid.addFace(0,1,2,1)
  grid.addFace(3,5,4,1)
  grid
 end

 def testCollideNormalFarApart
  grid  = facingGrid(5)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal 6, layer.nActiveNormal
  layer.terminateCollidingNormals
  assert_equal 6, layer.nActiveNormal
 end

 def testCollideNormalClose
  grid  = facingGrid(0.5)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal 6, layer.nActiveNormal
  layer.terminateCollidingNormals
  assert_equal 0, layer.nActiveNormal
 end

 def testDoNotCollideNormalAtSharpTrailingEdge
  wiggle = 1.0e-10
  tol    = 1.0e-8
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0,0,-wiggle])
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_in_delta 360, layer.edgeAngle(0,1), tol
  assert_equal 4, layer.nActiveNormal
  layer.terminateCollidingNormals
  assert_equal 4, layer.nActiveNormal
 end

 def testCollideTriangleFarApart
  grid  = facingGrid(5)
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_equal 6, layer.nActiveNormal
  layer.terminateCollidingTriangles(1.0)
  assert_equal 6, layer.nActiveNormal
 end

 def testCollideTriangleClose
  grid  = facingGrid(0.5)
  layer = Layer.new(grid).populateAdvancingFront([1]).setHeightOfAllNormals(0.1)
  assert_equal 6, layer.nActiveNormal
  layer.terminateCollidingTriangles(1.0)
  assert_equal 0, layer.nActiveNormal
 end

 def testDoNotCollideTriangleAtSharpTrailingEdge
  wiggle = 1.0e-10
  tol    = 1.0e-8
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0,0,-wiggle])
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_in_delta 360, layer.edgeAngle(0,1), tol
  assert_equal 4, layer.nActiveNormal
  layer.terminateCollidingTriangles(1.0)
  assert_equal 4, layer.nActiveNormal
 end

 def testDoNotCollideTriangleForFlatFace
  wiggle = 1.0e-10
  tol    = 1.0e-8
  grid  = flatTwoFaceGrid
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_in_delta 180, layer.edgeAngle(0,1), tol
  assert_equal 4, layer.nActiveNormal
  layer.terminateCollidingTriangles(1.0)
  assert_equal 4, layer.nActiveNormal
 end

 def testDoNotCollideTriangleFor90degCorner
  wiggle = 1.0e-10
  tol    = 1.0e-8
  grid  = flatTwoFaceGrid
  grid.setNodeXYZ(3,[0.5,0.5,1])
  layer = Layer.new(grid).populateAdvancingFront([1])
  assert_in_delta 90, layer.edgeAngle(0,1), tol
  assert_equal 4, layer.nActiveNormal
  layer.terminateCollidingTriangles(1.0)
  assert_equal 4, layer.nActiveNormal
 end

end
