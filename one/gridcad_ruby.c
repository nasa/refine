
#include "ruby.h"
#include "gridcad.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_projectNodeToEdge( VALUE self, VALUE node, VALUE edgeId )
{
  GET_GRID_FROM_SELF;
  return (gridProjectNodeToEdge( grid, NUM2INT(node), NUM2INT(edgeId) )==grid?self:Qnil);
}

VALUE grid_projectNodeToFace( VALUE self, VALUE node, VALUE faceId )
{
  GET_GRID_FROM_SELF;
  return (gridProjectNodeToFace( grid, NUM2INT(node), NUM2INT(faceId) )==grid?self:Qnil);
}

VALUE grid_evaluateEdgeAtT( VALUE self, VALUE node, VALUE t )
{
  GET_GRID_FROM_SELF;
  return (gridEvaluateEdgeAtT(grid, NUM2INT(node), NUM2DBL(t))==grid?self:Qnil);
}

VALUE grid_evaluateFaceAtUV( VALUE self, VALUE node, VALUE u, VALUE v )
{
  double uv[2];
  GET_GRID_FROM_SELF;
  uv[0] = NUM2DBL(u);
  uv[1] = NUM2DBL(v);
  return (gridEvaluateFaceAtUV(grid, NUM2INT(node), uv )==grid?self:Qnil);
}

VALUE grid_updateFaceParameters( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridUpdateFaceParameters( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE grid_projectNode( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridProjectNode( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE grid_smooth( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSmooth( grid, -1.0, -1.0 )==grid?self:Qnil);
}

VALUE grid_smoothNode( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSmoothNode( grid, NUM2INT(node), TRUE )==grid?self:Qnil);
}

VALUE grid_lineSearchT( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridLineSearchT( grid, NUM2INT(node), gridMinSurfaceSmoothCost(grid) )==grid?self:Qnil);
}

VALUE grid_linearProgramUV( VALUE self, VALUE node )
{
  GridBool callAgain;
  GET_GRID_FROM_SELF;
  return (gridLinearProgramUV( grid, NUM2INT(node), &callAgain )==grid?self:Qnil);
}

VALUE grid_smartLaplacian( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSmartLaplacian( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE grid_storeVolumeCostDerivatives( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (grid == gridStoreVolumeCostDerivatives( grid, 
						  NUM2INT(node) )?self:Qnil);
}

VALUE grid_storeFaceCostParameterDerivatives( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (grid == gridStoreFaceCostParameterDerivatives( grid, 
							NUM2INT(node) 
							)?self:Qnil);
}

VALUE grid_relaxNegativeCells( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridRelaxNegativeCells( grid, FALSE )==grid?self:Qnil);
}

VALUE grid_smoothNodeFaceAreaUV( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSmoothNodeFaceAreaUV( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE grid_untangleAreaUV( VALUE self, VALUE node )
{
  int recursive_depth = 0;
  GridBool allow_movement_near_ghost_nodes = FALSE;
  GET_GRID_FROM_SELF;
  return (gridUntangleAreaUV( grid, NUM2INT(node), 
			      recursive_depth, 
			      allow_movement_near_ghost_nodes )==grid?self:Qnil);
}

VALUE grid_untangleVolume( VALUE self, VALUE node )
{
  int recursive_depth = 0;
  GridBool allow_movement_near_ghost_nodes = FALSE;
  GET_GRID_FROM_SELF;
  return (gridUntangleVolume( grid, NUM2INT(node), 
			      recursive_depth, 
			      allow_movement_near_ghost_nodes )==grid?self:Qnil);
}

VALUE cGridCAD;

void Init_GridCAD() 
{
  cGridCAD = rb_define_module( "GridCAD" );
  rb_define_method( cGridCAD, "projectNodeToEdge", grid_projectNodeToEdge, 2 );
  rb_define_method( cGridCAD, "projectNodeToFace", grid_projectNodeToFace, 2 );

  rb_define_method( cGridCAD, "evaluateEdgeAtT", grid_evaluateEdgeAtT, 2 );
  rb_define_method( cGridCAD, "evaluateFaceAtUV", grid_evaluateFaceAtUV, 3 );
  rb_define_method( cGridCAD, "updateFaceParameters", 
		    grid_updateFaceParameters, 1 );

  rb_define_method( cGridCAD, "projectNode", grid_projectNode, 1 );

  rb_define_method( cGridCAD, "smooth", grid_smooth, 0 );
  rb_define_method( cGridCAD, "smoothNode", grid_smoothNode, 1 );

  rb_define_method( cGridCAD, "lineSearchT", grid_lineSearchT, 1 );
  rb_define_method( cGridCAD, "linearProgramUV", grid_linearProgramUV, 1 );

  rb_define_method( cGridCAD, "smartLaplacian", grid_smartLaplacian, 1 );

  rb_define_method( cGridCAD, "storeVolumeCostDerivatives",
		    grid_storeVolumeCostDerivatives, 1 );
  rb_define_method( cGridCAD, "storeFaceCostParameterDerivatives",
		    grid_storeFaceCostParameterDerivatives, 1 );

  rb_define_method( cGridCAD, "relaxNegativeCells", grid_relaxNegativeCells, 0 );
  rb_define_method( cGridCAD, "smoothNodeFaceAreaUV", grid_smoothNodeFaceAreaUV, 1 );
  rb_define_method( cGridCAD, "untangleAreaUV", grid_untangleAreaUV, 1 );
  rb_define_method( cGridCAD, "untangleVolume", grid_untangleVolume, 1 );
}
