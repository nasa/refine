
#include "ruby.h"
#include "gridmetric.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_edgeLength( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridEdgeLength( grid, NUM2INT(n0), NUM2INT(n1) ) );
}

VALUE grid_edgeRatio( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridEdgeRatio( grid, NUM2INT(n0), NUM2INT(n1) ) );
}

VALUE grid_averageEdgeLength( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridAverageEdgeLength( grid, NUM2INT(node) ) );
}

VALUE grid_largestRatioEdge( VALUE self, VALUE node )
{
  double ratio;
  int edgeNode;
  Grid *rGrid;
  GET_GRID_FROM_SELF;
  rGrid = gridLargestRatioEdge( grid, NUM2INT(node), &edgeNode, &ratio );
  return INT2NUM(edgeNode);
}

VALUE grid_smallestRatioEdge( VALUE self, VALUE node )
{
  double ratio;
  int edgeNode;
  Grid *rGrid;
  GET_GRID_FROM_SELF;
  rGrid = gridSmallestRatioEdge( grid, NUM2INT(node), &edgeNode, &ratio );
  return INT2NUM(edgeNode);
}

VALUE grid_spacing( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridSpacing( grid, NUM2INT(node) ) );
}

VALUE grid_resetSpacing( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridResetSpacing(grid)==grid?self:Qnil);
}

VALUE grid_scaleSpacing( VALUE self, VALUE node, VALUE scale )
{
  GET_GRID_FROM_SELF;
  return 
    (gridScaleSpacing(grid, NUM2INT(node), NUM2DBL(scale))==grid?self:Qnil);
}

VALUE grid_scaleSpacingSphere( VALUE self, 
			       VALUE x, VALUE y, VALUE z, VALUE r, 
			       VALUE scale )
{
  GET_GRID_FROM_SELF;
  return
    (gridScaleSpacingSphere(grid, 
			    NUM2DBL(x), NUM2DBL(y), NUM2DBL(z), NUM2DBL(r), 
			    NUM2DBL(scale))==grid?self:Qnil);
}

VALUE grid_setMap( VALUE self, VALUE node, 
		   VALUE m11, VALUE m12, VALUE m13,
		              VALUE m22, VALUE m23,
		                         VALUE m33)
{
  GET_GRID_FROM_SELF;
  return
    (gridSetMap(grid, NUM2INT(node), 
		NUM2DBL(m11), NUM2DBL(m12), NUM2DBL(m13), 
		              NUM2DBL(m22), NUM2DBL(m23), 
		                            NUM2DBL(m33) )==grid?self:Qnil);
}

VALUE grid_copySpacing( VALUE self, VALUE originalNode, VALUE newNode ) 
{
  GET_GRID_FROM_SELF;
  return (gridCopySpacing(grid, 
			  NUM2INT(originalNode), 
			  NUM2INT(newNode) )==grid?self:Qnil);
}

VALUE grid_setMapMatrixToAverageOfNodes( VALUE self, VALUE avgNode, 
					 VALUE node0, VALUE node1 ) 
{
  GET_GRID_FROM_SELF;
  return (gridSetMapMatrixToAverageOfNodes(grid, 
					   NUM2INT(avgNode), 
					   NUM2INT(node0), 
					   NUM2INT(node1) )==grid?self:Qnil);
}


VALUE grid_eigenSystem( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], eigenValues[3], v1[3], v2[3], v3[3];
  VALUE rb_eigenValues, rb_v1, rb_v2, rb_v3, rb_eigenSystem;
  Grid *rGrid;
  GET_GRID_FROM_SELF;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  rGrid = gridEigenSystem( grid, m, eigenValues, v1, v2, v3 );
  if ( rGrid == grid ){
    rb_eigenValues = rb_ary_new2(3);
    rb_v1 = rb_ary_new2(3);
    rb_v2 = rb_ary_new2(3);
    rb_v3 = rb_ary_new2(3);
    for(i=0;i<3;i++){
      rb_ary_store( rb_eigenValues, i, rb_float_new(eigenValues[i]) );
      rb_ary_store( rb_v1, i, rb_float_new(v1[i]) );
      rb_ary_store( rb_v2, i, rb_float_new(v2[i]) );
      rb_ary_store( rb_v3, i, rb_float_new(v3[i]) );
    }
    rb_eigenSystem = rb_ary_new2(4);
    rb_ary_store( rb_eigenSystem, 0, rb_eigenValues );
    rb_ary_store( rb_eigenSystem, 1, rb_v1 );
    rb_ary_store( rb_eigenSystem, 2, rb_v2 );
    rb_ary_store( rb_eigenSystem, 3, rb_v3 );
  }else{
    rb_eigenSystem = Qnil;
  }
  return rb_eigenSystem;
}

VALUE grid_convertMetricToJacobian( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], j[9];
  VALUE rb_j;
  Grid *rGrid;
  GET_GRID_FROM_SELF;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  rGrid = gridConvertMetricToJacobian( grid, m, j );
  if ( rGrid == grid ){
    rb_j = rb_ary_new2(9);
    for(i=0;i<9;i++) rb_ary_store( rb_j, i, rb_float_new(j[i]) );
  }else{
    rb_j = Qnil;
  }
  return rb_j;
}

VALUE grid_volume( VALUE self, VALUE rb_nodes )
{
  int i, nodes[4];
  GET_GRID_FROM_SELF;
  for ( i=0 ; i<4 ; i++ ) nodes[i] = NUM2INT(rb_ary_entry(rb_nodes,i));
  return rb_float_new( gridVolume( grid, nodes ) );
}

VALUE grid_ar( VALUE self, VALUE rb_nodes )
{
  int i, nodes[4];
  GET_GRID_FROM_SELF;
  for ( i=0 ; i<4 ; i++ ) nodes[i] = NUM2INT(rb_ary_entry(rb_nodes,i));
  return rb_float_new( gridAR( grid, nodes ) );
}

VALUE grid_nodeAR( VALUE self, VALUE node )
{
  double ar;
  GET_GRID_FROM_SELF;
  return (gridNodeAR( grid, NUM2INT(node), &ar )==grid?rb_float_new(ar):Qnil);
}

VALUE grid_cellARDerivative( VALUE self, VALUE rb_nodes )
{
  int i, nodes[4];
  double ar, dARdx[3];
  VALUE rb_ar;
  Grid *returnedGrid;
  GET_GRID_FROM_SELF;
  for ( i=0 ; i<4 ; i++ ) nodes[i] = NUM2INT(rb_ary_entry(rb_nodes,i));
  returnedGrid = gridCellARDerivative( grid, nodes, &ar, dARdx );
  if ( returnedGrid == grid ){
    rb_ar = rb_ary_new2(4);
    rb_ary_store( rb_ar, 0, rb_float_new(ar) );
    rb_ary_store( rb_ar, 1, rb_float_new(dARdx[0]) );
    rb_ary_store( rb_ar, 2, rb_float_new(dARdx[1]) );
    rb_ary_store( rb_ar, 3, rb_float_new(dARdx[2]) );
  }else{
    rb_ar = Qnil;
  }
  return rb_ar;
}

VALUE grid_nodeARDerivative( VALUE self, VALUE node )
{
  double ar, dARdx[3];
  VALUE rb_ar;
  Grid *returnedGrid;
  GET_GRID_FROM_SELF;
  returnedGrid = gridNodeARDerivative( grid, NUM2INT(node), &ar, dARdx );
  if ( returnedGrid == grid ){
    rb_ar = rb_ary_new2(4);
    rb_ary_store( rb_ar, 0, rb_float_new(ar) );
    rb_ary_store( rb_ar, 1, rb_float_new(dARdx[0]) );
    rb_ary_store( rb_ar, 2, rb_float_new(dARdx[1]) );
    rb_ary_store( rb_ar, 3, rb_float_new(dARdx[2]) );
  }else{
    rb_ar = Qnil;
  }
  return rb_ar;
}

VALUE grid_storeARDerivative( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (grid == gridStoreARDerivative( grid, NUM2INT(node) )?self:Qnil);
}

VALUE grid_storedARDegree( VALUE self )
{
  GET_GRID_FROM_SELF;
  return INT2NUM(gridStoredARDegree( grid ));
}

VALUE grid_minVolume( VALUE self )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridMinVolume( grid ) );
}

VALUE grid_minAR( VALUE self )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridMinAR( grid ) );
}

VALUE grid_rightHandedFace( VALUE self, VALUE face )
{
  GET_GRID_FROM_SELF;
  return (gridRightHandedFace(grid, NUM2INT(face))?Qtrue:Qfalse);
}

VALUE grid_rightHandedBoundary( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridRightHandedBoundary(grid)?Qtrue:Qfalse);
}

VALUE grid_faceArea( VALUE self, VALUE n0, VALUE n1, VALUE n2 )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridFaceArea(grid, 
				    NUM2INT(n0), NUM2INT(n1), NUM2INT(n2) ) );
}

VALUE grid_faceAR( VALUE self, VALUE n0, VALUE n1, VALUE n2 )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridFaceAR(grid, NUM2INT(n0), NUM2INT(n1), NUM2INT(n2) ) );
}

VALUE grid_minFaceMR( VALUE self )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridMinFaceMR(grid) );
}

VALUE grid_faceMR( VALUE self, VALUE n0, VALUE n1, VALUE n2 )
{
  GET_GRID_FROM_SELF;
  return rb_float_new( gridFaceMR(grid, NUM2INT(n0), NUM2INT(n1), NUM2INT(n2) ) );
}

VALUE grid_faceMRDerivative( VALUE self, VALUE rb_nodes )
{
  int i, nodes[3];
  double mr, dMRdx[3];
  VALUE rb_mr;
  Grid *rGrid;
  GET_GRID_FROM_SELF;
  for ( i=0 ; i<3 ; i++ ) nodes[i] = NUM2INT(rb_ary_entry(rb_nodes,i));
  rGrid = gridFaceMRDerivative( grid, nodes, &mr, dMRdx );
  if ( rGrid == grid ){
    rb_mr = rb_ary_new2(4);
    rb_ary_store( rb_mr, 0, rb_float_new(mr) );
    rb_ary_store( rb_mr, 1, rb_float_new(dMRdx[0]) );
    rb_ary_store( rb_mr, 2, rb_float_new(dMRdx[1]) );
    rb_ary_store( rb_mr, 3, rb_float_new(dMRdx[2]) );
  }else{
    rb_mr = Qnil;
  }
  return rb_mr;
}

VALUE grid_FaceMRDerivative( VALUE self, 
			     VALUE x1, VALUE y1, VALUE z1,
			     VALUE x2, VALUE y2, VALUE z2,
			     VALUE x3, VALUE y3, VALUE z3)
{
  double mr, dMRdx[3];
  VALUE rb_mr;
  GET_GRID_FROM_SELF;
  
  FaceMRDerivative(NUM2DBL(x1), NUM2DBL(y1), NUM2DBL(z1),
		   NUM2DBL(x2), NUM2DBL(y2), NUM2DBL(z2),
		   NUM2DBL(x3), NUM2DBL(y3), NUM2DBL(z3), 
		   &mr, dMRdx );

  rb_mr = rb_ary_new2(4);
  rb_ary_store( rb_mr, 0, rb_float_new(mr) );
  rb_ary_store( rb_mr, 1, rb_float_new(dMRdx[0]) );
  rb_ary_store( rb_mr, 2, rb_float_new(dMRdx[1]) );
  rb_ary_store( rb_mr, 3, rb_float_new(dMRdx[2]) );
  return rb_mr;
}

VALUE grid_nodeFaceMR( VALUE self, VALUE node )
{
  double ar;
  GET_GRID_FROM_SELF;
  return (gridNodeFaceMR( grid, NUM2INT(node), &ar )==grid?rb_float_new(ar):Qnil);
}

VALUE grid_nodeFaceMRDerivative( VALUE self, VALUE node )
{
  double mr, dMRdx[3];
  VALUE rb_mr;
  Grid *rGrid;
  GET_GRID_FROM_SELF;
  rGrid = gridNodeFaceMRDerivative( grid, NUM2INT(node), &mr, dMRdx );
  if ( rGrid == grid ){
    rb_mr = rb_ary_new2(4);
    rb_ary_store( rb_mr, 0, rb_float_new(mr) );
    rb_ary_store( rb_mr, 1, rb_float_new(dMRdx[0]) );
    rb_ary_store( rb_mr, 2, rb_float_new(dMRdx[1]) );
    rb_ary_store( rb_mr, 3, rb_float_new(dMRdx[2]) );
  }else{
    rb_mr = Qnil;
  }
  return rb_mr;
}

VALUE grid_cellMeanRatio( VALUE self, VALUE rb_n0, VALUE rb_n1,  
			  VALUE rb_n2, VALUE rb_n3 )
{
  int i;
  double n0[3], n1[3], n2[3], n3[3];
  GET_GRID_FROM_SELF;
  for (i=0;i<3;i++) {
    n0[i] = NUM2DBL(rb_ary_entry(rb_n0,i));
    n1[i] = NUM2DBL(rb_ary_entry(rb_n1,i));
    n2[i] = NUM2DBL(rb_ary_entry(rb_n2,i));
    n3[i] = NUM2DBL(rb_ary_entry(rb_n3,i));
  }
  return rb_float_new(gridCellMeanRatio(n0,n1,n2,n3));
}

VALUE grid_cellMeanRatioDerivative( VALUE self, VALUE rb_n0, VALUE rb_n1,  
				    VALUE rb_n2, VALUE rb_n3 )
{
  int i;
  double n0[3], n1[3], n2[3], n3[3], mr, dMRdx[3];
  VALUE rb_mr;
  GET_GRID_FROM_SELF;
  for (i=0;i<3;i++) {
    n0[i] = NUM2DBL(rb_ary_entry(rb_n0,i));
    n1[i] = NUM2DBL(rb_ary_entry(rb_n1,i));
    n2[i] = NUM2DBL(rb_ary_entry(rb_n2,i));
    n3[i] = NUM2DBL(rb_ary_entry(rb_n3,i));
  }
  gridCellMeanRatioDerivative(n0,n1,n2,n3,&mr,dMRdx);
  rb_mr = rb_ary_new2(4);
  rb_ary_store(rb_mr,0,rb_float_new(mr));
  for (i=0;i<3;i++) rb_ary_store(rb_mr,i+1,rb_float_new(dMRdx[i]));
  return rb_mr;
}

VALUE cGridMetric;

void Init_GridMetric() 
{
  cGridMetric = rb_define_module( "GridMetric" );
  rb_define_method( cGridMetric, "edgeLength", grid_edgeLength, 2 );
  rb_define_method( cGridMetric, "edgeRatio", grid_edgeRatio, 2 );
  rb_define_method( cGridMetric, "averageEdgeLength", grid_averageEdgeLength, 1 );
  rb_define_method( cGridMetric, "largestRatioEdge", grid_largestRatioEdge, 1 );  rb_define_method( cGridMetric, "smallestRatioEdge", grid_smallestRatioEdge, 1 );
  rb_define_method( cGridMetric, "spacing", grid_spacing, 1 );
  rb_define_method( cGridMetric, "resetSpacing", grid_resetSpacing, 0 );
  rb_define_method( cGridMetric, "scaleSpacing", grid_scaleSpacing, 2 );
  rb_define_method( cGridMetric, "scaleSpacingSphere", grid_scaleSpacingSphere, 5 );
  rb_define_method( cGridMetric, "copySpacing", grid_copySpacing, 2 );
  rb_define_method( cGridMetric, "setMapMatrixToAverageOfNodes", 
		    grid_setMapMatrixToAverageOfNodes, 3 );

  rb_define_method( cGridMetric, "eigenSystem", grid_eigenSystem, 1 );
  rb_define_method( cGridMetric, "convertMetricToJacobian", 
		    grid_convertMetricToJacobian, 1 );
  rb_define_method( cGridMetric, "volume", grid_volume, 1 );
  rb_define_method( cGridMetric, "ar", grid_ar, 1 );
  rb_define_method( cGridMetric, "nodeAR", grid_nodeAR, 1 );
  rb_define_method( cGridMetric, "cellARDerivative", grid_cellARDerivative, 1 );  rb_define_method( cGridMetric, "nodeARDerivative", grid_nodeARDerivative, 1 );  rb_define_method( cGridMetric, "storeARDerivative", grid_storeARDerivative, 1 );
  rb_define_method( cGridMetric, "storedARDegree", grid_storedARDegree, 0 );
  rb_define_method( cGridMetric, "minVolume", grid_minVolume, 0 );
  rb_define_method( cGridMetric, "minAR", grid_minAR, 0 );

  rb_define_method( cGridMetric, "rightHandedFace", grid_rightHandedFace, 1 );
  rb_define_method( cGridMetric, "rightHandedBoundary", grid_rightHandedBoundary, 0);
  rb_define_method( cGridMetric, "faceArea", grid_faceArea, 3);
  rb_define_method( cGridMetric, "faceAR", grid_faceAR, 3);
  rb_define_method( cGridMetric, "faceMR", grid_faceMR, 3);
  rb_define_method( cGridMetric, "minFaceMR", grid_minFaceMR, 0);
  rb_define_method( cGridMetric, "faceMRDerivative", grid_faceMRDerivative, 1);
  rb_define_method( cGridMetric, "FaceMRDerivative", grid_FaceMRDerivative, 9);
  rb_define_method( cGridMetric, "nodeFaceMR", grid_nodeFaceMR, 1 );
  rb_define_method( cGridMetric, "nodeFaceMRDerivative", grid_nodeFaceMRDerivative, 1);
  rb_define_method( cGridMetric, "cellMeanRatio", grid_cellMeanRatio, 4 );
  rb_define_method( cGridMetric, "cellMeanRatioDerivative", grid_cellMeanRatioDerivative, 4 );
}
