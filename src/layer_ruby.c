
#include "ruby.h"
#include "grid.h"
#include "layer.h"

#define GET_LAYER_FROM_SELF Layer *layer; Data_Get_Struct(self, Layer, layer);

static void layer_free( void *layer )
{
  layerFree( layer );
}

VALUE layer_new( VALUE class, VALUE rb_grid )
{
  Layer *layer;
  Grid *grid;
  VALUE obj;
  Data_Get_Struct(rb_grid, Grid, grid);
  layer = layerCreate( grid );
  obj = Data_Wrap_Struct( class, 0, layer_free, layer ); // GC mark for grid?
  return obj;
}

VALUE layer_ntriangle( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNTriangle(layer) );
}

VALUE layer_maxtriangle( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerMaxTriangle(layer) );
}

VALUE layer_nblend( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNBlend(layer) );
}

VALUE layer_nnormal( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNNormal(layer) );
}

VALUE layer_maxnormal( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerMaxNormal(layer) );
}

VALUE layer_maxnode( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerMaxNode(layer) );
}

VALUE layer_addTriangle( VALUE self, VALUE n0, VALUE n1, VALUE n2 )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerAddTriangle(layer,NUM2INT(n0),NUM2INT(n1),NUM2INT(n2))?self:Qnil );
}

VALUE layer_populateAdvancingFront( VALUE self, VALUE rb_bc )
{
  Layer *rLayer;
  int i, nbc, *bc;
  GET_LAYER_FROM_SELF;
  nbc=0;
  while ( Qnil != rb_ary_entry(rb_bc,nbc) ) nbc++;
  bc = malloc( nbc*sizeof(int));
  for (i=0;i<nbc;i++) bc[i] = NUM2INT(rb_ary_entry(rb_bc,i));
  rLayer = layerPopulateAdvancingFront(layer,nbc,bc);
  free(bc);
  return ( layer == rLayer?self:Qnil );
}

VALUE layer_triangle( VALUE self, VALUE triangle )
{
  int nodes[3];
  VALUE rb_triangle;
  GET_LAYER_FROM_SELF;
  if (layer != layerTriangle(layer, NUM2INT(triangle), nodes )) return Qnil;
  rb_triangle = rb_ary_new2(3);
  rb_ary_store( rb_triangle, 0, INT2NUM(nodes[0]) );
  rb_ary_store( rb_triangle, 1, INT2NUM(nodes[1]) );
  rb_ary_store( rb_triangle, 2, INT2NUM(nodes[2]) );
  return rb_triangle;
}

VALUE layer_addParentGeomFace( VALUE self, VALUE faceId )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerAddParentGeomFace(layer,NUM2INT(faceId))?self:Qnil );
}

VALUE layer_parentGeomFace( VALUE self, VALUE faceId )
{
  GET_LAYER_FROM_SELF;
  return ( layerParentGeomFace(layer,NUM2INT(faceId))?Qtrue:Qfalse );
}

VALUE layer_triangleDirection( VALUE self, VALUE triangle )
{
  int i;
  double direction[3];
  VALUE rb_direction;
  GET_LAYER_FROM_SELF;
  if (layer == layerTriangleDirection(layer,NUM2INT(triangle),direction)){
    rb_direction = rb_ary_new2(3);
    for ( i=0 ; i < 3 ; i++ ) 
      rb_ary_store( rb_direction, i, rb_float_new(direction[i]) );
  }else{
    rb_direction = Qnil;
  }
  return rb_direction;
}

VALUE layer_addNormal( VALUE self, VALUE globalNodeId )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerAddNormal(layer,NUM2INT(globalNodeId)));
}

VALUE layer_uniqueNormalId( VALUE self, VALUE globalNodeId )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerUniqueNormalId(layer,NUM2INT(globalNodeId)));
}

VALUE layer_triangleNormals( VALUE self, VALUE triangle )
{
  int normals[3];
  VALUE rb_normals;
  GET_LAYER_FROM_SELF;
  if (layer != layerTriangleNormals(layer, NUM2INT(triangle), normals )) return Qnil;
  rb_normals = rb_ary_new2(3);
  rb_ary_store( rb_normals, 0, INT2NUM(normals[0]) );
  rb_ary_store( rb_normals, 1, INT2NUM(normals[1]) );
  rb_ary_store( rb_normals, 2, INT2NUM(normals[2]) );
  return rb_normals;
}

VALUE layer_normalRoot( VALUE self, VALUE normal )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerNormalRoot(layer,NUM2INT(normal)));
}

VALUE layer_normalDeg( VALUE self, VALUE normal )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerNormalDeg(layer,NUM2INT(normal)));
}

VALUE layer_normalTriangles( VALUE self, VALUE normal )
{
  int i, ntriangle, *triangle;
  VALUE rb_triangle;
  GET_LAYER_FROM_SELF;
  ntriangle = layerNormalDeg(layer,NUM2INT(normal));
  if (ntriangle<=0) return Qnil;
  triangle = malloc(ntriangle*sizeof(int));
  if (layer == layerNormalTriangles(layer,NUM2INT(normal),ntriangle,triangle)){
    rb_triangle = rb_ary_new2(ntriangle);
    for ( i=0 ; i < ntriangle ; i++ ) 
      rb_ary_store( rb_triangle, i, INT2NUM(triangle[i]) );
  }else{
    rb_triangle = Qnil;
  }
  free(triangle);
  return rb_triangle;
}

VALUE layer_previousTriangle( VALUE self, VALUE normal, VALUE triangle )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerPreviousTriangle(layer,NUM2INT(normal),NUM2INT(triangle)) );
}

VALUE layer_nextTriangle( VALUE self, VALUE normal, VALUE triangle )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNextTriangle(layer,NUM2INT(normal),NUM2INT(triangle)) );
}

VALUE layer_edgeAngle( VALUE self, VALUE triangle0, VALUE triangle1 )
{
  GET_LAYER_FROM_SELF;
  return rb_float_new( layerEdgeAngle(layer,
				      NUM2INT(triangle0),
				      NUM2INT(triangle1) ) );
}

VALUE layer_normalDirection( VALUE self, VALUE normal )
{
  int i;
  double direction[3];
  VALUE rb_direction;
  GET_LAYER_FROM_SELF;
  if (layer == layerNormalDirection(layer,NUM2INT(normal),direction)){
    rb_direction = rb_ary_new2(3);
    for ( i=0 ; i < 3 ; i++ ) 
      rb_ary_store( rb_direction, i, rb_float_new(direction[i]) );
  }else{
    rb_direction = Qnil;
  }
  return rb_direction;
}

VALUE layer_setNormalHeight( VALUE self, VALUE normal, VALUE height )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerSetNormalHeight( layer,
					  NUM2INT(normal),
					  NUM2DBL(height) )?self:Qnil );
}

VALUE layer_getNormalHeight( VALUE self, VALUE normal )
{
  double height;
  GET_LAYER_FROM_SELF;
  if (layer != layerGetNormalHeight( layer, NUM2INT(normal), &height ) )
    return Qnil;
  return rb_float_new(height);
}

VALUE layer_assignPolynomialNormalHeight( VALUE self, VALUE constant, VALUE slope, VALUE exponent, VALUE rb_origin, VALUE rb_direction )
{
  int i;
  double origin[3], direction[3];
  GET_LAYER_FROM_SELF;
  for(i=0;i<3;i++){
    origin[i]=NUM2DBL(rb_ary_entry(rb_origin,i));
    direction[i]=NUM2DBL(rb_ary_entry(rb_direction,i));
  }
  return ( layer == layerAssignPolynomialNormalHeight( layer, 
                                                   NUM2DBL(constant),
                                                   NUM2DBL(slope),
                                                   NUM2DBL(exponent),
                                                   origin,
                                                   direction)?self:Qnil );
}

VALUE layer_visibleNormals( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerVisibleNormals(layer,-1.0,-1.0)?self:Qnil );
}

VALUE layer_projectNormalsToConstraints( VALUE self  )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerProjectNormalsToConstraints(layer)?self:Qnil );
}

VALUE layer_constrainNormal( VALUE self, VALUE edgeface )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerConstrainNormal(layer,NUM2INT(edgeface))?self:Qnil );
}

VALUE layer_constrainingGeometry( VALUE self, VALUE edgeface )
{
  GET_LAYER_FROM_SELF;
  return ( layerConstrainingGeometry(layer,NUM2INT(edgeface))?Qtrue:Qfalse );
}

VALUE layer_constrained( VALUE self, VALUE normal )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerConstrained(layer,NUM2INT(normal)));
}

VALUE layer_constrainTriangleSide( VALUE self, VALUE normal0, VALUE normal1, 
				VALUE bc )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerConstrainTriangleSide(layer,
					    NUM2INT(normal0),
					    NUM2INT(normal1),
					    NUM2INT(bc) )?self:Qnil );
}

VALUE layer_constrainedSide( VALUE self, VALUE triangle, VALUE side )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerConstrainedSide(layer,NUM2INT(triangle),NUM2INT(side)));
}

VALUE layer_nConstrainedSides( VALUE self, VALUE faceId )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerNConstrainedSides(layer,NUM2INT(faceId)));
}

VALUE layer_findParentGeomEdges( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerFindParentGeomEdges(layer)?self:Qnil );
}

VALUE layer_setParentGeomEdge( VALUE self, VALUE normal0, VALUE normal1, 
			   VALUE edgeId )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerSetParentGeomEdge(layer,
				       NUM2INT(normal0),
				       NUM2INT(normal1),
				       NUM2INT(edgeId) )?self:Qnil );
}

VALUE layer_parentGeomEdge( VALUE self, VALUE triangle, VALUE side )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerParentGeomEdge(layer,NUM2INT(triangle),NUM2INT(side)));
}

VALUE layer_nParentGeomEdgeSegments( VALUE self, VALUE edgeId )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerNParentGeomEdgeSegments(layer,NUM2INT(edgeId)));
}

VALUE layer_terminateNormal( VALUE self, VALUE normal )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerTerminateNormal(layer,NUM2INT(normal))?self:Qnil );
}

VALUE layer_normalTerminated( VALUE self, VALUE normal )
{
  GET_LAYER_FROM_SELF;
  return ( layerNormalTerminated(layer,NUM2INT(normal))?Qtrue:Qfalse );
}

VALUE layer_cellInLayer( VALUE self, VALUE cell )
{
  GET_LAYER_FROM_SELF;
  return ( layerCellInLayer(layer,NUM2INT(cell))?Qtrue:Qfalse );
}

VALUE layer_faceInLayer( VALUE self, VALUE face )
{
  GET_LAYER_FROM_SELF;
  return ( layerFaceInLayer(layer,NUM2INT(face))?Qtrue:Qfalse );
}

VALUE layer_edgeInLayer( VALUE self, VALUE edge )
{
  GET_LAYER_FROM_SELF;
  return ( layerEdgeInLayer(layer,NUM2INT(edge))?Qtrue:Qfalse );
}

VALUE layer_advance( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerAdvance(layer)?self:Qnil );
}

VALUE layer_advanceConstantHeight( VALUE self, VALUE height )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerAdvanceConstantHeight(layer,NUM2DBL(height))?self:Qnil );
}

VALUE layer_wiggle( VALUE self, VALUE height )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerWiggle(layer,NUM2DBL(height))?self:Qnil );
}

VALUE layer_tetrahedraOnly( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layerTetrahedraOnly(layer)?Qtrue:Qfalse);
}

VALUE layer_toggleMixedElementMode( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerToggleMixedElementMode(layer)?self:Qnil );
}

VALUE layer_blend( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerBlend(layer)?self:Qnil );
}

VALUE layer_splitBlend( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerSplitBlend(layer)?self:Qnil );
}

VALUE layer_blendNormals( VALUE self, VALUE blend )
{
  int i, normals[4];
  VALUE rb_normals;
  GET_LAYER_FROM_SELF;
  if (layer != layerBlendNormals(layer, NUM2INT(blend), normals )) return Qnil;
  rb_normals = rb_ary_new2(4);
  for(i=0;i<4;i++) rb_ary_store( rb_normals, i, INT2NUM(normals[i]) );
  return rb_normals;
}

VALUE cLayer;

void Init_Layer() 
{
  cLayer = rb_define_class( "Layer", rb_cObject );
  rb_define_singleton_method( cLayer, "new", layer_new, 1 );
  rb_define_method( cLayer, "ntriangle", layer_ntriangle, 0 );
  rb_define_method( cLayer, "maxtriangle", layer_maxtriangle, 0 );
  rb_define_method( cLayer, "nblend", layer_nblend, 0 );
  rb_define_method( cLayer, "maxnormal", layer_maxnormal, 0 );
  rb_define_method( cLayer, "nnormal", layer_nnormal, 0 );
  rb_define_method( cLayer, "maxnode", layer_maxnode, 0 );
  rb_define_method( cLayer, "addTriangle", layer_addTriangle, 3 );
  rb_define_method( cLayer, "populateAdvancingFront", layer_populateAdvancingFront, 1 );
  rb_define_method( cLayer, "triangle", layer_triangle, 1 );
  rb_define_method( cLayer, "triangleDirection", layer_triangleDirection, 1 );
  rb_define_method( cLayer, "addNormal", layer_addNormal, 1 );
  rb_define_method( cLayer, "uniqueNormalId", layer_uniqueNormalId, 1 );
  rb_define_method( cLayer, "addParentGeomFace", layer_addParentGeomFace, 1 );
  rb_define_method( cLayer, "parentGeomFace", layer_parentGeomFace, 1 );
  rb_define_method( cLayer, "triangleNormals", layer_triangleNormals, 1 );
  rb_define_method( cLayer, "normalRoot", layer_normalRoot, 1 );
  rb_define_method( cLayer, "normalDeg", layer_normalDeg, 1 );
  rb_define_method( cLayer, "normalTriangles", layer_normalTriangles, 1 );
  rb_define_method( cLayer, "previousTriangle", layer_previousTriangle, 2 );
  rb_define_method( cLayer, "nextTriangle", layer_nextTriangle, 2 );
  rb_define_method( cLayer, "edgeAngle", layer_edgeAngle, 2 );
  rb_define_method( cLayer, "normalDirection", layer_normalDirection, 1 );
  rb_define_method( cLayer, "setNormalHeight", layer_setNormalHeight, 2 );
  rb_define_method( cLayer, "getNormalHeight", layer_getNormalHeight, 1 );
  rb_define_method( cLayer, "assignPolynomialNormalHeight",
                    layer_assignPolynomialNormalHeight, 5 );
  rb_define_method( cLayer, "visibleNormals", layer_visibleNormals, 0 );
  rb_define_method( cLayer, "projectNormalsToConstraints", 
		    layer_projectNormalsToConstraints, 0 );
  rb_define_method( cLayer, "constrainNormal", layer_constrainNormal, 1 );
  rb_define_method( cLayer, "constrainingGeometry", layer_constrainingGeometry, 1 );
  rb_define_method( cLayer, "constrained", layer_constrained, 1 );
  rb_define_method( cLayer, "constrainTriangleSide", layer_constrainTriangleSide, 3 );
  rb_define_method( cLayer, "constrainedSide", layer_constrainedSide, 2 );
  rb_define_method( cLayer, "nConstrainedSides", layer_nConstrainedSides, 1 );
  rb_define_method( cLayer, "findParentGeomEdges", layer_findParentGeomEdges, 0 );
  rb_define_method( cLayer, "setParentGeomEdge", layer_setParentGeomEdge, 3 );
  rb_define_method( cLayer, "parentGeomEdge", layer_parentGeomEdge, 2 );
  rb_define_method( cLayer, "nParentGeomEdgeSegments", layer_nParentGeomEdgeSegments, 1 );

  rb_define_method( cLayer, "terminateNormal", layer_terminateNormal, 1 );
  rb_define_method( cLayer, "normalTerminated", layer_normalTerminated, 1 );

  rb_define_method( cLayer, "cellInLayer", layer_cellInLayer, 1 );
  rb_define_method( cLayer, "faceInLayer", layer_faceInLayer, 1 );
  rb_define_method( cLayer, "edgeInLayer", layer_edgeInLayer, 1 );

  rb_define_method( cLayer, "advanceConstantHeight", 
		    layer_advanceConstantHeight, 1 );
  rb_define_method( cLayer, "advance", layer_advance, 0 );
  rb_define_method( cLayer, "wiggle", layer_wiggle, 1 );

  rb_define_method( cLayer, "tetrahedraOnly", layer_tetrahedraOnly, 0 ); 
  rb_define_method( cLayer, "toggleMixedElementMode", layer_toggleMixedElementMode, 0 ); 

  rb_define_method( cLayer, "blend", layer_blend, 0 );
  rb_define_method( cLayer, "splitBlend", layer_splitBlend, 0 );
  rb_define_method( cLayer, "blendNormals", layer_blendNormals, 1 );

}
