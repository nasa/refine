
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

VALUE layer_nfront( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNFront(layer) );
}

VALUE layer_nnormal( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNNormal(layer) );
}

VALUE layer_maxnode( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerMaxNode(layer) );
}

VALUE layer_makeFront( VALUE self, VALUE rb_bc )
{
  Layer *rLayer;
  int i, nbc, *bc;
  GET_LAYER_FROM_SELF;
  nbc=0;
  while ( Qnil != rb_ary_entry(rb_bc,nbc) ) nbc++;
  bc = malloc( nbc*sizeof(int));
  for (i=0;i<nbc;i++) bc[i] = NUM2INT(rb_ary_entry(rb_bc,i));
  rLayer = layerMakeFront(layer,nbc,bc);
  free(bc);
  return ( layer == rLayer?self:Qnil );
}

VALUE layer_front( VALUE self, VALUE front )
{
  int nodes[3];
  VALUE rb_front;
  GET_LAYER_FROM_SELF;
  if (layer != layerFront(layer, NUM2INT(front), nodes )) return Qnil;
  rb_front = rb_ary_new2(3);
  rb_ary_store( rb_front, 0, INT2NUM(nodes[0]) );
  rb_ary_store( rb_front, 1, INT2NUM(nodes[1]) );
  rb_ary_store( rb_front, 2, INT2NUM(nodes[2]) );
  return rb_front;
}

VALUE layer_frontDirection( VALUE self, VALUE front )
{
  int i;
  double direction[3];
  VALUE rb_direction;
  GET_LAYER_FROM_SELF;
  if (layer == layerFrontDirection(layer,NUM2INT(front),direction)){
    rb_direction = rb_ary_new2(3);
    for ( i=0 ; i < 3 ; i++ ) 
      rb_ary_store( rb_direction, i, rb_float_new(direction[i]) );
  }else{
    rb_direction = Qnil;
  }
  return rb_direction;
}

VALUE layer_makeNormal( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerMakeNormal(layer)?self:Qnil );
}

VALUE layer_frontNormals( VALUE self, VALUE front )
{
  int normals[3];
  VALUE rb_normals;
  GET_LAYER_FROM_SELF;
  if (layer != layerFrontNormals(layer, NUM2INT(front), normals )) return Qnil;
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

VALUE layer_normalFronts( VALUE self, VALUE normal )
{
  int i, nfront, *front;
  VALUE rb_front;
  GET_LAYER_FROM_SELF;
  nfront = layerNormalDeg(layer,NUM2INT(normal));
  if (nfront<=0) return Qnil;
  front = malloc(nfront*sizeof(int));
  if (layer == layerNormalFronts(layer,NUM2INT(normal),nfront,front)){
    rb_front = rb_ary_new2(nfront);
    for ( i=0 ; i < nfront ; i++ ) 
      rb_ary_store( rb_front, i, INT2NUM(front[i]) );
  }else{
    rb_front = Qnil;
  }
  free(front);
  return rb_front;
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

VALUE layer_constrainNormal( VALUE self, VALUE bc )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerConstrainNormal(layer,NUM2INT(bc))?self:Qnil );
}

VALUE layer_constrained( VALUE self, VALUE normal )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM(layerConstrained(layer,NUM2INT(normal)));
}

VALUE layer_advance( VALUE self, VALUE height )
{
  GET_LAYER_FROM_SELF;
  return ( layer == layerAdvance(layer,NUM2DBL(height))?self:Qnil );
}

VALUE cLayer;

void Init_Layer() 
{
  cLayer = rb_define_class( "Layer", rb_cObject );
  rb_define_singleton_method( cLayer, "new", layer_new, 1 );
  rb_define_method( cLayer, "nfront", layer_nfront, 0 );
  rb_define_method( cLayer, "nnormal", layer_nnormal, 0 );
  rb_define_method( cLayer, "maxnode", layer_maxnode, 0 );
  rb_define_method( cLayer, "makeFront", layer_makeFront, 1 );
  rb_define_method( cLayer, "front", layer_front, 1 );
  rb_define_method( cLayer, "frontDirection", layer_frontDirection, 1 );
  rb_define_method( cLayer, "makeNormal", layer_makeNormal, 0 );
  rb_define_method( cLayer, "frontNormals", layer_frontNormals, 1 );
  rb_define_method( cLayer, "normalRoot", layer_normalRoot, 1 );
  rb_define_method( cLayer, "normalDeg", layer_normalDeg, 1 );
  rb_define_method( cLayer, "normalFronts", layer_normalFronts, 1 );
  rb_define_method( cLayer, "normalDirection", layer_normalDirection, 1 );
  rb_define_method( cLayer, "constrainNormal", layer_constrainNormal, 1 );
  rb_define_method( cLayer, "constrained", layer_constrained, 1 );
  rb_define_method( cLayer, "advance", layer_advance, 1 );
}
