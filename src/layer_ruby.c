
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
  for (i=0;i<nbc;i++) bc[i] = NUM2INT(rb_ary_entry(rb_bc,1));
  rLayer = layerMakeFront(layer,nbc,bc);
  free(bc);
  return ( layer == rLayer?self:Qnil );
}

VALUE cLayer;

void Init_Layer() 
{
  cLayer = rb_define_class( "Layer", rb_cObject );
  rb_define_singleton_method( cLayer, "new", layer_new, 1 );
  rb_define_method( cLayer, "nfront", layer_nfront, 0 );
  rb_define_method( cLayer, "maxnode", layer_maxnode, 0 );
  rb_define_method( cLayer, "makeFront", layer_makeFront, 1 );
}
