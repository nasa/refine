
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
  obj = Data_Wrap_Struct( class, 0, layer_free, layer );
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

VALUE cLayer;

void Init_Layer() 
{
  cLayer = rb_define_class( "Layer", rb_cObject );
  rb_define_singleton_method( cLayer, "new", layer_new, 1 );
  rb_define_method( cLayer, "nfront", layer_nfront, 0 );
  rb_define_method( cLayer, "maxnode", layer_maxnode, 0 );
}
