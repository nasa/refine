
#include "ruby.h"
#include "layer.h"

#define GET_LAYER_FROM_SELF Layer *layer; Data_Get_Struct(self, Layer, layer);

static void layer_free( void *layer )
{
  layerFree( layer );
}

VALUE layer_new( VALUE class )
{
  Layer *layer;
  VALUE obj;
  layer = layerCreate( );
  obj = Data_Wrap_Struct( class, 0, layer_free, layer );
  return obj;
}

VALUE layer_nfront( VALUE self )
{
  GET_LAYER_FROM_SELF;
  return INT2NUM( layerNFront(layer) );
}

VALUE cLayer;

void Init_Layer() 
{
  cLayer = rb_define_class( "Layer", rb_cObject );
  rb_define_singleton_method( cLayer, "new", layer_new, 0 );
  rb_define_method( cLayer, "nfront", layer_nfront, 0 );
}
