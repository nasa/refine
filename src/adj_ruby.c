
#include "ruby.h"
#include "adj.h"

#define GET_ADJ_FROM_SELF Adj *adj; Data_Get_Struct( self, Adj, adj );

static void adj_free( void *adj )
{
  adjFree( adj );
}

VALUE adj_new( VALUE class, VALUE nnode, VALUE perNode )
{
  Adj *adj;
  VALUE obj;
  adj = adjCreate( NUM2INT(nnode), NUM2INT(perNode) );
  obj = Data_Wrap_Struct( class, 0, adj_free, adj );
  return obj;
}

VALUE adj_nnode( VALUE self )
{
  GET_ADJ_FROM_SELF;
  return INT2NUM( adjNNode(adj) );
}

VALUE cAdj;

void Init_Adj() 
{
  cAdj = rb_define_class( "Adj", rb_cObject );
  rb_define_singleton_method( cAdj, "new", adj_new, 2 );
  rb_define_method( cAdj, "nnode", adj_nnode, 0 );
}
