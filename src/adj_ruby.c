
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

VALUE adj_nadj( VALUE self )
{
  GET_ADJ_FROM_SELF;
  return INT2NUM( adjNAdj(adj) );
}

VALUE adj_chunkSize( VALUE self )
{
  GET_ADJ_FROM_SELF;
  return INT2NUM( adjChunkSize(adj) );
}

VALUE adj_realloc( VALUE self, VALUE nnode )
{
  GET_ADJ_FROM_SELF;
  return (adjRealloc( adj, NUM2INT(nnode) )==NULL?Qnil:self);
}

VALUE adj_register( VALUE self, VALUE node, VALUE item )
{
  GET_ADJ_FROM_SELF;
  return (adjRegister( adj, NUM2INT(node), NUM2INT(item) )==NULL?Qnil:self);
}

VALUE adj_remove( VALUE self, VALUE node, VALUE item )
{
  GET_ADJ_FROM_SELF;
  return ( adjRemove( adj, NUM2INT(node), NUM2INT(item) )==NULL?Qnil:self );
}

VALUE adj_valid( VALUE self )
{
  GET_ADJ_FROM_SELF;
  return (adjValid(adjGetCurrent(adj))?Qtrue:Qfalse);
}

VALUE adj_more( VALUE self )
{
  GET_ADJ_FROM_SELF;
  return (adjMore(adjGetCurrent(adj))?Qtrue:Qfalse);
}

VALUE adj_first( VALUE self, VALUE node )
{
  GET_ADJ_FROM_SELF;
  adjSetCurrent( adj, adjFirst(adj, NUM2INT(node)) );
  return self;
}

VALUE adj_item( VALUE self )
{
  GET_ADJ_FROM_SELF;
  return INT2NUM( adjItem(adjGetCurrent(adj)) );
}

VALUE adj_next( VALUE self )
{
  GET_ADJ_FROM_SELF;
  adjSetCurrent(adj,adjNext(adjGetCurrent(adj)));
  return self;
}

VALUE adj_exists( VALUE self, VALUE node, VALUE item )
{
  GET_ADJ_FROM_SELF;
  return( adjExists( adj, NUM2INT(node), NUM2INT(item) )?Qtrue:Qfalse );
}

VALUE adj_degree( VALUE self, VALUE node )
{
  GET_ADJ_FROM_SELF;
  return INT2NUM( adjDegree(adj, NUM2INT(node) ) );
}

VALUE cAdj;

void Init_Adj() 
{
  cAdj = rb_define_class( "Adj", rb_cObject );
  rb_define_singleton_method( cAdj, "new", adj_new, 2 );
  rb_define_method( cAdj, "nnode", adj_nnode, 0 );
  rb_define_method( cAdj, "nadj", adj_nadj, 0 );
  rb_define_method( cAdj, "chunkSize", adj_chunkSize, 0 );
  rb_define_method( cAdj, "realloc", adj_realloc, 1 );
  rb_define_method( cAdj, "register", adj_register, 2 );
  rb_define_method( cAdj, "remove", adj_remove, 2 );
  rb_define_method( cAdj, "valid", adj_valid, 0 );
  rb_define_method( cAdj, "more", adj_more, 0 );
  rb_define_method( cAdj, "first", adj_first, 1 );
  rb_define_method( cAdj, "item", adj_item, 0);
  rb_define_method( cAdj, "next", adj_next, 0);
  rb_define_method( cAdj, "exists", adj_exists, 2 );
  rb_define_method( cAdj, "degree", adj_degree, 1 );
}
