
#include "ruby.h"
#include "grid.h"
#include "gridedger.h"

#define GET_GE_FROM_SELF GridEdger *ge; Data_Get_Struct(self, GridEdger, ge);

static void gridedger_mark( void *voidGridEdger )
{
  GridEdger *ge = (GridEdger *)voidGridEdger;
  VALUE grid =(VALUE)(ge->gridRubyVALUEusedForGC); 
  rb_gc_mark(grid);
}

static void gridedger_free( void *voidGridEdger )
{
  GridEdger *ge = (GridEdger *)voidGridEdger;
  gridedgerFree( ge );
}

static VALUE gridedger_new( VALUE class, VALUE rb_grid, VALUE edgeId )
{
  GridEdger *ge;
  Grid *grid;
  VALUE obj;
  Data_Get_Struct(rb_grid, Grid, grid);
  ge = gridedgerCreate( grid, NUM2INT( edgeId ) );
  ge->gridRubyVALUEusedForGC = (void *)rb_grid;
  obj = Data_Wrap_Struct( class, gridedger_mark, gridedger_free, ge );
  return obj;
}

VALUE gridedger_edgeId( VALUE self )
{
  GET_GE_FROM_SELF;
  return INT2NUM( gridedgerEdgeId( ge ) );
}

VALUE gridedger_idealNodes( VALUE self )
{
  GET_GE_FROM_SELF;
  return INT2NUM( gridedgerIdealNodes( ge ) );
}

VALUE gridedger_idealNodeT( VALUE self, VALUE rb_node )
{
  int node;
  double s;
  GET_GE_FROM_SELF;
  node = NUM2INT( rb_node );
  return (ge == gridedgerIdealNodeT( ge, node, &s )?rb_float_new(s):Qnil);
}

VALUE gridedger_discreteSegmentAndRatio( VALUE self, VALUE rb_segment )
{
  double segment;
  int discrete_segment;
  double segment_ratio;
  VALUE rb_result;
  GET_GE_FROM_SELF;
  segment = NUM2DBL( rb_segment );
  if ( ge != gridedgerDiscreteSegmentAndRatio( ge, segment, 
					       &discrete_segment, 
					       &segment_ratio ) ) return Qnil;
  rb_result = rb_ary_new2(2);
  rb_ary_store( rb_result, 0, INT2NUM(discrete_segment) );
  rb_ary_store( rb_result, 1, rb_float_new(segment_ratio) );
  return rb_result;
}

VALUE gridedger_supportingSegment( VALUE self, VALUE rb_t )
{
  double t, segment;
  GET_GE_FROM_SELF;
  t = NUM2DBL( rb_t );
  return (ge == gridedgerSupportingSegment( ge, t, &segment )?rb_float_new(segment):Qnil);
}

VALUE gridedger_segmentT( VALUE self, VALUE rb_segment )
{
  double t, segment;
  GET_GE_FROM_SELF;
  segment = NUM2DBL( rb_segment );
  return (ge == gridedgerSegmentT( ge, segment, &t )?rb_float_new(t):Qnil);
}

VALUE gridedger_segmentMap( VALUE self, VALUE rb_segment )
{
  double segment;
  int i;
  double map[6];
  VALUE rb_map;
  GET_GE_FROM_SELF;
  segment = NUM2DBL( rb_segment );
  if ( ge != gridedgerSegmentMap( ge, segment, map ) ) return Qnil;
  rb_map = rb_ary_new2(6);
  for(i=0;i<6;i++) rb_ary_store( rb_map, i, rb_float_new(map[i]) );
  return rb_map;
}

VALUE gridedger_lengthToS( VALUE self, VALUE rb_segment, VALUE rb_length )
{
  double segment, length, s;
  GET_GE_FROM_SELF;
  segment = NUM2DBL( rb_segment );
  length  = NUM2DBL( rb_length );
  return (ge == gridedgerLengthToS( ge, segment, length, &s )?rb_float_new(s):Qnil);
}

VALUE gridedger_discretize( VALUE self, VALUE rb_length )
{
  double length;
  GET_GE_FROM_SELF;
  length = NUM2DBL( rb_length );
  return (ge == gridedgerDiscretize( ge, length )?self:Qnil);
}

VALUE gridedger_discretizeEvenly( VALUE self )
{
  GET_GE_FROM_SELF;
  return (ge == gridedgerDiscretizeEvenly( ge )?self:Qnil);
}

VALUE gridedger_discretizeOnce( VALUE self )
{
  GET_GE_FROM_SELF;
  return (ge == gridedgerDiscretizeOnce( ge )?self:Qnil);
}

VALUE gridedger_insert( VALUE self )
{
  GET_GE_FROM_SELF;
  return (ge == gridedgerInsert( ge )?self:Qnil);
}

VALUE gridedger_removeUnused( VALUE self )
{
  GET_GE_FROM_SELF;
  return (ge == gridedgerRemoveUnused( ge )?self:Qnil);
}

VALUE cGridEdger;

void Init_GridEdger() 
{
  cGridEdger = rb_define_class( "GridEdger", rb_cObject );
  rb_define_singleton_method( cGridEdger, "new", gridedger_new, 2 );
  rb_define_method( cGridEdger, "edgeId", gridedger_edgeId, 0 );
  rb_define_method( cGridEdger, "idealNodes", gridedger_idealNodes, 0 );
  rb_define_method( cGridEdger, "idealNodeT", gridedger_idealNodeT, 1 );
  rb_define_method( cGridEdger, "supportingSegment", 
		    gridedger_supportingSegment, 1 );
  rb_define_method( cGridEdger, "discreteSegmentAndRatio", 
		    gridedger_discreteSegmentAndRatio, 1 );
  rb_define_method( cGridEdger, "segmentT", gridedger_segmentT, 1 );
  rb_define_method( cGridEdger, "segmentMap", gridedger_segmentMap, 1 );
  rb_define_method( cGridEdger, "lengthToS", gridedger_lengthToS, 2 );
  rb_define_method( cGridEdger, "discretize", gridedger_discretize, 1 );
  rb_define_method( cGridEdger, "discretizeEvenly", 
		    gridedger_discretizeEvenly, 0 );
  rb_define_method( cGridEdger, "discretizeOnce", 
		    gridedger_discretizeOnce, 0 );
  rb_define_method( cGridEdger, "insert", gridedger_insert, 0 );
  rb_define_method( cGridEdger, "removeUnused", gridedger_removeUnused, 0 );
}
