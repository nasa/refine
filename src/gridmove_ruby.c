
#include "ruby.h"
#include "grid.h"
#include "gridmove.h"

#define GET_GM_FROM_SELF GridMove *gm; Data_Get_Struct(self, GridMove, gm);

static void gridmove_mark( void *voidGridMove )
{
  GridMove *gm = (GridMove *)voidGridMove;
  VALUE grid =(VALUE)(gm->gridRubyVALUEusedForGC); 
  rb_gc_mark(grid);
}

static void gridmove_free( void *voidGridMove )
{
  GridMove *gm = (GridMove *)voidGridMove;
  gridmoveFree( gm );
}

VALUE gridmove_new( VALUE class, VALUE rb_grid )
{
  GridMove *gm;
  Grid *grid;
  VALUE obj;
  Data_Get_Struct(rb_grid, Grid, grid);
  gm = gridmoveCreate( grid );
  gm->gridRubyVALUEusedForGC = (void *)rb_grid;
  obj = Data_Wrap_Struct( class, gridmove_mark, gridmove_free, gm );
  return obj;
}

VALUE gridmove_displace( VALUE self, VALUE node, VALUE rb_displace )
{
  int i;
  double displace[3];
  GET_GM_FROM_SELF;
  for(i=0;i<3;i++) displace[i]=NUM2DBL(rb_ary_entry(rb_displace,i));
  return ( gm == gridmoveDisplace( gm, NUM2INT(node), displace)?self:Qnil );
}

VALUE gridmove_displacement( VALUE self, VALUE node )
{
  int i;
  double displacement[3];
  VALUE rb_displacement;
  GET_GM_FROM_SELF;
  if (gm == gridmoveDisplacement(gm,NUM2INT(node),displacement)){
    rb_displacement = rb_ary_new2(3);
    for ( i=0 ; i < 3 ; i++ ) 
      rb_ary_store( rb_displacement, i, rb_float_new(displacement[i]) );
  }else{
    rb_displacement = Qnil;
  }
  return rb_displacement;
}

VALUE gridmove_specified( VALUE self, VALUE node )
{
  GET_GM_FROM_SELF;
  return ( gridmoveSpecified( gm, NUM2INT(node) )?Qtrue:Qfalse );
}

VALUE gridmove_springRelaxation( VALUE self, VALUE nsteps, VALUE subIterations )
{
  GET_GM_FROM_SELF;
  return ( gm == gridmoveSpringRelaxation( gm, NUM2INT(nsteps), 
					   NUM2INT(subIterations) )?self:Qnil );
}

VALUE gridmove_springs( VALUE self )
{
  VALUE rb_springs;
  int i, nsprings, *springs;
  GET_GM_FROM_SELF;
  if ( gm == gridmoveSprings( gm, &nsprings, &springs ) ) {
    rb_springs = rb_ary_new2(nsprings*2);
    for ( i=0 ; i < nsprings*2 ; i++ ) 
      rb_ary_store( rb_springs, i, INT2NUM(springs[i]) );
    free(springs);
  }else{
    rb_springs = Qnil;
  }
  return rb_springs;
}

VALUE gridmove_applyDisplacements( VALUE self )
{
  GET_GM_FROM_SELF;
  return ( gm == gridmoveApplyDisplacements( gm )?self:Qnil );
}

static VALUE rb_ary_length(VALUE ary)
{
  return INT2NUM(RARRAY(ary)->len);
}

VALUE gridmove_cellFaceNormal( VALUE self, VALUE rb_xyz, VALUE rb_nodes, 
			       VALUE normalIndex )
{
  int len, i;
  double *xyz;
  int nodes[4];
  double normals[4][3];
  VALUE rb_normal;
  GET_GM_FROM_SELF;
  len = NUM2INT(rb_ary_length(rb_xyz));
  xyz = malloc(len*sizeof(double));
  for(i=0;i<len;i++) xyz[i] = NUM2DBL(rb_ary_entry(rb_xyz, i));
  for(i=0;i<4;i++) nodes[i] = NUM2INT(rb_ary_entry(rb_nodes, i));
  if ( gm == gridmoveCellFaceNormals( gm, xyz, nodes, normals ) ) {
    rb_normal = rb_ary_new2(3);
    for ( i=0 ; i < 3 ; i++ ) 
      rb_ary_store( rb_normal, i, 
		    rb_float_new( normals[NUM2INT(normalIndex)][i]) );
  }else{
    rb_normal = Qnil;
  }
  free(xyz);
  return rb_normal;
}

static VALUE gridmove_rowStart(VALUE self, VALUE row)
{
  GET_GM_FROM_SELF;
  return INT2NUM(gridmoveRowStart(gm,NUM2INT(row)));
}

static VALUE gridmove_nnz(VALUE self)
{
  GET_GM_FROM_SELF;
  return INT2NUM(gridmoveNNZ(gm));
}

VALUE gridmove_rowNodes( VALUE self, VALUE rb_row )
{
  int row, length, start, end;
  int entry;
  VALUE rb_nodes;
  GET_GM_FROM_SELF;
  row = NUM2INT(rb_row);
  start = gridmoveRowStart(gm,row);
  end = gridmoveRowStart(gm,row+1);
  if ( EMPTY == start || EMPTY == end ) return Qnil;
  length = end - start;  
  rb_nodes = rb_ary_new2(length);
  for (entry = 0 ; entry < length ; entry++)
    rb_ary_store( rb_nodes, entry, INT2NUM(gridmoveRowNode(gm,start+entry)));
  return rb_nodes;
}

static VALUE gridmove_rowEntry(VALUE self, VALUE row, VALUE node)
{
  GET_GM_FROM_SELF;
  return INT2NUM(gridmoveRowEntry(gm,NUM2INT(row),NUM2INT(node)));
}

VALUE cGridMove;

void Init_GridMove() 
{
  cGridMove = rb_define_class( "GridMove", rb_cObject );
  rb_define_singleton_method( cGridMove, "new", gridmove_new, 1 );
  rb_define_method( cGridMove, "displace", gridmove_displace, 2 );
  rb_define_method( cGridMove, "displacement", gridmove_displacement, 1 );
  rb_define_method( cGridMove, "specified", gridmove_specified, 1 );
  rb_define_method( cGridMove, "springRelaxation", gridmove_springRelaxation, 2 );
  rb_define_method( cGridMove, "springs", gridmove_springs, 0 );
  rb_define_method( cGridMove, "applyDisplacements", gridmove_applyDisplacements, 0 );

  rb_define_method( cGridMove, "cellFaceNormal", gridmove_cellFaceNormal, 3 );
  
  rb_define_method( cGridMove, "rowStart", gridmove_rowStart, 1 );
  rb_define_method( cGridMove, "nnz", gridmove_nnz, 0 );
  rb_define_method( cGridMove, "rowNodes", gridmove_rowNodes, 1 );
  rb_define_method( cGridMove, "rowEntry", gridmove_rowEntry, 2 );
}
