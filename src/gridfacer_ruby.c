
#include "ruby.h"
#include "grid.h"
#include "gridfacer.h"

#define GET_GF_FROM_SELF GridFacer *gf; Data_Get_Struct(self, GridFacer, gf);

static void gridfacer_mark( void *voidGridFacer )
{
  GridFacer *gf = (GridFacer *)voidGridFacer;
  VALUE grid =(VALUE)(gf->gridRubyVALUEusedForGC); 
  rb_gc_mark(grid);
}

static void gridfacer_free( void *voidGridFacer )
{
  GridFacer *gf = (GridFacer *)voidGridFacer;
  gridfacerFree( gf );
}

static VALUE gridfacer_new( VALUE class, VALUE rb_grid, VALUE faceId )
{
  GridFacer *gf;
  Grid *grid;
  VALUE obj;
  Data_Get_Struct(rb_grid, Grid, grid);
  gf = gridfacerCreate( grid, NUM2INT( faceId ) );
  gf->gridRubyVALUEusedForGC = (void *)rb_grid;
  obj = Data_Wrap_Struct( class, gridfacer_mark, gridfacer_free, gf );
  return obj;
}

VALUE gridfacer_faceId( VALUE self )
{
  GET_GF_FROM_SELF;
  return INT2NUM( gridfacerFaceId( gf ) );
}

VALUE cGridFacer;

void Init_GridFacer() 
{
  cGridFacer = rb_define_class( "GridFacer", rb_cObject );
  rb_define_singleton_method( cGridFacer, "new", gridfacer_new, 2 );
  rb_define_method( cGridFacer, "faceId", gridfacer_faceId, 0 );
}
