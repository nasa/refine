
#include "ruby.h"
#include "ref_cell.h"

#define GET_REF_CELL_FROM_SELF REF_CELL ref_cell; Data_Get_Struct( self, REF_CELL, ref_cell );

static void ref_cell_sweeper( void *ref_cell )
{
  ref_cell_free( (REF_CELL)ref_cell );
}

VALUE ref_cell_new( VALUE class, VALUE nnode )
{
  REF_CELL ref_cell;
  VALUE obj;
  ref_cell_create( NUM2INT(nnode), &ref_cell );
  obj = Data_Wrap_Struct( class, 0, ref_cell_sweeper, ref_cell );
  return obj;
}

VALUE ref_cell_class;

void Init_ref_cell()
{
  ref_cell_class = rb_define_class( "Ref_Cell", rb_cObject );
  rb_define_singleton_method( ref_cell_class, "new", ref_cell_new, 1 );
}
