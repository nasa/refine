
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

static VALUE gridedger_new( VALUE class, VALUE rb_grid )
{
  GridEdger *ge;
  Grid *grid;
  VALUE obj;
  Data_Get_Struct(rb_grid, Grid, grid);
  ge = gridedgerCreate( grid );
  ge->gridRubyVALUEusedForGC = (void *)rb_grid;
  obj = Data_Wrap_Struct( class, gridedger_mark, gridedger_free, ge );
  return obj;
}

VALUE cGridEdger;

void Init_GridEdger() 
{
  cGridEdger = rb_define_class( "GridEdger", rb_cObject );
  rb_define_singleton_method( cGridEdger, "new", gridedger_new, 1 );
}
