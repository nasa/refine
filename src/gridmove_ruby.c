
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

VALUE gridmove_move( VALUE self )
{
  GET_GM_FROM_SELF;
  return ( gm == gridmoveMove( gm )?self:Qnil );
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

VALUE cGridMove;

void Init_GridMove() 
{
  cGridMove = rb_define_class( "GridMove", rb_cObject );
  rb_define_singleton_method( cGridMove, "new", gridmove_new, 1 );
  rb_define_method( cGridMove, "displace", gridmove_displace, 2 );
  rb_define_method( cGridMove, "displacement", gridmove_displacement, 1 );
  rb_define_method( cGridMove, "move", gridmove_move, 0 );
  rb_define_method( cGridMove, "springs", gridmove_springs, 0 );
}
