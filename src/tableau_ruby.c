
#include "ruby.h"
#include "tableau.h"

#define GET_TABLEAU_FROM_SELF Tableau *tableau; Data_Get_Struct( self, Tableau, tableau );

static void tableau_free( void *tableau )
{
  tableauFree( tableau );
}

VALUE tableau_new( VALUE class, VALUE rb_constraints, VALUE rb_dimension )
{
  Tableau *tableau;
  VALUE obj;
  tableau = tableauCreate( NUM2INT(rb_constraints), NUM2INT(rb_dimension) );
  obj = Data_Wrap_Struct( class, 0, tableau_free, tableau );
  return obj;
}

VALUE tableau_constraints( VALUE self )
{
  GET_TABLEAU_FROM_SELF;
  return INT2NUM( tableauConstraints(tableau) );
}

VALUE tableau_dimension( VALUE self )
{
  GET_TABLEAU_FROM_SELF;
  return INT2NUM( tableauDimension(tableau) );
}

VALUE cTableau;

void Init_Tableau() 
{
  cTableau = rb_define_class( "Tableau", rb_cObject );
  rb_define_singleton_method( cTableau, "new", tableau_new, 2 );

  rb_define_method( cTableau, "constraints", tableau_constraints, 0 );
  rb_define_method( cTableau, "dimension", tableau_dimension, 0 );
}
