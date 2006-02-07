
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

VALUE tableau_constraintMatrix( VALUE self, VALUE rb_data )
{
  int i, ndata;
  double *data;
  VALUE rb_result;
  GET_TABLEAU_FROM_SELF;
  ndata = RARRAY(rb_data)->len;
  data = (double *)malloc(ndata*sizeof(double));
  for ( i=0 ; i < ndata ; i++ ) data[i] = NUM2DBL( rb_ary_entry( rb_data, i) );
  rb_result = (tableau == tableauConstraintMatrix(tableau, data)?self:Qnil);
  free(data);
  return rb_result;
}

VALUE tableau_cost( VALUE self, VALUE rb_data )
{
  int i, ndata;
  double *data;
  VALUE rb_result;
  GET_TABLEAU_FROM_SELF;
  ndata = RARRAY(rb_data)->len;
  data = (double *)malloc(ndata*sizeof(double));
  for ( i=0 ; i < ndata ; i++ ) data[i] = NUM2DBL( rb_ary_entry( rb_data, i) );
  rb_result = (tableau == tableauCost(tableau, data)?self:Qnil);
  free(data);
  return rb_result;
}

VALUE tableau_constraint( VALUE self, VALUE rb_data )
{
  int i, ndata;
  double *data;
  VALUE rb_result;
  GET_TABLEAU_FROM_SELF;
  ndata = RARRAY(rb_data)->len;
  data = (double *)malloc(ndata*sizeof(double));
  for ( i=0 ; i < ndata ; i++ ) data[i] = NUM2DBL( rb_ary_entry( rb_data, i) );
  rb_result = (tableau == tableauConstraint(tableau, data)?self:Qnil);
  free(data);
  return rb_result;
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

VALUE tableau_basis( VALUE self )
{
  VALUE rb_basis;
  int i, nbasis;
  int *basis;
  GET_TABLEAU_FROM_SELF;
  nbasis = tableauConstraints(tableau);
  basis = (int *)malloc(nbasis*sizeof(int));
  if (tableau != tableauBasis(tableau, basis) ) {
    free(basis);
    return Qnil;
  }
  rb_basis = rb_ary_new2(nbasis);
  for ( i=0 ; i < nbasis ; i++ ) 
    rb_ary_store( rb_basis, i, INT2NUM(basis[i]) );
  return rb_basis;
}

VALUE tableau_solve( VALUE self )
{
  GET_TABLEAU_FROM_SELF;
  return ( tableau == tableauSolve(tableau) ? self : Qnil );
}

VALUE cTableau;

void Init_Tableau() 
{
  cTableau = rb_define_class( "Tableau", rb_cObject );
  rb_define_singleton_method( cTableau, "new", tableau_new, 2 );

  rb_define_method( cTableau, "constraintMatrix", tableau_constraintMatrix, 1);
  rb_define_method( cTableau, "constraint", tableau_constraint, 1);
  rb_define_method( cTableau, "cost", tableau_cost, 1);

  rb_define_method( cTableau, "constraints", tableau_constraints, 0 );
  rb_define_method( cTableau, "dimension", tableau_dimension, 0 );

  rb_define_method( cTableau, "basis", tableau_basis, 0 );

  rb_define_method( cTableau, "solve", tableau_solve, 0 );
}
