
#include "ruby.h"
#include "plan.h"

#define GET_PLAN_FROM_SELF Plan *plan; Data_Get_Struct( self, Plan, plan );

static void plan_free( void *plan )
{
  planFree( plan );
}

VALUE plan_new( VALUE class, VALUE initial_size, VALUE chunk_size )
{
  Plan *plan;
  VALUE obj;
  plan = planCreate( NUM2INT(initial_size), NUM2INT(chunk_size) );
  obj = Data_Wrap_Struct( class, 0, plan_free, plan );
  return obj;
}

VALUE plan_size( VALUE self )
{
  GET_PLAN_FROM_SELF;
  return INT2NUM( planSize(plan) );
}

VALUE plan_chunkSize( VALUE self )
{
  GET_PLAN_FROM_SELF;
  return INT2NUM( planChunkSize(plan) );
}

VALUE cPlan;

void Init_Plan() 
{
  cPlan = rb_define_class( "Plan", rb_cObject );
  rb_define_singleton_method( cPlan, "new", plan_new, 3 );
  rb_define_method( cPlan, "size", plan_size, 0 );
  rb_define_method( cPlan, "chunkSize", plan_chunkSize, 0 );
}
