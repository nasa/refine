
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

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

VALUE plan_max_size( VALUE self )
{
  GET_PLAN_FROM_SELF;
  return INT2NUM( planMaxSize(plan) );
}

VALUE plan_chunk_size( VALUE self )
{
  GET_PLAN_FROM_SELF;
  return INT2NUM( planChunkSize(plan) );
}

VALUE plan_add_item_with_priority( VALUE self, VALUE item, VALUE priority )
{
  GET_PLAN_FROM_SELF;
  return (planAddItemWithPriority(plan,NUM2INT(item),NUM2DBL(priority))==NULL?Qnil:self);
}

VALUE plan_derive_rankings_from_priorities( VALUE self )
{
  GET_PLAN_FROM_SELF;
  return (planDeriveRankingsFromPriorities( plan )==NULL?Qnil:self);  
}

VALUE plan_item_with_this_ranking( VALUE self, VALUE ranking )
{
  GET_PLAN_FROM_SELF;
  return INT2NUM(planItemWithThisRanking(plan,NUM2INT(ranking)));
}

VALUE plan_priority_with_this_ranking( VALUE self, VALUE ranking )
{
  GET_PLAN_FROM_SELF;
  return rb_float_new(planPriorityWithThisRanking(plan,NUM2INT(ranking)));
}

VALUE cPlan;

void Init_Plan() 
{
  cPlan = rb_define_class( "Plan", rb_cObject );
  rb_define_singleton_method( cPlan, "new", plan_new, 2 );
  rb_define_method( cPlan, "size", plan_size, 0 );
  rb_define_method( cPlan, "max_size", plan_max_size, 0 );
  rb_define_method( cPlan, "chunk_size", plan_chunk_size, 0 );

  rb_define_method( cPlan, "add_item_with_priority", 
		    plan_add_item_with_priority, 2 );
  rb_define_method( cPlan, "derive_rankings_from_priorities", 
		    plan_derive_rankings_from_priorities, 0 );
  rb_define_method( cPlan, "item_with_this_ranking", 
		    plan_item_with_this_ranking, 1 );
  rb_define_method( cPlan, "priority_with_this_ranking", 
		    plan_priority_with_this_ranking, 1 );
}
