
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
#include "gridinsert.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_adapt( VALUE self, VALUE minLength, VALUE maxLength)
{
  GET_GRID_FROM_SELF;
  return (gridAdapt( grid, 
		     NUM2DBL(minLength), NUM2DBL(maxLength) )==grid?self:Qnil);
}

VALUE grid_splitEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return INT2NUM(gridSplitEdge( grid, NUM2INT(n0),  NUM2INT(n1) ));
}

VALUE grid_splitEdgeRatio( VALUE self, VALUE n0, VALUE n1, VALUE ratio )
{
  GET_GRID_FROM_SELF;
  return INT2NUM(gridSplitEdgeRatio( grid, NULL, NUM2INT(n0), NUM2INT(n1), 
				     NUM2DBL(ratio) ));
}

VALUE grid_splitFaceAt( VALUE self, VALUE rb_face_nodes, VALUE rb_xyz )
{
  int i;
  int face_nodes[3];
  double xyz[3];
  GET_GRID_FROM_SELF;
  for (i=0;i<3;i++) face_nodes[i] = NUM2INT(rb_ary_entry(rb_face_nodes,i));
  for (i=0;i<3;i++) xyz[i]        = NUM2DBL(rb_ary_entry(rb_xyz,i));
  return INT2NUM(gridSplitFaceAt( grid, face_nodes, xyz ));
}

VALUE grid_splitCellAt( VALUE self, VALUE cell, VALUE rb_xyz )
{
  int i;
  double xyz[3];
  GET_GRID_FROM_SELF;
  for (i=0;i<3;i++) xyz[i] =  NUM2DBL(rb_ary_entry(rb_xyz,i));
  return INT2NUM(gridSplitCellAt( grid, NUM2INT(cell), xyz ));
}

VALUE grid_collapseEdge( VALUE self, VALUE n0, VALUE n1, VALUE ratio )
{
  GET_GRID_FROM_SELF;
  return (gridCollapseEdge( grid, NULL, 
			    NUM2INT(n0), NUM2INT(n1), 
			    NUM2DBL(ratio) )==grid?self:Qnil);
}

VALUE cGridInsert;

void Init_GridInsert() 
{
  cGridInsert = rb_define_module( "GridInsert" );
  rb_define_method( cGridInsert, "adapt", grid_adapt, 2 );
  rb_define_method( cGridInsert, "splitEdge", grid_splitEdge, 2 );
  rb_define_method( cGridInsert, "splitEdgeRatio", grid_splitEdgeRatio, 3 );
  rb_define_method( cGridInsert, "splitFaceAt", grid_splitFaceAt, 2 );
  rb_define_method( cGridInsert, "splitCellAt", grid_splitCellAt, 2 );
  rb_define_method( cGridInsert, "collapseEdge", grid_collapseEdge, 3 );
}
