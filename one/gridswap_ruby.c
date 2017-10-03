
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
#include "gridswap.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_swapFace( VALUE self, VALUE n0, VALUE n1, VALUE n2)
{
  GET_GRID_FROM_SELF;
  return (gridSwapFace( grid, NULL, NUM2INT(n0),  
			NUM2INT(n1), NUM2INT(n2) )==grid?self:Qnil);
}

VALUE grid_swapEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return (gridSwapEdge( grid, NULL, NUM2INT(n0), NUM2INT(n1) )==grid?self:Qnil);
}

VALUE grid_swap( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSwap( grid, -1.0 )==grid?self:Qnil);
}

VALUE grid_removeTwoFaceCell( VALUE self, VALUE cell )
{
  GET_GRID_FROM_SELF;
  return (gridRemoveTwoFaceCell( grid, NULL, NUM2INT(cell) )==grid?self:Qnil);
}

VALUE grid_removeThreeFaceCell( VALUE self, VALUE cell )
{
  GET_GRID_FROM_SELF;
  return (gridRemoveThreeFaceCell( grid, NULL, NUM2INT(cell) )==grid?self:Qnil);
}

VALUE cGridSwap;

void Init_GridSwap() 
{
  cGridSwap = rb_define_module( "GridSwap" );
  rb_define_method( cGridSwap, "swapFace", grid_swapFace, 3 );
  rb_define_method( cGridSwap, "swapEdge", grid_swapEdge, 2 );
  rb_define_method( cGridSwap, "swap", grid_swap, 0 );
  rb_define_method( cGridSwap, "removeTwoFaceCell", grid_removeTwoFaceCell, 1 );
  rb_define_method( cGridSwap, "removeThreeFaceCell", grid_removeThreeFaceCell, 1 );
}
