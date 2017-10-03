
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "FAKEGeom.h"
#include "grid.h"

void inspect_faux(void);

int main( int argc, char *argv[] )
{
  GridBool status;
  int vol = 1;
  int faceId;
  double uv[2], xyz[3], xyznew[3];

  if ( argc < 5 ) return 1;

  faceId = atoi(argv[1]);
  xyz[0] = atof(argv[2]);
  xyz[1] = atof(argv[3]);
  xyz[2] = atof(argv[4]);

  status = CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew );

  printf(" status %d\n", (int)status);
  printf(" face %d\n", faceId);

  printf(" xyz %f %f %f\n", xyz[0], xyz[1], xyz[2]);

  printf(" uv %f %f\n", uv[0], uv[1]);

  printf(" new %f %f %f\n", xyznew[0], xyznew[1], xyznew[2]);
			  
  inspect_faux();

  return 0;
}

