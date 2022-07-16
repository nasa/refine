
/* Copyright 2006, 2014, 2021 United States Government as represented
 * by the Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine version 3 unstructured grid adaptation platform is
 * licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include "ref_oct.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_mpi.h"

int main(int argc, char *argv[]) {
REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  { /* create */
    REF_OCT ref_oct;
    RSS(ref_oct_create(&ref_oct), "make oct");
    RSS(ref_oct_free(ref_oct), "search oct");
  }

RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
