
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_defs.h"

#include "ref_args.h"
#include "ref_mpi.h"

#include "ref_grid.h"

static void usage(const char *name) {
  printf("usage: \n %s [--help] <command> [<args>]\n", name);
  printf("\n");
  printf("These are common ref commands:\n");
  printf("  boostrap Create initial grid from EGADS file\n");
}
static void bootstrap_help(const char *name) {
  printf("usage: \n %s boostrap project.egads\n", name);
  printf("\n");
}

static REF_STATUS bootstrap(REF_MPI ref_mpi, int argc, char *argv[]) {
  size_t end_of_string;
  char project[1004];
  REF_GRID ref_grid = NULL;
  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref boostrap is no parallel");
  }
  if (argc < 3) goto shutdown;
  end_of_string = strlen(argv[2]);
  if (7 > end_of_string ||
      strncmp(&(argv[2][end_of_string - 6]), ".egads", 6) != 0)
    goto shutdown;
  strncpy(project, argv[2], end_of_string - 7);
  project[end_of_string - 7] = '\0';

  RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
  printf("loading %s.egads\n", project);
  RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[2]), "ld egads");

  RSS(ref_grid_free(ref_grid), "create");

  return REF_SUCCESS;
shutdown:
  bootstrap_help(argv[0]);
  return REF_FAILURE;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT help_pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");
  ref_mpi_stopwatch_start(ref_mpi);

  RXS(ref_args_find(argc, argv, "--help", &help_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY == help_pos) {
    RXS(ref_args_find(argc, argv, "-h", &help_pos), REF_NOT_FOUND,
        "arg search");
  }

  if (1 == argc || 1 == help_pos) {
    usage(argv[0]);
    goto shutdown;
  }

  if (strncmp(argv[1], "b", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(bootstrap(ref_mpi, argc, argv), "bootstrap");
    } else {
      bootstrap_help(argv[0]);
      goto shutdown;
    }
  } else {
    usage(argv[0]);
    goto shutdown;
  }

  ref_mpi_stopwatch_stop(ref_mpi, "done.");
shutdown:
  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
