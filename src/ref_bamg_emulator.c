
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
#include <time.h>

#include "ref_adapt.h"
#include "ref_args.h"
#include "ref_export.h"
#include "ref_grid.h"
#include "ref_histogram.h"
#include "ref_import.h"
#include "ref_malloc.h"
#include "ref_metric.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_part.h"
#include "ref_validation.h"

int main(int argc, char *argv[]) {
  char *input_filename = NULL;
  char *output_filename = NULL;
  char *metric_filename = NULL;
  REF_BOOL noop;
  REF_INT location;
  REF_INT i, passes;

  REF_MPI ref_mpi;
  REF_GRID ref_grid = NULL;
  REF_GRID background_grid = NULL;
  char timestamp[512];
  char command[1024];

  time_t rawtime = time(0);
  struct tm *now = localtime(&rawtime);
  RAS(0 < strftime(timestamp, 512, "%Y-%m-%d_%H:%M:%S", now),
      "unable to format time stamp");

  RSS(ref_args_find(argc, argv, "-b", &location), "-b argument missing");
  printf(" %s ", argv[location]);
  RAS(location < argc - 1, "-b missing");
  input_filename = argv[1 + location];
  printf("'%s'\n", input_filename);

  RXS(ref_args_find(argc, argv, "-noop", &location), REF_NOT_FOUND, "noop");
  noop = (location != REF_EMPTY);

  if (!noop) {
    RSS(ref_args_find(argc, argv, "-M", &location), "-M argument missing");
    printf(" %s ", argv[location]);
    RAS(location < argc - 1, "-M missing");
    metric_filename = argv[1 + location];
    printf("'%s'\n", metric_filename);
  }

  RSS(ref_args_find(argc, argv, "-o", &location), "-o argument missing");
  printf(" %s ", argv[location]);
  RAS(location < argc - 1, "-o missing");
  output_filename = argv[1 + location];
  printf("'%s'\n", output_filename);

  printf("start up at %s\n", timestamp);
  RAS(0 < sprintf(command, "cp %s ref-bu-%s-in.msh", input_filename, timestamp),
      "in");
  printf("%s\n", command);
  REIS(0, system(command), "cp command failed");

  RSS(ref_mpi_create(&ref_mpi), "create");
  RSS(ref_import_by_extension(&ref_grid, ref_mpi, input_filename), "in");
  RSS(ref_import_by_extension(&background_grid, ref_mpi, input_filename), "in");
  ref_grid_inspect(ref_grid);

  if (noop) {
    REF_DBL *metric;
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_imply_from(metric, ref_grid), "imply");
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "to");
    ref_free(metric);
  } else {
    RSS(ref_part_bamg_metric(ref_grid, metric_filename), "metric");
    if (NULL != background_grid)
      RSS(ref_part_bamg_metric(background_grid, metric_filename), "metric");
    RAS(0 < sprintf(command, "cp %s ref-bu-%s.metric", metric_filename,
                    timestamp),
        "in");
    printf("%s\n", command);
    REIS(0, system(command), "cp command failed");
  }

  RSS(ref_validation_cell_volume(ref_grid), "vol");
  RSS(ref_histogram_quality(ref_grid), "gram");
  RSS(ref_histogram_ratio(ref_grid), "gram");

  passes = 20;
  for (i = 0; i < passes; i++) {
    REF_BOOL all_done;
    printf(" pass %d of %d\n", i, passes);
    RSS(ref_adapt_pass(ref_grid, &all_done), "pass");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "pass");
    if (NULL != background_grid) {
      RSS(ref_metric_interpolate_twod(ref_grid, background_grid), "interp");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "interp");
    }
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    RSS(ref_migrate_to_balance(ref_grid), "balance");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "balance");
  }

  RSS(ref_export_by_extension(ref_grid, output_filename), "out");

  RAS(0 < sprintf(command, "cp %s ref-bu-%s-out.msh", output_filename,
                  timestamp),
      "in");
  printf("%s\n", command);
  REIS(0, system(command), "cp command failed");

  RSS(ref_grid_free(background_grid), "free");
  RSS(ref_grid_free(ref_grid), "free");
  RSS(ref_mpi_free(ref_mpi), "free");
  return 0;
}
