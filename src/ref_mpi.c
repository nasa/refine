
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"

#ifdef HAVE_MPI

#define ref_type_mpi_type(macro_ref_type, macro_mpi_type) \
  switch (macro_ref_type) {                               \
    case REF_INT_TYPE:                                    \
      (macro_mpi_type) = MPI_INT;                         \
      break;                                              \
    case REF_LONG_TYPE:                                   \
      (macro_mpi_type) = MPI_LONG;                        \
      break;                                              \
    case REF_DBL_TYPE:                                    \
      (macro_mpi_type) = MPI_DOUBLE;                      \
      break;                                              \
    case REF_BYTE_TYPE:                                   \
      (macro_mpi_type) = MPI_UNSIGNED_CHAR;               \
      break;                                              \
    default:                                              \
      (macro_mpi_type) = MPI_DATATYPE_NULL;               \
      RSS(REF_IMPLEMENT, "data type");                    \
  }

#define ref_mpi_comm(ref_mpi) (*((MPI_Comm *)(ref_mpi->comm)))

#define ref_mpi_debugging(ref_mpi) (ref_mpi_once(ref_mpi) && (ref_mpi)->debug)

#define ref_mpi_where_am_i(ref_mpi) \
  if (ref_mpi_debugging(ref_mpi))   \
    printf("%4d: %s: %d: %s:\n", (ref_mpi)->id, __FILE__, __LINE__, __func__);

#endif

REF_STATUS ref_mpi_create_from_comm(REF_MPI *ref_mpi_ptr, void *comm_ptr) {
  REF_MPI ref_mpi;

  ref_malloc(*ref_mpi_ptr, 1, REF_MPI_STRUCT);
  ref_mpi = (*ref_mpi_ptr);

  ref_mpi->id = 0;
  ref_mpi->n = 1;

  ref_mpi->comm = NULL;

  ref_mpi->first_time = (REF_DBL)clock() / ((REF_DBL)CLOCKS_PER_SEC);
  ref_mpi->start_time = ref_mpi->first_time;

  ref_mpi->debug = REF_FALSE;

#ifdef HAVE_MPI
  {
    int running;
    ref_malloc(ref_mpi->comm, 1, MPI_Comm);
    ref_mpi_comm(ref_mpi) = *((MPI_Comm *)(comm_ptr));
    REIS(MPI_SUCCESS, MPI_Initialized(&running), "running?");
    if (running) {
      MPI_Comm_size(ref_mpi_comm(ref_mpi), &(ref_mpi->n));
      MPI_Comm_rank(ref_mpi_comm(ref_mpi), &(ref_mpi->id));
    }
  }
  ref_mpi->first_time = (REF_DBL)MPI_Wtime();
  ref_mpi->start_time = ref_mpi->first_time;
#else
  SUPRESS_UNUSED_COMPILER_WARNING(comm_ptr);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_create(REF_MPI *ref_mpi_ptr) {
#ifdef HAVE_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  RSS(ref_mpi_create_from_comm(ref_mpi_ptr, (void *)(&(comm))),
      "create from MPI_WORLD_COMM");
#else
  RSS(ref_mpi_create_from_comm(ref_mpi_ptr, NULL), "create from NULL comm");
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_half_comm(REF_MPI ref_mpi, REF_MPI *split_mpi_ptr) {
#ifdef HAVE_MPI
  MPI_Comm split_comm;
  int color, rank;
  rank = ref_mpi_rank(ref_mpi);
  color = rank % 2;
  MPI_Comm_split(ref_mpi_comm(ref_mpi), color, rank, &split_comm);
  RSS(ref_mpi_create_from_comm(split_mpi_ptr, &split_comm), "split comm");
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  RSS(ref_mpi_create_from_comm(split_mpi_ptr, NULL), "create from NULL comm");
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_front_comm(REF_MPI ref_mpi, REF_MPI *split_mpi_ptr,
                              REF_INT n) {
#ifdef HAVE_MPI
  MPI_Comm split_comm;
  int color, rank;
  rank = ref_mpi_rank(ref_mpi);
  if (rank < n) {
    color = 0;
  } else {
    color = 1;
  }
  MPI_Comm_split(ref_mpi_comm(ref_mpi), color, rank, &split_comm);
  RSS(ref_mpi_create_from_comm(split_mpi_ptr, &split_comm), "split comm");
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  RSS(ref_mpi_create_from_comm(split_mpi_ptr, NULL), "create from NULL comm");
  SUPRESS_UNUSED_COMPILER_WARNING(n);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_join_comm(REF_MPI split_mpi) {
#ifdef HAVE_MPI
  MPI_Comm_free(&ref_mpi_comm(split_mpi));
#else
  SUPRESS_UNUSED_COMPILER_WARNING(split_mpi);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_mpi_free(REF_MPI ref_mpi) {
  if (NULL == (void *)ref_mpi) return REF_NULL;
  ref_free(ref_mpi->comm);
  ref_free(ref_mpi);
  return REF_SUCCESS;
}

REF_STATUS ref_mpi_deep_copy(REF_MPI *ref_mpi_ptr, REF_MPI original) {
  REF_MPI ref_mpi;

  ref_malloc(*ref_mpi_ptr, 1, REF_MPI_STRUCT);
  ref_mpi = (*ref_mpi_ptr);

  ref_mpi->id = original->id;
  ref_mpi->n = original->n;

  ref_mpi->first_time = original->first_time;
  ref_mpi->start_time = original->start_time;

  ref_mpi->debug = original->debug;

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_start(int argc, char *argv[]) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(argc);
  SUPRESS_UNUSED_COMPILER_WARNING(argv);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_stop() {
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_int_size_type(REF_SIZE size, REF_TYPE *type) {
  *type = REF_UNKNOWN_TYPE;
  switch (size) {
    case sizeof(REF_INT):
      *type = REF_INT_TYPE;
      break;
    case sizeof(REF_LONG):
      *type = REF_LONG_TYPE;
      break;
    default:
      RSB(REF_IMPLEMENT, "data size", {
        printf("size %lu not %lu %lu\n", (unsigned long)size,
               (unsigned long)sizeof(REF_INT), (unsigned long)sizeof(REF_LONG));
      });
  }
  return REF_SUCCESS;
}

REF_STATUS ref_mpi_stopwatch_start(REF_MPI ref_mpi) {
#ifdef HAVE_MPI
  if (ref_mpi_para(ref_mpi)) MPI_Barrier(ref_mpi_comm(ref_mpi));
  ref_mpi->start_time = (REF_DBL)MPI_Wtime();
#else
  ref_mpi->start_time = (REF_DBL)clock() / ((REF_DBL)CLOCKS_PER_SEC);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_stopwatch_stop(REF_MPI ref_mpi, const char *message) {
#ifdef HAVE_MPI
  REF_DBL before_barrier, after_barrier, elapsed;
  REF_DBL first, last;
  before_barrier = (REF_DBL)MPI_Wtime() - ref_mpi->start_time;
  if (ref_mpi_para(ref_mpi)) MPI_Barrier(ref_mpi_comm(ref_mpi));
  after_barrier = (REF_DBL)MPI_Wtime();
  elapsed = after_barrier - ref_mpi->first_time;
  after_barrier = after_barrier - ref_mpi->start_time;
  RSS(ref_mpi_min(ref_mpi, &before_barrier, &first, REF_DBL_TYPE), "min");
  RSS(ref_mpi_min(ref_mpi, &after_barrier, &last, REF_DBL_TYPE), "max");
  if (ref_mpi_once(ref_mpi)) {
    printf("%9.4f: %10.6f (%10.6f) %6.2f%% %s\n", elapsed, last, first,
           (last > 0.0 ? first / last * 100.0 : 100.0), message);
    fflush(stdout);
  }
  RSS(ref_mpi_stopwatch_start(ref_mpi), "restart");
#else
  printf("%9.4f: %10.6f (%10.6f) %6.2f%% %s\n",
         (REF_DBL)clock() / ((REF_DBL)CLOCKS_PER_SEC) - ref_mpi->first_time,
         (REF_DBL)clock() / ((REF_DBL)CLOCKS_PER_SEC) - ref_mpi->start_time,
         (REF_DBL)clock() / ((REF_DBL)CLOCKS_PER_SEC) - ref_mpi->start_time,
         110.0, message);
  fflush(stdout);
  RSS(ref_mpi_stopwatch_start(ref_mpi), "restart");
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_bcast(REF_MPI ref_mpi, void *data, REF_INT n,
                         REF_TYPE type) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_type_mpi_type(type, datatype);

  MPI_Bcast(data, n, datatype, 0, ref_mpi_comm(ref_mpi));
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(data);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_send(REF_MPI ref_mpi, void *data, REF_INT n, REF_TYPE type,
                        REF_INT dest) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;
  REF_INT tag;

  ref_type_mpi_type(type, datatype);

  tag = ref_mpi_n(ref_mpi) * dest + ref_mpi_rank(ref_mpi);

  MPI_Send(data, n, datatype, dest, tag, ref_mpi_comm(ref_mpi));

  return REF_SUCCESS;
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(data);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(dest);
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_mpi_recv(REF_MPI ref_mpi, void *data, REF_INT n, REF_TYPE type,
                        REF_INT source) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;
  REF_INT tag;
  MPI_Status status;

  ref_type_mpi_type(type, datatype);

  tag = ref_mpi_n(ref_mpi) * ref_mpi_rank(ref_mpi) + source;

  MPI_Recv(data, n, datatype, source, tag, ref_mpi_comm(ref_mpi), &status);

  return REF_SUCCESS;
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(data);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(source);
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_mpi_alltoall(REF_MPI ref_mpi, void *send, void *recv,
                            REF_TYPE type) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);

  ref_type_mpi_type(type, datatype);

  MPI_Alltoall(send, 1, datatype, recv, 1, datatype, ref_mpi_comm(ref_mpi));
#else
  switch (type) {
    case REF_INT_TYPE:
      ((REF_INT *)recv)[0] = ((REF_INT *)send)[0];
      break;
    case REF_LONG_TYPE:
      ((REF_LONG *)recv)[0] = ((REF_LONG *)send)[0];
      break;
    case REF_DBL_TYPE:
      ((REF_DBL *)recv)[0] = ((REF_DBL *)send)[0];
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_alltoallv(REF_MPI ref_mpi, void *send, REF_INT *send_size,
                             void *recv, REF_INT *recv_size, REF_INT n,
                             REF_TYPE type) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;
  REF_INT *send_disp;
  REF_INT *recv_disp;
  REF_INT *send_size_n;
  REF_INT *recv_size_n;
  REF_INT part;

  ref_type_mpi_type(type, datatype);

  ref_malloc(send_size_n, ref_mpi_n(ref_mpi), REF_INT);
  each_ref_mpi_part(ref_mpi, part) {
    RAS(0 <= send_size[part], "negative send_size");
    RAS(ref_math_int_multipliable(n, send_size[part]),
        "int overflow send_size_n");
    send_size_n[part] = n * send_size[part];
  }

  ref_malloc(recv_size_n, ref_mpi_n(ref_mpi), REF_INT);
  each_ref_mpi_part(ref_mpi, part) {
    RAS(0 <= recv_size[part], "negative recv_size");
    RAS(ref_math_int_multipliable(n, recv_size[part]),
        "int overflow recv_size_n");
    recv_size_n[part] = n * recv_size[part];
  }

  ref_malloc(send_disp, ref_mpi_n(ref_mpi), REF_INT);
  send_disp[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    RAS(ref_math_int_addable(send_disp[part - 1], send_size_n[part - 1]),
        "int overflow send_disp");
    send_disp[part] = send_disp[part - 1] + send_size_n[part - 1];
  }

  ref_malloc(recv_disp, ref_mpi_n(ref_mpi), REF_INT);
  recv_disp[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    RAS(ref_math_int_addable(recv_disp[part - 1], recv_size_n[part - 1]),
        "int overflow recv_disp");
    recv_disp[part] = recv_disp[part - 1] + recv_size_n[part - 1];
  }

  MPI_Alltoallv(send, send_size_n, send_disp, datatype, recv, recv_size_n,
                recv_disp, datatype, ref_mpi_comm(ref_mpi));

  free(recv_disp);
  free(send_disp);
  free(recv_size_n);
  free(send_size_n);

  return REF_SUCCESS;
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(send);
  SUPRESS_UNUSED_COMPILER_WARNING(send_size);
  SUPRESS_UNUSED_COMPILER_WARNING(recv);
  SUPRESS_UNUSED_COMPILER_WARNING(recv_size);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_mpi_min(REF_MPI ref_mpi, void *input, void *output,
                       REF_TYPE type) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  if (!ref_mpi_para(ref_mpi)) {
    switch (type) {
      case REF_INT_TYPE:
        *(REF_INT *)output = *(REF_INT *)input;
        break;
      case REF_DBL_TYPE:
        *(REF_DBL *)output = *(REF_DBL *)input;
        break;
      default:
        RSS(REF_IMPLEMENT, "data type");
    }
    return REF_SUCCESS;
  }
  ref_type_mpi_type(type, datatype);

  ref_mpi_where_am_i(ref_mpi);
  MPI_Reduce(input, output, 1, datatype, MPI_MIN, 0, ref_mpi_comm(ref_mpi));

#else
  switch (type) {
    case REF_INT_TYPE:
      *(REF_INT *)output = *(REF_INT *)input;
      break;
    case REF_DBL_TYPE:
      *(REF_DBL *)output = *(REF_DBL *)input;
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_all_or(REF_MPI ref_mpi, REF_BOOL *boolean) {
#ifdef HAVE_MPI
  REF_BOOL output;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_mpi_where_am_i(ref_mpi);
  MPI_Allreduce(boolean, &output, 1, MPI_INT, MPI_SUM, ref_mpi_comm(ref_mpi));
  *boolean = MIN(output, 1);

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(boolean);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_max(REF_MPI ref_mpi, void *input, void *output,
                       REF_TYPE type) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  if (!ref_mpi_para(ref_mpi)) {
    switch (type) {
      case REF_INT_TYPE:
        *(REF_INT *)output = *(REF_INT *)input;
        break;
      case REF_DBL_TYPE:
        *(REF_DBL *)output = *(REF_DBL *)input;
        break;
      default:
        RSS(REF_IMPLEMENT, "data type");
    }
    return REF_SUCCESS;
  }
  ref_type_mpi_type(type, datatype);

  ref_mpi_where_am_i(ref_mpi);
  MPI_Reduce(input, output, 1, datatype, MPI_MAX, 0, ref_mpi_comm(ref_mpi));

#else
  switch (type) {
    case REF_INT_TYPE:
      *(REF_INT *)output = *(REF_INT *)input;
      break;
    case REF_DBL_TYPE:
      *(REF_DBL *)output = *(REF_DBL *)input;
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_sum(REF_MPI ref_mpi, void *input, void *output, REF_INT n,
                       REF_TYPE type) {
  REF_INT i;
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  if (!ref_mpi_para(ref_mpi)) {
    switch (type) {
      case REF_INT_TYPE:
        for (i = 0; i < n; i++) ((REF_INT *)output)[i] = ((REF_INT *)input)[i];
        break;
      case REF_LONG_TYPE:
        for (i = 0; i < n; i++)
          ((REF_LONG *)output)[i] = ((REF_LONG *)input)[i];
        break;
      case REF_DBL_TYPE:
        for (i = 0; i < n; i++) ((REF_DBL *)output)[i] = ((REF_DBL *)input)[i];
        break;
      default:
        RSS(REF_IMPLEMENT, "data type");
    }
  } else {
    ref_type_mpi_type(type, datatype);
    ref_mpi_where_am_i(ref_mpi);
    MPI_Reduce(input, output, n, datatype, MPI_SUM, 0, ref_mpi_comm(ref_mpi));
  }

#else
  switch (type) {
    case REF_INT_TYPE:
      for (i = 0; i < n; i++) ((REF_INT *)output)[i] = ((REF_INT *)input)[i];
      break;
    case REF_LONG_TYPE:
      for (i = 0; i < n; i++) ((REF_LONG *)output)[i] = ((REF_LONG *)input)[i];
      break;
    case REF_DBL_TYPE:
      for (i = 0; i < n; i++) ((REF_DBL *)output)[i] = ((REF_DBL *)input)[i];
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_allsum(REF_MPI ref_mpi, void *value, REF_INT n,
                          REF_TYPE type) {
  REF_INT i;
  void *temp = NULL;
  switch (type) {
    case REF_INT_TYPE:
      ref_malloc(temp, n, REF_INT);
      for (i = 0; i < n; i++) ((REF_INT *)temp)[i] = ((REF_INT *)value)[i];
      break;
    case REF_LONG_TYPE:
      ref_malloc(temp, n, REF_LONG);
      for (i = 0; i < n; i++) ((REF_LONG *)temp)[i] = ((REF_LONG *)value)[i];
      break;
    case REF_DBL_TYPE:
      ref_malloc(temp, n, REF_DBL);
      for (i = 0; i < n; i++) ((REF_DBL *)temp)[i] = ((REF_DBL *)value)[i];
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  RSS(ref_mpi_sum(ref_mpi, temp, value, n, type), "sum");
  ref_free(temp);
  RSS(ref_mpi_bcast(ref_mpi, value, n, type), "bcast")
  return REF_SUCCESS;
}

REF_STATUS ref_mpi_allgather(REF_MPI ref_mpi, void *scalar, void *array,
                             REF_TYPE type) {
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  if (!ref_mpi_para(ref_mpi)) {
    switch (type) {
      case REF_INT_TYPE:
        *(REF_INT *)array = *(REF_INT *)scalar;
        break;
      case REF_LONG_TYPE:
        *(REF_LONG *)array = *(REF_LONG *)scalar;
        break;
      case REF_DBL_TYPE:
        *(REF_DBL *)array = *(REF_DBL *)scalar;
        break;
      default:
        RSS(REF_IMPLEMENT, "data type");
    }
    return REF_SUCCESS;
  }

  ref_type_mpi_type(type, datatype);

  MPI_Allgather(scalar, 1, datatype, array, 1, datatype, ref_mpi_comm(ref_mpi));

#else
  switch (type) {
    case REF_INT_TYPE:
      *(REF_INT *)array = *(REF_INT *)scalar;
      break;
    case REF_LONG_TYPE:
      *(REF_LONG *)array = *(REF_LONG *)scalar;
      break;
    case REF_DBL_TYPE:
      *(REF_DBL *)array = *(REF_DBL *)scalar;
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_allgatherv(REF_MPI ref_mpi, void *local_array,
                              REF_INT *counts, void *concatenated_array,
                              REF_TYPE type) {
#ifdef HAVE_MPI
  REF_INT proc;
  REF_INT *displs;
  MPI_Datatype datatype;
  REF_INT i;
  ref_type_mpi_type(type, datatype);

  if (!ref_mpi_para(ref_mpi)) {
    switch (type) {
      case REF_INT_TYPE:
        for (i = 0; i < counts[0]; i++)
          ((REF_INT *)concatenated_array)[i] = ((REF_INT *)local_array)[i];
        break;
      case REF_LONG_TYPE:
        for (i = 0; i < counts[0]; i++)
          ((REF_LONG *)concatenated_array)[i] = ((REF_LONG *)local_array)[i];
        break;
      case REF_DBL_TYPE:
        for (i = 0; i < counts[0]; i++)
          ((REF_DBL *)concatenated_array)[i] = ((REF_DBL *)local_array)[i];
        break;
      default:
        RSS(REF_IMPLEMENT, "data type");
    }
    return REF_SUCCESS;
  }

  ref_malloc(displs, ref_mpi_n(ref_mpi), REF_INT);

  displs[0] = 0;
  each_ref_mpi_worker(ref_mpi, proc) displs[proc] =
      displs[proc - 1] + counts[proc - 1];

  MPI_Allgatherv(local_array, counts[ref_mpi_rank(ref_mpi)], datatype,
                 concatenated_array, counts, displs, datatype,
                 ref_mpi_comm(ref_mpi));

  ref_free(displs);

#else
  REF_INT i;
  switch (type) {
    case REF_INT_TYPE:
      for (i = 0; i < counts[0]; i++)
        ((REF_INT *)concatenated_array)[i] = ((REF_INT *)local_array)[i];
      break;
    case REF_LONG_TYPE:
      for (i = 0; i < counts[0]; i++)
        ((REF_LONG *)concatenated_array)[i] = ((REF_LONG *)local_array)[i];
      break;
    case REF_DBL_TYPE:
      for (i = 0; i < counts[0]; i++)
        ((REF_DBL *)concatenated_array)[i] = ((REF_DBL *)local_array)[i];
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }
  SUPRESS_UNUSED_COMPILER_WARNING(ref_mpi);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_allconcat(REF_MPI ref_mpi, REF_INT ldim, REF_INT my_size,
                             void *my_array, REF_INT *total_size,
                             REF_INT **source, void **concatenated,
                             REF_TYPE type) {
  REF_INT proc, i, tot;
  REF_INT *counts;

  ref_malloc(counts, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, (void *)(&my_size), (void *)counts,
                        REF_INT_TYPE),
      "gather size");

  *total_size = 0;
  each_ref_mpi_part(ref_mpi, proc)(*total_size) += counts[proc];
  ref_malloc(*source, *total_size, REF_INT);

  tot = 0;
  each_ref_mpi_part(ref_mpi, proc) for (i = 0; i < counts[proc]; i++) {
    (*source)[tot] = proc;
    tot++;
  }
  REIS(*total_size, tot, "count mismatch");

  switch (type) {
    case REF_INT_TYPE:
      ref_malloc(*((REF_INT **)concatenated), ldim * (*total_size), REF_INT);
      break;
    case REF_DBL_TYPE:
      ref_malloc(*((REF_DBL **)concatenated), ldim * (*total_size), REF_DBL);
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }

  /* pad ldim */
  each_ref_mpi_part(ref_mpi, proc) counts[proc] *= ldim;

  RSS(ref_mpi_allgatherv(ref_mpi, my_array, counts, *concatenated, type),
      "gather values");

  ref_free(counts);

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_allminwho(REF_MPI ref_mpi, REF_DBL *val, REF_INT *who,
                             REF_INT n) {
  REF_INT i;
  if (!ref_mpi_para(ref_mpi)) {
    for (i = 0; i < n; i++) who[i] = ref_mpi_rank(ref_mpi);
    return REF_SUCCESS;
  }
#ifdef HAVE_MPI
  {
    typedef struct REF_MPI_DBLWHO REF_MPI_DBLWHO;
    struct REF_MPI_DBLWHO {
      double val;
      int rank;
    };
    REF_MPI_DBLWHO *in, *out;
    ref_malloc(in, n, REF_MPI_DBLWHO);
    ref_malloc(out, n, REF_MPI_DBLWHO);

    for (i = 0; i < n; i++) {
      in[i].val = val[i];
      in[i].rank = ref_mpi_rank(ref_mpi);
    }
    ref_mpi_where_am_i(ref_mpi);
    if (ref_mpi_debugging(ref_mpi)) printf("MPI_MINLOC %d\n", n);
    MPI_Allreduce(in, out, n, MPI_DOUBLE_INT, MPI_MINLOC,
                  ref_mpi_comm(ref_mpi));
    for (i = 0; i < n; i++) {
      val[i] = out[i].val;
      who[i] = out[i].rank;
    }

    ref_free(out);
    ref_free(in);
  }
#else
  SUPRESS_UNUSED_COMPILER_WARNING(val);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_mpi_blindsend(REF_MPI ref_mpi, REF_INT *proc, void *send,
                             REF_INT ldim, REF_INT nsend, void **recv,
                             REF_INT *nrecv, REF_TYPE type) {
  REF_INT i, l, part;
  REF_INT *a_size, *b_size;
  REF_INT *a_next;
  REF_INT a_total, b_total;

  void *a_data;

  a_data = NULL;
  *recv = NULL;

  if (!ref_mpi_para(ref_mpi)) {
    *nrecv = nsend;
    switch (type) {
      case REF_INT_TYPE:
        ref_malloc(*((REF_INT **)recv), ldim * nsend, REF_INT);
        for (i = 0; i < ldim * nsend; i++) {
          (*((REF_INT **)recv))[i] = ((REF_INT *)send)[i];
        }
        break;
      case REF_LONG_TYPE:
        ref_malloc(*((REF_LONG **)recv), ldim * nsend, REF_LONG);
        for (i = 0; i < ldim * nsend; i++) {
          (*((REF_LONG **)recv))[i] = ((REF_LONG *)send)[i];
        }
        break;
      case REF_DBL_TYPE:
        ref_malloc(*((REF_DBL **)recv), ldim * nsend, REF_DBL);
        for (i = 0; i < ldim * nsend; i++) {
          (*((REF_DBL **)recv))[i] = ((REF_DBL *)send)[i];
        }
        break;
      default:
        RSS(REF_IMPLEMENT, "data type");
    }
    return REF_SUCCESS;
  }

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  for (i = 0; i < nsend; i++) a_size[proc[i]]++;

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  a_data = NULL;
  *recv = NULL;
  switch (type) {
    case REF_INT_TYPE:
      ref_malloc(a_data, ldim * a_total, REF_INT);
      ref_malloc(*((REF_INT **)recv), ldim * b_total, REF_INT);
      for (i = 0; i < nsend; i++) {
        for (l = 0; l < ldim; l++)
          ((REF_INT *)a_data)[l + ldim * a_next[proc[i]]] =
              ((REF_INT *)send)[l + ldim * i];
        a_next[proc[i]]++;
      }
      break;
    case REF_LONG_TYPE:
      ref_malloc(a_data, ldim * a_total, REF_LONG);
      ref_malloc(*((REF_LONG **)recv), ldim * b_total, REF_LONG);
      for (i = 0; i < nsend; i++) {
        for (l = 0; l < ldim; l++)
          ((REF_LONG *)a_data)[l + ldim * a_next[proc[i]]] =
              ((REF_LONG *)send)[l + ldim * i];
        a_next[proc[i]]++;
      }
      break;
    case REF_DBL_TYPE:
      ref_malloc(a_data, ldim * a_total, REF_DBL);
      ref_malloc(*((REF_DBL **)recv), ldim * b_total, REF_DBL);
      for (i = 0; i < nsend; i++) {
        for (l = 0; l < ldim; l++)
          ((REF_DBL *)a_data)[l + ldim * a_next[proc[i]]] =
              ((REF_DBL *)send)[l + ldim * i];
        a_next[proc[i]]++;
      }
      break;
    default:
      RSS(REF_IMPLEMENT, "data type");
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_data, a_size, *recv, b_size, ldim, type),
      "alltoallv global");

  *nrecv = b_total;

  ref_free(a_next);
  ref_free(a_data);
  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

static REF_INT find_destination(REF_INT n, REF_INT *shares, REF_INT gid) {
  REF_INT part;
  for (part = 0; part < n; part++) {
    if (gid < shares[part]) return part;
    gid -= shares[part];
  }
  return n - 1;
}

REF_STATUS ref_mpi_balance(REF_MPI ref_mpi, REF_INT ldim, REF_INT nitem,
                           void *items, REF_INT first_rank, REF_INT last_rank,
                           REF_INT *nbalanced, void **balanced, REF_TYPE type) {
  REF_INT *shares, total, share, remainder;
  REF_INT *haves;
  REF_INT *destination, offset, part, i;
  REF_INT active;

  active = last_rank - first_rank + 1;

  total = nitem;
  RSS(ref_mpi_allsum(ref_mpi, &total, 1, REF_INT_TYPE), "total items");
  if (first_rank <= ref_mpi_rank(ref_mpi) &&
      ref_mpi_rank(ref_mpi) <= last_rank) {
    share = total / active;
    remainder = total - share * active;
  } else {
    share = 0;
    remainder = 0;
  }
  if (MAX(0, ref_mpi_rank(ref_mpi) - first_rank) < remainder) {
    share++;
  }
  *nbalanced = share;

  ref_malloc(haves, ref_mpi_n(ref_mpi), REF_INT);
  ref_malloc(shares, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, &share, shares, REF_INT_TYPE), "all share");
  RSS(ref_mpi_allgather(ref_mpi, &nitem, haves, REF_INT_TYPE), "all share");

  ref_malloc_init(destination, nitem, REF_INT, 0);
  offset = 0;
  for (part = 0; part < ref_mpi_rank(ref_mpi); part++) {
    offset += haves[part];
  }
  for (i = 0; i < nitem; i++) {
    destination[i] = find_destination(ref_mpi_n(ref_mpi), shares, offset + i);
  }

  RSS(ref_mpi_blindsend(ref_mpi, destination, items, ldim, nitem, balanced,
                        nbalanced, type),
      "blind send node");
  REIB(share, *nbalanced, "share mismatch", {
    printf(
        "rank %d first %d last %d active %d total %d rem %d share %d nitem "
        "%d\n",
        ref_mpi_rank(ref_mpi), first_rank, last_rank, active, total, remainder,
        share, nitem);
    for (i = 0; i < ref_mpi_n(ref_mpi); i++) printf(" %d", shares[i]);
    printf("\n");
  });

  ref_free(destination);
  ref_free(shares);
  ref_free(haves);

  return REF_SUCCESS;
}
