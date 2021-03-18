

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

#ifndef REF_MPI_H
#define REF_MPI_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_MPI_STRUCT REF_MPI_STRUCT;
typedef REF_MPI_STRUCT *REF_MPI;
typedef int REF_TYPE;
#define REF_UNKNOWN_TYPE (0)
#define REF_INT_TYPE (1)
#define REF_LONG_TYPE (2)
#define REF_DBL_TYPE (3)
#define REF_BYTE_TYPE (4)
END_C_DECLORATION

BEGIN_C_DECLORATION
struct REF_MPI_STRUCT {
  REF_INT n;
  REF_INT id;
  void *comm;
  REF_DBL start_time;
  REF_DBL first_time;
  REF_BOOL debug;
};

#define ref_mpi_n(ref_mpi) ((ref_mpi)->n)
#define ref_mpi_rank(ref_mpi) ((ref_mpi)->id)
#define ref_mpi_para(ref_mpi) ((ref_mpi)->n > 1)
#define ref_mpi_once(ref_mpi) (0 == (ref_mpi)->id)

#define each_ref_mpi_part(ref_mpi, part) \
  for ((part) = 0; (part) < ref_mpi_n(ref_mpi); (part)++)

#define each_ref_mpi_worker(ref_mpi, part) \
  for ((part) = 1; (part) < ref_mpi_n(ref_mpi); (part)++)

/* pass void pointer to MPI_Comm */
REF_STATUS ref_mpi_create_from_comm(REF_MPI *ref_mpi, void *comm);

/* assumes MPI_COMM_WORLD */
REF_STATUS ref_mpi_create(REF_MPI *ref_mpi);

REF_STATUS ref_mpi_half_comm(REF_MPI ref_mpi, REF_MPI *split_mpi);
REF_STATUS ref_mpi_front_comm(REF_MPI ref_mpi, REF_MPI *split_mpi, REF_INT n);
REF_STATUS ref_mpi_join_comm(REF_MPI split_mpi);

REF_STATUS ref_mpi_free(REF_MPI ref_mpi);
REF_STATUS ref_mpi_deep_copy(REF_MPI *ref_mpi, REF_MPI original);

REF_STATUS ref_mpi_start(int argc, char *argv[]);
REF_STATUS ref_mpi_stop(void);

REF_STATUS ref_mpi_int_size_type(REF_SIZE size, REF_TYPE *type);

REF_STATUS ref_mpi_elapsed(REF_DBL *seconds);
REF_STATUS ref_mpi_stopwatch_start(REF_MPI ref_mpi);
REF_STATUS ref_mpi_stopwatch_stop(REF_MPI ref_mpi, const char *message);
REF_STATUS ref_mpi_stopwatch_delta(REF_MPI ref_mpi, REF_DBL *delta);

REF_STATUS ref_mpi_bcast(REF_MPI ref_mpi, void *data, REF_INT n, REF_TYPE type);

REF_STATUS ref_mpi_send(REF_MPI ref_mpi, void *data, REF_INT n, REF_TYPE type,
                        REF_INT dest);
REF_STATUS ref_mpi_recv(REF_MPI ref_mpi, void *data, REF_INT n, REF_TYPE type,
                        REF_INT source);

REF_STATUS ref_mpi_alltoall(REF_MPI ref_mpi, void *send, void *recv,
                            REF_TYPE type);
REF_STATUS ref_mpi_alltoallv_native(REF_MPI ref_mpi, void *send,
                                    REF_INT *send_size, void *recv,
                                    REF_INT *recv_size, REF_INT n,
                                    REF_TYPE type);
REF_STATUS ref_mpi_alltoallv(REF_MPI ref_mpi, void *send, REF_INT *send_size,
                             void *recv, REF_INT *recv_size, REF_INT n,
                             REF_TYPE type);

REF_STATUS ref_mpi_all_or(REF_MPI ref_mpi, REF_BOOL *boolean);
REF_STATUS ref_mpi_min(REF_MPI ref_mpi, void *input, void *output,
                       REF_TYPE type);
REF_STATUS ref_mpi_max(REF_MPI ref_mpi, void *input, void *output,
                       REF_TYPE type);
REF_STATUS ref_mpi_sum(REF_MPI ref_mpi, void *input, void *output, REF_INT n,
                       REF_TYPE type);
REF_STATUS ref_mpi_allsum(REF_MPI ref_mpi, void *value, REF_INT n,
                          REF_TYPE type);

REF_STATUS ref_mpi_allgather(REF_MPI ref_mpi, void *scalar, void *array,
                             REF_TYPE type);

REF_STATUS ref_mpi_allgatherv(REF_MPI ref_mpi, void *local_array,
                              REF_INT *counts, void *concatenated_array,
                              REF_TYPE type);

REF_STATUS ref_mpi_allconcat(REF_MPI ref_mpi, REF_INT ldim, REF_INT my_size,
                             void *my_array, REF_INT *total_size,
                             REF_INT **source, void **concatenated,
                             REF_TYPE type);

REF_STATUS ref_mpi_allminwho(REF_MPI ref_mpi, REF_DBL *val, REF_INT *who,
                             REF_INT n);

REF_STATUS ref_mpi_blindsend(REF_MPI ref_mpi, REF_INT *proc, void *send,
                             REF_INT ldim, REF_INT nsend, void **recv,
                             REF_INT *nrecv, REF_TYPE type);
REF_STATUS ref_mpi_balance(REF_MPI ref_mpi, REF_INT ldim, REF_INT nitem,
                           void *items, REF_INT first_rank, REF_INT last_rank,
                           REF_INT *nbalanced, void **balanced, REF_TYPE type);

END_C_DECLORATION

#endif /* REF_MPI_H */
