
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

#include "ref_mpi.h"

#include "ref_malloc.h"

int main( int argc, char *argv[] )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_start( argc, argv ), "start" );
  RSS( ref_mpi_create( &ref_mpi ), "make mpi" );

  if ( ref_mpi_para(ref_mpi) && ref_mpi_once(ref_mpi) )
    printf("%s number of processors %d \n",argv[0],ref_mpi_m(ref_mpi));
  
  if ( !ref_mpi_para(ref_mpi) )
    { /* no mpi or mpi with one proc */
      REIS( 1, ref_mpi_m(ref_mpi), "n" );
      REIS( 0, ref_mpi_rank(ref_mpi), "rank" );
      RAS( ref_mpi_once(ref_mpi), "master" );
    }

  /* bcast */
  {
    REF_INT bc;

    bc = REF_EMPTY;
    if ( ref_mpi_once(ref_mpi) ) bc = 5;
    RSS( ref_mpi_bcast( ref_mpi, &bc, 1, REF_INT_TYPE ), "bcast" );
    REIS( 5, bc, "bc wrong" );
  }

  /* alltoall */
  {
    REF_INT part;
    REF_INT *a_size, *b_size;

    ref_malloc_init( a_size, ref_mpi_m(ref_mpi), REF_INT, REF_EMPTY );
    ref_malloc_init( b_size, ref_mpi_m(ref_mpi), REF_INT, REF_EMPTY );
      
    each_ref_mpi_part(ref_mpi, part)
      a_size[part] = part;

    RSS( ref_mpi_alltoall( ref_mpi,
			   a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

    each_ref_mpi_part(ref_mpi, part)
      REIS(part, a_size[part], "a_size changed" );
    each_ref_mpi_part(ref_mpi, part)
      REIS(ref_mpi_rank(ref_mpi), b_size[part], "b_size wrong" );

    ref_free( b_size );
    ref_free( a_size );
  }

  /* allconcat */
  {
    REF_INT total, part;
    REF_INT *separate, *together;

    ref_malloc( separate, 2, REF_INT );
      
    separate[0] = ref_mpi_rank(ref_mpi);
    separate[1] = REF_EMPTY;

    RSS( ref_mpi_allconcat( ref_mpi,
			    2, (void *)separate,
			    &total, (void **)&together, REF_INT_TYPE ), "allconcat");

    REIS( 2*ref_mpi_m(ref_mpi), total, "expected size" );
    each_ref_mpi_part( ref_mpi, part )
      {
	REIS( part, together[0+2*part], "const" );
	REIS( REF_EMPTY, together[1+2*part], "const" );
      }

    ref_free( together );
    ref_free( separate );
  }

  /* allminwho */
  {
    REF_INT part;
    REF_DBL *val;
    REF_INT *who;

    ref_malloc( val, ref_mpi_m(ref_mpi), REF_DBL );
    ref_malloc( who, ref_mpi_m(ref_mpi), REF_INT );
      
    each_ref_mpi_part( ref_mpi, part )
      val[part] = 1.0;
    val[ref_mpi_rank(ref_mpi)] = 0.0;

    RSS( ref_mpi_allminwho( ref_mpi, val, who, ref_mpi_m(ref_mpi) ), "minloc");

    each_ref_mpi_part( ref_mpi, part )
      {
	REIS( part, who[part], "who" );
	RWDS( 0.0, val[part], -1.0, "min" );
      }

    ref_free( who );
    ref_free( val );
  }

  /* blindsend int with 0,0 */
  {
    REF_INT i, part, size, n;
    REF_INT *send;
    REF_INT *proc;
    REF_INT *recv;

    size = 0;
    each_ref_mpi_part( ref_mpi, part )
      size += MIN(part,10);

    ref_malloc( proc, size, REF_INT );
    ref_malloc( send, size, REF_INT );
    
    size = 0;
    each_ref_mpi_part( ref_mpi, part )
      for ( i=0;i<MIN(part,10);i++ )
	{
	  proc[size] = part;
	  send[size] = ref_mpi_rank(ref_mpi);
	  size++;
	}

    RSS( ref_mpi_blindsend( ref_mpi, proc, (void*)send, size,
			    (void **)&recv, &n, REF_INT_TYPE ),
	 "blindsend");

    size = 0;
    each_ref_mpi_part( ref_mpi, part )
      for ( i=0;i<MIN(ref_mpi_rank(ref_mpi),10);i++ )
	{
	  REIS(part,recv[size],"recv mismatch");
	  size++;
	}
    REIS( size, n, "size mismatch" );
    
    ref_free( recv );
    ref_free( send );
    ref_free( proc );
  }

  /* blindsend dbl with 0,1 */
  {
    REF_INT i, part, size, n;
    REF_DBL *send;
    REF_INT *proc;
    REF_DBL *recv;

    size = 0;
    each_ref_mpi_part( ref_mpi, part )
      size += MAX(1,MIN(part,10));

    ref_malloc( proc, size, REF_INT );
    ref_malloc( send, size, REF_DBL );
    
    size = 0;
    each_ref_mpi_part( ref_mpi, part )
      for ( i=0;i<MAX(1,MIN(part,10));i++ )
	{
	  proc[size] = part;
	  send[size] = (REF_DBL)ref_mpi_rank(ref_mpi);
	  size++;
	}

    RSS( ref_mpi_blindsend( ref_mpi, proc, (void*)send, size,
			    (void **)&recv, &n, REF_DBL_TYPE ),
	 "blindsend");

    size = 0;
    each_ref_mpi_part( ref_mpi, part )
      for ( i=0;i<MAX(1,MIN(ref_mpi_rank(ref_mpi),10));i++ )
	{
	  RWDS((REF_DBL)part,recv[size],-1,"recv mismatch");
	  size++;
	}
    REIS( size, n, "size mismatch" );
    
    ref_free( recv );
    ref_free( send );
    ref_free( proc );
  }

  RSS( ref_mpi_free( ref_mpi ), "mpi free" );
  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
