
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

#include "ref_agents.h"

#include "ref_mpi.h"

int main( int argc, char *argv[] )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_start( argc, argv ), "start" );
  RSS( ref_mpi_create( &ref_mpi ), "make mpi" );

  { /* create */
    REF_AGENTS ref_agents;
    RSS( ref_agents_create( &ref_agents, ref_mpi ), "make agents" );
    RSS( ref_agents_free( ref_agents ), "agents free" );
  }

  { /* remove max id */
    REF_INT id, max;
    REF_AGENTS ref_agents;
    RSS(ref_agents_create(&ref_agents,ref_mpi),"create");

    max = ref_agents->max;
    for ( id = 0; id < max ; id++ )
      RSS(ref_agents_push(ref_agents,id),"add");

    RSS(ref_agents_remove(ref_agents,max-1),"remove last agents");

    RSS(ref_agents_free(ref_agents),"free");
  }

  RSS( ref_mpi_free( ref_mpi ), "mpi free" );
  RSS( ref_mpi_stop( ), "stop" );
  return 0;
}
