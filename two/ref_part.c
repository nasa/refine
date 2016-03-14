
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "ref_part.h"
#include "ref_mpi.h"
#include "ref_endian.h"

#include "ref_malloc.h"

#include "ref_migrate.h"

#include "ref_export.h"

#include "ref_twod.h"

REF_STATUS ref_part_b8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL swapped_dbl;
  REF_INT part;
  REF_INT n;
  REF_DBL *xyz;

  long conn_offset, faceid_offset;

  REF_GRID ref_grid;
  REF_NODE ref_node;

  REF_BOOL instrument = REF_FALSE;

  if (instrument) ref_mpi_stopwatch_start();

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  /* header */

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"r");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      RES( 1, fread( &nnode, sizeof(REF_INT), 1, file ), "nnode" );
      RES( 1, fread( &ntri, sizeof(REF_INT), 1, file ), "ntri" );
      RES( 1, fread( &nqua, sizeof(REF_INT), 1, file ), "nqua" );
      RES( 1, fread( &ntet, sizeof(REF_INT), 1, file ), "ntet" );
      RES( 1, fread( &npyr, sizeof(REF_INT), 1, file ), "npyr" );
      RES( 1, fread( &npri, sizeof(REF_INT), 1, file ), "npri" );
      RES( 1, fread( &nhex, sizeof(REF_INT), 1, file ), "nhex" );

      SWAP_INT(nnode);
      SWAP_INT(ntri);
      SWAP_INT(nqua);
      SWAP_INT(ntet);
      SWAP_INT(npyr);
      SWAP_INT(npri);
      SWAP_INT(nhex);
    }

  RSS( ref_mpi_bcast( &nnode, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &ntri, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &nqua, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &ntet, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &npyr, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &npri, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &nhex, 1, REF_INT_TYPE ), "bcast" ); 

  /* guess twod status */

  if ( 0 == ntet && 0 == npyr && 0 != npri && 0 == nhex )
    ref_grid_twod(ref_grid) = REF_TRUE;

  /* nodes */

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init nnodesg");

  if ( ref_mpi_master )
    {
      part = 0;
      for (node=0;node<ref_part_first( nnode, ref_mpi_n, 1 ); node++)
	{
	  RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
	  ref_node_part(ref_node,new_node) = ref_mpi_id;
	  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "x" );
	  SWAP_DBL(swapped_dbl);
	  ref_node_xyz( ref_node, 0, new_node ) = swapped_dbl;
	  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "y" );
	  SWAP_DBL(swapped_dbl);
	  ref_node_xyz( ref_node, 1, new_node ) = swapped_dbl;
	  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "z" );
	  SWAP_DBL(swapped_dbl);
	  ref_node_xyz( ref_node, 2, new_node ) = swapped_dbl;
	}
      for ( part = 1; part<ref_mpi_n ; part++ )
	{
	  n = ref_part_first( nnode, ref_mpi_n, part+1 )
            - ref_part_first( nnode, ref_mpi_n, part );
	  RSS( ref_mpi_send( &n, 1, REF_INT_TYPE, part ), "send" );
	  if ( n > 0 )
	    {
	      ref_malloc( xyz, 3*n, REF_DBL);
	      for (node=0;node<n; node++)
		{
		  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "x" );
		  SWAP_DBL(swapped_dbl);
		  xyz[0+3*node] = swapped_dbl;
		  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "y" );
		  SWAP_DBL(swapped_dbl);
		  xyz[1+3*node] = swapped_dbl;
		  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "z" );
		  SWAP_DBL(swapped_dbl);
		  xyz[2+3*node] = swapped_dbl;
		}
	      RSS( ref_mpi_send( xyz, 3*n, REF_DBL_TYPE, part ), "send" );
	      free(xyz);
	    }
	}
    }
  else
    {
      RSS( ref_mpi_recv( &n, 1, REF_INT_TYPE, 0 ), "recv" );
      if ( n > 0 )
	{
	  ref_malloc( xyz, 3*n, REF_DBL);
	  RSS( ref_mpi_recv( xyz, 3*n, REF_DBL_TYPE, 0 ), "recv" );
	  for (node=0;node<n; node++)
	    {
	      RSS( ref_node_add(ref_node, 
				node+ref_part_first( nnode, 
						     ref_mpi_n, 
						     ref_mpi_id ),
				&new_node ), "new_node");
	      ref_node_part(ref_node,new_node) = ref_mpi_id;
	      ref_node_xyz( ref_node, 0, new_node ) = xyz[0+3*node];
	      ref_node_xyz( ref_node, 1, new_node ) = xyz[1+3*node];
	      ref_node_xyz( ref_node, 2, new_node ) = xyz[2+3*node];
	    }
	  free(xyz);
	}

    }

  if (instrument) ref_mpi_stopwatch_stop("nodes");

  if ( 0 < ntri )
    {
      conn_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode;
      faceid_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)3*(long)ntri
	+ (long)4*(long)4*(long)nqua;
      RSS( ref_part_b8_ugrid_cell( ref_grid_tri(ref_grid), ntri, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "tri" );
    }

  if ( 0 < nqua )
    {
      conn_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)3*(long)ntri;
      faceid_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)4*(long)ntri
	+ (long)4*(long)4*(long)nqua;
      RSS( ref_part_b8_ugrid_cell( ref_grid_qua(ref_grid), nqua, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "qua" );
    }

  if (instrument) ref_mpi_stopwatch_stop("bound");

  if ( 0 < ntet )
    {
      conn_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)4*(long)ntri
	+ (long)4*(long)5*(long)nqua;
      faceid_offset = (long)REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_tet(ref_grid), ntet, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "tet" );
    }

  if ( 0 < npyr )
    {
      conn_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)4*(long)ntri
	+ (long)4*(long)5*(long)nqua
	+ (long)4*(long)4*(long)ntet;
      faceid_offset = (long)REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_pyr(ref_grid), npyr, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "pyr" );
    }

  if ( 0 < npri )
    {
      conn_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)4*(long)ntri
	+ (long)4*(long)5*(long)nqua
	+ (long)4*(long)4*(long)ntet
	+ (long)4*(long)5*(long)npyr;
      faceid_offset = (long)REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_pri(ref_grid), npri, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "pri" );
    }

  if ( 0 < nhex )
    {
      conn_offset = (long)4*(long)7
	+ (long)8*(long)3*(long)nnode
	+ (long)4*(long)4*(long)ntri
	+ (long)4*(long)5*(long)nqua
	+ (long)4*(long)4*(long)ntet
	+ (long)4*(long)5*(long)npyr
	+ (long)4*(long)6*(long)npri;
      faceid_offset = REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_hex(ref_grid), nhex, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "hex" );
    }

  if ( ref_mpi_master ) fclose(file);

  /* ghost xyz */

  RSS( ref_node_ghost_real( ref_node ), "ghost real");

  if (instrument) ref_mpi_stopwatch_stop("volume");

  return REF_SUCCESS;
}

REF_STATUS ref_part_b8_ugrid_cell( REF_CELL ref_cell, REF_INT ncell,
				   REF_NODE ref_node, REF_INT nnode,
				   FILE *file, 
				   long conn_offset, 
				   long faceid_offset )
{
  REF_INT ncell_read;
  REF_INT chunk;
  REF_INT end_of_message = REF_EMPTY;
  REF_INT elements_to_receive;
  REF_INT *c2n;
  REF_INT *tag;
  REF_INT *c2t;
  REF_INT *sent_c2n;
  REF_INT *dest;
  REF_INT *sent_part;
  REF_INT *elements_to_send;
  REF_INT *start_to_send;
  REF_INT node_per, size_per;
  REF_INT section_size;
  REF_INT cell;
  REF_INT part, node;
  REF_INT ncell_keep;
  REF_INT new_location;

  chunk = MAX(1000000, ncell/ref_mpi_n);

  size_per = ref_cell_size_per(ref_cell);
  node_per = ref_cell_node_per(ref_cell);

  ref_malloc( sent_c2n, size_per*chunk, REF_INT );

  if ( ref_mpi_master )
    {

      ref_malloc( elements_to_send, ref_mpi_n, REF_INT );
      ref_malloc( start_to_send, ref_mpi_n, REF_INT );
      ref_malloc( c2n, size_per*chunk, REF_INT );
      ref_malloc( tag, chunk, REF_INT );
      ref_malloc( c2t, size_per*chunk, REF_INT );
      ref_malloc( dest, chunk, REF_INT );

      ncell_read = 0;
      while ( ncell_read < ncell )
	{
	  section_size = MIN(chunk,ncell-ncell_read);

	  REIS(0,
	       fseek(file,conn_offset+(long)4*(long)node_per*(long)ncell_read,
		     SEEK_SET), "seek conn failed");
	  RES((size_t)(section_size*node_per),
	      fread( c2n,
		     sizeof(REF_INT), section_size*node_per, file ), "cn" );
	  for (cell=0;cell<section_size*node_per;cell++)
	    SWAP_INT(c2n[cell]);
	  for (cell=0;cell<section_size*node_per;cell++)
	    c2n[cell]--;
	  if ( ref_cell_last_node_is_an_id(ref_cell) )
	    {
	      REIS(0,fseek(file,faceid_offset+(long)4*(long)ncell_read,
			   SEEK_SET),"seek tag failed");
	      RES((size_t)(section_size),
		  fread( tag,
			 sizeof(REF_INT), section_size, file ), "tag" );
	      for (cell=0;cell<section_size;cell++)
		SWAP_INT(tag[cell]);

	      /* sort into right locations */
	      for (cell=0;cell<section_size;cell++)
		for (node=0;node<node_per;node++)
		  c2t[node+cell*size_per] = c2n[node+cell*node_per];
	      for (cell=0;cell<section_size;cell++)
		c2t[node_per+cell*size_per] = tag[cell];
	      for (cell=0;cell<section_size*size_per;cell++)
		c2n[cell]=c2t[cell];
	    }

	  ncell_read += section_size;

	  for (cell=0;cell<section_size;cell++)
	    dest[cell] = ref_part_implicit( nnode, ref_mpi_n, 
					    c2n[size_per*cell] );

	  for ( part = 0; part< ref_mpi_n;part++ )
	    elements_to_send[part] = 0;
	  for (cell=0;cell<section_size;cell++)
	    elements_to_send[dest[cell]]++;

	  start_to_send[0]=0;
	  for ( part = 0; part< ref_mpi_n-1;part++ )
	    start_to_send[part+1] = start_to_send[part]+elements_to_send[part];

	  for ( part = 0; part< ref_mpi_n;part++ )
	    elements_to_send[part] = 0;
	  for (cell=0;cell<section_size;cell++)
	    {
	      new_location = start_to_send[dest[cell]]
		+ elements_to_send[dest[cell]];
	      for (node=0;node<size_per;node++)
		sent_c2n[node+size_per*new_location] = c2n[node+size_per*cell];
	      elements_to_send[dest[cell]]++;
	    }

	  /* master keepers */

	  ncell_keep = elements_to_send[0];
	  if ( 0 < ncell_keep )
	    { 
	      ref_malloc_init( sent_part, size_per*ncell_keep, 
			       REF_INT, REF_EMPTY );

	      for (cell=0;cell<ncell_keep;cell++)
		for (node=0;node<node_per;node++)
		  sent_part[node+size_per*cell] =
		    ref_part_implicit( nnode, ref_mpi_n, 
				       sent_c2n[node+size_per*cell] );

	      RSS( ref_cell_add_many_global( ref_cell, ref_node,
					     ncell_keep, 
					     sent_c2n, sent_part,
					     ref_mpi_id ),"many glob");

	      ref_free(sent_part);
	    }

	  for ( part = 1; part< ref_mpi_n;part++ )
	    if ( 0 < elements_to_send[part] ) 
	      {
		RSS( ref_mpi_send( &(elements_to_send[part]), 
				   1, REF_INT_TYPE, part ), "send" );
		RSS( ref_mpi_send( &(sent_c2n[size_per*start_to_send[part]]), 
				   size_per*elements_to_send[part], 
				   REF_INT_TYPE, part ), "send" );
	      }

	}

      ref_free(dest);
      ref_free(c2t);
      ref_free(tag);
      ref_free(c2n);
      ref_free(start_to_send);
      ref_free(elements_to_send);

      /* signal we are done */
      for ( part = 1; part<ref_mpi_n ; part++ )
	RSS( ref_mpi_send( &end_of_message, 1, REF_INT_TYPE, part ), "send" );

    }  
  else
    {
      do {
	RSS( ref_mpi_recv( &elements_to_receive, 1, REF_INT_TYPE, 0 ), "recv" );
	if ( elements_to_receive > 0 )
	  {
	    RSS( ref_mpi_recv( sent_c2n, size_per*elements_to_receive, 
			       REF_INT_TYPE, 0 ), "send" );

	    ref_malloc_init( sent_part, size_per*elements_to_receive, 
			     REF_INT, REF_EMPTY );

	    for (cell=0;cell<elements_to_receive;cell++)
	      for (node=0;node<node_per;node++)
		sent_part[node+size_per*cell] =
		  ref_part_implicit( nnode, ref_mpi_n, 
				     sent_c2n[node+size_per*cell] );

	    RSS( ref_cell_add_many_global( ref_cell, ref_node,
					   elements_to_receive, 
					   sent_c2n, sent_part,
					   ref_mpi_id ), "many glob");

	    ref_free( sent_part );
	  }
      } while ( elements_to_receive != end_of_message );
    }

  free(sent_c2n);

  RSS( ref_migrate_shufflin_cell( ref_node, ref_cell ), "fill ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_part_metric( REF_NODE ref_node, char *filename )
{
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_INT nnode_read, section_size;
  REF_INT node, local, global, im;
  size_t end_of_string;
  REF_BOOL sol_format, found_keyword;
  REF_INT nnode, ntype, type;
  REF_INT status;
  char line[1024];

  file = NULL;
  sol_format = REF_FALSE;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"r");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      end_of_string = strlen(filename);
      if( strcmp(&filename[end_of_string-4],".sol") == 0 ) 
	{
	  sol_format = REF_TRUE;
	  found_keyword = REF_FALSE;
	  while (!feof(file))
	    {
	      status = fscanf( file, "%s", line);
	      if ( EOF == status ) break;
	      REIS( 1, status, "line read failed");
	      
	      if ( 0 == strcmp("SolAtVertices",line))
		{
		  REIS( 1, fscanf(file, "%d", &nnode), "read nnode" );
		  REIS( ref_node_n_global(ref_node), nnode, 
			"wrong vertex number in .sol" );
		  REIS( 2, fscanf(file, "%d %d", 
				  &ntype, &type),"read header" );
		  REIS( 1, ntype, "expected one type in .sol" );
		  REIS( 3, type, "expected type GmfSymMat in .sol" );
		  fscanf(file, "%*[^1234567890-+.]"); /* remove blank line */
		  found_keyword = REF_TRUE;
		  break;
		}
	    }
	  RAS(found_keyword,"SolAtVertices keyword missing from .sol metric");
	}
    }

  chunk = MAX(100000, ref_node_n_global(ref_node)/ref_mpi_n);
  chunk = MIN( chunk, ref_node_n_global(ref_node) );

  ref_malloc_init( metric, 6*chunk, REF_DBL, -1.0 );
  
  nnode_read = 0;
  while ( nnode_read < ref_node_n_global(ref_node) )
    {
      section_size = MIN(chunk,ref_node_n_global(ref_node)-nnode_read);
      if ( ref_mpi_master )
	{
	  for (node=0;node<section_size;node++)
	  if ( sol_format ) 
	    {
	      REIS( 6, fscanf( file, "%lf %lf %lf %lf %lf %lf", 
			       &(metric[0+6*node]),
			       &(metric[1+6*node]),
			       &(metric[3+6*node]), /* transposed 3,2 */
			       &(metric[2+6*node]),
			       &(metric[4+6*node]),
			       &(metric[5+6*node]) ), "metric read error" );
	    }
	  else
	    {
	      REIS( 6, fscanf( file, "%lf %lf %lf %lf %lf %lf", 
			       &(metric[0+6*node]),
			       &(metric[1+6*node]),
			       &(metric[2+6*node]),
			       &(metric[3+6*node]),
			       &(metric[4+6*node]),
			       &(metric[5+6*node]) ), "metric read error" );
	    }
	  RSS( ref_mpi_bcast( metric, 6*chunk, REF_DBL_TYPE ), "bcast" );
	}
      else
	{
	  RSS( ref_mpi_bcast( metric, 6*chunk, REF_DBL_TYPE ), "bcast" );
	}
      for (node=0;node<section_size;node++)
	{
	  global = node + nnode_read;
	  RXS( ref_node_local( ref_node, global, &local ), 
	       REF_NOT_FOUND, "local" );
	  if ( REF_EMPTY != local )
	    for(im=0;im<6;im++)
	      ref_node_metric(ref_node,im,local) = metric[im+6*node];
	}
      nnode_read += section_size;
    }

  ref_free( metric );

  return REF_SUCCESS;
}

REF_STATUS ref_part_bamg_metric( REF_GRID ref_grid, char *filename )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_INT file_nnode, nnode_read, section_size;
  REF_INT node, local, opposite, global;
  REF_INT nterm, nnode;

  nnode = ref_node_n_global(ref_node)/2;
  
  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"r");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );
      REIS( 2, fscanf(file, "%d %d", 
		      &file_nnode, &nterm),"read header" );
      REIS( nnode, file_nnode, "wrong node count" );
      REIS( 3, nterm, "expected 3 term 2x2 anisotropic M" );
    }

  chunk = MAX(100000, nnode/ref_mpi_n);
  chunk = MIN( chunk, nnode );

  ref_malloc_init( metric, 3*chunk, REF_DBL, -1.0 );
  
  nnode_read = 0;
  while ( nnode_read < nnode )
    {
      section_size = MIN(chunk,ref_node_n_global(ref_node)-nnode_read);
      if ( ref_mpi_master )
	{
	  for (node=0;node<section_size;node++)
	    {
	      REIS( 3, fscanf( file, "%lf %lf %lf", 
			       &(metric[0+3*node]),
			       &(metric[1+3*node]),
			       &(metric[2+3*node]) ), "metric read error" );
	    }
	  RSS( ref_mpi_bcast( metric, 3*chunk, REF_DBL_TYPE ), "bcast" );
	}
      else
	{
	  RSS( ref_mpi_bcast( metric, 3*chunk, REF_DBL_TYPE ), "bcast" );
	}
      for (node=0;node<section_size;node++)
	{
	  global = node + nnode_read;
	  RXS( ref_node_local( ref_node, global, &local ), 
	       REF_NOT_FOUND, "local" );
	  if ( REF_EMPTY != local )
	    {
	      ref_node_metric(ref_node,0,local) = metric[0+3*node];
	      ref_node_metric(ref_node,1,local) = 0.0;
	      ref_node_metric(ref_node,2,local) = metric[1+3*node];
	      ref_node_metric(ref_node,3,local) = 1.0;
	      ref_node_metric(ref_node,4,local) = 0.0;
	      ref_node_metric(ref_node,5,local) = metric[2+3*node];
	      RSS( ref_twod_opposite_node( ref_grid_pri(ref_grid),
					   local, &opposite ),
		   "opposite twod node on other plane missing" );
	      ref_node_metric(ref_node,0,opposite) = metric[0+3*node];
	      ref_node_metric(ref_node,1,opposite) = 0.0;
	      ref_node_metric(ref_node,2,opposite) = metric[1+3*node];
	      ref_node_metric(ref_node,3,opposite) = 1.0;
	      ref_node_metric(ref_node,4,opposite) = 0.0;
	      ref_node_metric(ref_node,5,opposite) = metric[2+3*node];
	    }
	}
      nnode_read += section_size;
    }

  ref_free( metric );

  return REF_SUCCESS;
}

REF_STATUS ref_part_ratio( REF_NODE ref_node, REF_DBL *ratio, char *filename )
{
  FILE *file;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT nnode_read, section_size;
  REF_INT node, local, global;

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"r");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );
    }

  chunk = MAX(100000, ref_node_n_global(ref_node)/ref_mpi_n);
  chunk = MIN( chunk, ref_node_n_global(ref_node) );

  ref_malloc_init( data, chunk, REF_DBL, -1.0 );
  
  nnode_read = 0;
  while ( nnode_read < ref_node_n_global(ref_node) )
    {
      section_size = MIN(chunk,ref_node_n_global(ref_node)-nnode_read);
      if ( ref_mpi_master )
	{
	  for (node=0;node<section_size;node++)
	    REIS( 1, fscanf( file, "%lf", 
			     &(data[node]) ), "data ratio read error" );
	  RSS( ref_mpi_bcast( data, chunk, REF_DBL_TYPE ), "bcast" );
	}
      else
	{
	  RSS( ref_mpi_bcast( data, chunk, REF_DBL_TYPE ), "bcast" );
	}
      for (node=0;node<section_size;node++)
	{
	  global = node + nnode_read;
	  RXS( ref_node_local( ref_node, global, &local ), 
	       REF_NOT_FOUND, "local" );
	  if ( REF_EMPTY != local )
	    ratio[local] = data[node];
	}
      nnode_read += section_size;
    }

  ref_free( data );

  return REF_SUCCESS;
}

