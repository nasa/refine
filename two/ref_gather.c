
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_gather.h"
#include "ref_export.h"

#include "ref_endian.h"
#include "ref_malloc.h"
#include "ref_mpi.h"

REF_STATUS ref_gather_plot( REF_GRID ref_grid, const char *filename  )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);

  RSS( ref_node_synchronize_globals( ref_node ), "sync" );

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen("ref_gather_script.m","w");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      fprintf(file, "filename='%s';\n",filename);
      fprintf(file, "nodes=[\n");
    }

  RSS( ref_gather_node_tec_part( ref_node, file ), "nodes");

  if ( ref_mpi_master )
    {
      fprintf(file, "];\n");
      fprintf(file, "elements=[\n");
    }

  RSS( ref_gather_cell_tec( ref_node, ref_grid_tri(ref_grid), file ), "nodes");

  if ( ref_mpi_master )
    {
      fprintf(file, "];\n");
      fprintf(file, "x=[];z=[];\n");
      fprintf(file, "for elem=1:size(elements,1)\n");
      fprintf(file, "  x=[ x\n");
      fprintf(file, "      nodes(elements(elem,1),1)\n");
      fprintf(file, "      nodes(elements(elem,2),1)\n");
      fprintf(file, "      nodes(elements(elem,3),1)\n");
      fprintf(file, "      nodes(elements(elem,1),1)\n");
      fprintf(file, "      NaN\n");
      fprintf(file, "    ];\n");
      fprintf(file, "  z=[ z\n");
      fprintf(file, "      nodes(elements(elem,1),3)\n");
      fprintf(file, "      nodes(elements(elem,2),3)\n");
      fprintf(file, "      nodes(elements(elem,3),3)\n");
      fprintf(file, "      nodes(elements(elem,1),3)\n");
      fprintf(file, "      NaN\n");
      fprintf(file, "    ];\n");
      fprintf(file, "end\n");
      fprintf(file, "axis('square');\n");
      fprintf(file, "plot(x,z);\n");
      fprintf(file, "print(filename,'-deps');\n");
      fprintf(file, "\n");
      fclose(file);

      system("octave -q ref_gather_script.m && rm ref_gather_script.m");

    }

  return REF_SUCCESS;
}

static REF_BOOL recording_movie = REF_FALSE;

REF_STATUS ref_gather_tec_movie_record_button( REF_BOOL on_or_off )
{
  recording_movie = on_or_off;
  return REF_SUCCESS;
}

static FILE *movie_file = NULL;
static REF_DBL movie_time = 0.0;

REF_STATUS ref_gather_tec_movie_frame( REF_GRID ref_grid,
				       const char *zone_title )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nnode,ntri;

  if ( !recording_movie ) return REF_SUCCESS;

  RSS( ref_node_synchronize_globals( ref_node ), "sync" );

  nnode = ref_node_n_global(ref_node);

  RSS( ref_gather_ncell( ref_node, ref_grid_tri(ref_grid), &ntri ), "ntri");

  if ( ref_mpi_master )
    {
      if ( NULL == (void *)movie_file )
	{ 
	  movie_file = fopen("ref_gather_movie.tec","w");
	  if ( NULL == (void *)movie_file ) 
	    printf("unable to open ref_gather_movie.tec\n");
	  RNS(movie_file, "unable to open file" );
	  
	  fprintf(movie_file, "title=\"tecplot refine partion file\"\n");
	  fprintf(movie_file, "variables = \"x\" \"y\" \"z\" \"p\" \"a\"\n");
	}
      if ( NULL == zone_title )
	{
	  fprintf(movie_file,
		  "zone t=\"part\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s, solutiontime=%f\n",
		  nnode, ntri, "point", "fetriangle", movie_time );
	}
      else
	{
	  fprintf(movie_file,
		  "zone t=\"%s\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s, solutiontime=%f\n",
		  zone_title, nnode, ntri, "point", "fetriangle", movie_time );
	}
      movie_time += 1.0;
    }

  RSS( ref_gather_node_tec_part( ref_node, movie_file ), "nodes");
  RSS( ref_gather_cell_tec( ref_node, ref_grid_tri(ref_grid), movie_file ), "nodes");

  return REF_SUCCESS;
}

REF_STATUS ref_gather_tec_part( REF_GRID ref_grid, const char *filename  )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nnode,ntri;

  RSS( ref_node_synchronize_globals( ref_node ), "sync" );

  nnode = ref_node_n_global(ref_node);

  RSS( ref_gather_ncell( ref_node, ref_grid_tri(ref_grid), &ntri ), "ntri");

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"w");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      fprintf(file, "title=\"tecplot refine partion file\"\n");
      fprintf(file, "variables = \"x\" \"y\" \"z\" \"p\" \"a\"\n");
      fprintf(file,
	"zone t=part, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	      nnode, ntri, "point", "fetriangle" );
    }

  RSS( ref_gather_node_tec_part( ref_node, file ), "nodes");
  RSS( ref_gather_cell_tec( ref_node, ref_grid_tri(ref_grid), file ), "nodes");

  if ( ref_mpi_master ) fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_meshb( REF_GRID ref_grid, const char *filename  )
{
  REF_BOOL verbose = REF_TRUE;
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT code, version, dim;
  int next_position, keyword_code;
  REF_BOOL swap_endian = REF_FALSE;
  REF_BOOL has_id = REF_TRUE;
  
  RAS( !ref_grid_twod(ref_grid), "only 3D" );
  
  RSS( ref_node_synchronize_globals( ref_node ), "sync" );
  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"w");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      code = 1;
      REIS(1, fwrite(&code,sizeof(int),1,file),"code");
      version = 2;
      REIS(1, fwrite(&version,sizeof(int),1,file),"version");
      next_position = 4+4+4+ftell(file);
      keyword_code = 3;
      REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"dim code");
      REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
      dim = 3;
      REIS(1, fwrite(&dim,sizeof(int),1,file),"dim");
    }

  if ( ref_mpi_master )
    {
      next_position = 4+4+4+ref_node_n_global(ref_node)*(3*8+4)+ftell(file);
      keyword_code = 4;
      REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
      REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
      REIS(1, fwrite(&(ref_node_n(ref_node)),sizeof(int),1,file),"nnode");
      if (verbose) printf("vertex kw %d next %d n %d\n",
			  keyword_code,next_position,
			  ref_node_n_global(ref_node));

    }
  RSS( ref_gather_node( ref_node, swap_endian, has_id, file ), "nodes");
  if ( ref_mpi_master )
    REIS( next_position, ftell(file), "vertex inconsistent");

  return REF_SUCCESS;
}

REF_STATUS ref_gather_b8_ugrid( REF_GRID ref_grid, const char *filename  )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nnode,ntri,nqua,ntet,npyr,npri,nhex;
  REF_CELL ref_cell;
  REF_INT group;
  REF_INT faceid, min_faceid, max_faceid;
  REF_BOOL swap_endian = REF_TRUE;
  REF_BOOL has_id = REF_FALSE;

  RSS( ref_node_synchronize_globals( ref_node ), "sync" );

  nnode = ref_node_n_global(ref_node);

  RSS( ref_gather_ncell( ref_node, ref_grid_tri(ref_grid), &ntri ), "ntri");
  RSS( ref_gather_ncell( ref_node, ref_grid_qua(ref_grid), &nqua ), "nqua");

  RSS( ref_gather_ncell( ref_node, ref_grid_tet(ref_grid), &ntet ), "ntet");
  RSS( ref_gather_ncell( ref_node, ref_grid_pyr(ref_grid), &npyr ), "npyr");
  RSS( ref_gather_ncell( ref_node, ref_grid_pri(ref_grid), &npri ), "npri");
  RSS( ref_gather_ncell( ref_node, ref_grid_hex(ref_grid), &nhex ), "nhex");
  
  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"w");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      SWAP_INT(nnode);
      SWAP_INT(ntri);
      SWAP_INT(nqua);
      SWAP_INT(ntet);
      SWAP_INT(npyr);
      SWAP_INT(npri);
      SWAP_INT(nhex);

      REIS(1, fwrite(&nnode,sizeof(REF_INT),1,file),"nnode");

      REIS(1, fwrite(&ntri,sizeof(REF_INT),1,file),"ntri");
      REIS(1, fwrite(&nqua,sizeof(REF_INT),1,file),"nqua");

      REIS(1, fwrite(&ntet,sizeof(REF_INT),1,file),"ntet");
      REIS(1, fwrite(&npyr,sizeof(REF_INT),1,file),"npyr");
      REIS(1, fwrite(&npri,sizeof(REF_INT),1,file),"npri");
      REIS(1, fwrite(&nhex,sizeof(REF_INT),1,file),"nhex");
    }

  RSS( ref_gather_node( ref_node, swap_endian, has_id, file ), "nodes");

  RSS( ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid), "range");

  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    RSS( ref_gather_cell( ref_node,ref_grid_tri(ref_grid), 
			  REF_FALSE, faceid, file ), "tri c2n");
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    RSS( ref_gather_cell( ref_node,ref_grid_qua(ref_grid), 
			  REF_FALSE, faceid, file ), "qua c2n");

  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    RSS( ref_gather_cell( ref_node,ref_grid_tri(ref_grid), 
			  REF_TRUE, faceid, file ), "tri faceid");
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    RSS( ref_gather_cell( ref_node,ref_grid_qua(ref_grid), 
			  REF_TRUE, faceid, file ), "qua faceid");
  faceid = REF_EMPTY;
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    RSS( ref_gather_cell( ref_node, ref_cell, 
			  REF_FALSE, faceid, file ), "cell c2n");

  if ( ref_mpi_master ) fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_metric( REF_GRID ref_grid, const char *filename  )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);

  RSS( ref_node_synchronize_globals( ref_node ), "sync" );

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"w");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );
    }

  RSS( ref_gather_node_metric( ref_node, file ), "nodes");

  if ( ref_mpi_master ) fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_ncell( REF_NODE ref_node, REF_CELL ref_cell, 
			     REF_INT *ncell )
{
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell_local;

  ncell_local = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) )
      ncell_local++;

  RSS( ref_mpi_sum( &ncell_local, ncell, 1, REF_INT_TYPE ), "sum");

  return REF_SUCCESS;
}

REF_STATUS ref_gather_node( REF_NODE ref_node,
			    REF_BOOL swap_endian, REF_BOOL has_id, FILE *file )
{
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_DBL swapped_dbl;
  REF_INT nnode_written, first, n, i;
  REF_INT global, local;
  REF_INT id = 0;
  REF_STATUS status;

  chunk = ref_node_n_global(ref_node)/ref_mpi_n + 1;

  ref_malloc( local_xyzm, 4*chunk, REF_DBL );
  ref_malloc( xyzm, 4*chunk, REF_DBL );

  nnode_written = 0;
  while ( nnode_written < ref_node_n_global(ref_node) )
    {

      first = nnode_written;
      n = MIN( chunk, ref_node_n_global(ref_node)-nnode_written );

      nnode_written += n;

      for (i=0;i<4*chunk;i++)
	local_xyzm[i] = 0.0;

      for (i=0;i<n;i++)
	{
	  global = first + i;
	  status = ref_node_local( ref_node, global, &local );
	  RXS( status, REF_NOT_FOUND, "node local failed" );
	  if ( REF_SUCCESS == status &&
	       ref_mpi_id == ref_node_part(ref_node,local) )
	    {
	      local_xyzm[0+4*i] = ref_node_xyz(ref_node,0,local);
	      local_xyzm[1+4*i] = ref_node_xyz(ref_node,1,local);
	      local_xyzm[2+4*i] = ref_node_xyz(ref_node,2,local);
	      local_xyzm[3+4*i] = 1.0;
	    }
	  else
	    {
	      local_xyzm[0+4*i] = 0.0;
	      local_xyzm[1+4*i] = 0.0;
	      local_xyzm[2+4*i] = 0.0;
	      local_xyzm[3+4*i] = 0.0;
	    }
	}

      RSS( ref_mpi_sum( local_xyzm, xyzm, 4*n, REF_DBL_TYPE ), "sum" );
      
      if ( ref_mpi_master )
	for ( i=0; i<n; i++ )
	  {
	    if ( ABS( xyzm[3+4*i] - 1.0 ) > 0.1 )
	      {
		printf("error gather node %d %f\n",first+i, xyzm[3+4*i]);
	      }
	    swapped_dbl = xyzm[0+4*i]; if (swap_endian) SWAP_DBL(swapped_dbl);
	    REIS(1, fwrite(&swapped_dbl,sizeof(REF_DBL),1,file),"x");
	    swapped_dbl = xyzm[1+4*i]; if (swap_endian) SWAP_DBL(swapped_dbl);
	    REIS(1, fwrite(&swapped_dbl,sizeof(REF_DBL),1,file),"y");
	    swapped_dbl = xyzm[2+4*i]; if (swap_endian) SWAP_DBL(swapped_dbl);
	    REIS(1, fwrite(&swapped_dbl,sizeof(REF_DBL),1,file),"z");
	    if (has_id) REIS(1, fwrite(&id,sizeof(REF_INT),1,file),"id");
	  }
    }

  ref_free( xyzm );
  ref_free( local_xyzm );

  return REF_SUCCESS;
}
REF_STATUS ref_gather_node_tec_part( REF_NODE ref_node, FILE *file )
{
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_INT nnode_written, first, n, i;
  REF_INT global, local;
  REF_STATUS status;
  REF_INT dim=6;

  chunk = ref_node_n_global(ref_node)/ref_mpi_n + 1;

  ref_malloc( local_xyzm, dim*chunk, REF_DBL );
  ref_malloc( xyzm, dim*chunk, REF_DBL );

  nnode_written = 0;
  while ( nnode_written < ref_node_n_global(ref_node) )
    {
      first = nnode_written;
      n = MIN( chunk, ref_node_n_global(ref_node)-nnode_written );

      nnode_written += n;

      for (i=0;i<dim*chunk;i++)
	local_xyzm[i] = 0.0;

      for (i=0;i<n;i++)
	{
	  global = first + i;
	  status = ref_node_local( ref_node, global, &local );
	  RXS( status, REF_NOT_FOUND, "node local failed" );
	  if ( REF_SUCCESS == status &&
	       ref_mpi_id == ref_node_part(ref_node,local) )
	    {
	      local_xyzm[0+dim*i] = ref_node_xyz(ref_node,0,local);
	      local_xyzm[1+dim*i] = ref_node_xyz(ref_node,1,local);
	      local_xyzm[2+dim*i] = ref_node_xyz(ref_node,2,local);
	      local_xyzm[3+dim*i] = (REF_DBL)ref_node_part(ref_node,local);
	      local_xyzm[4+dim*i] = (REF_DBL)ref_node_age(ref_node,local);
	      local_xyzm[5+dim*i] = 1.0;
	    }
	  else
	    {
	      local_xyzm[0+dim*i] = 0.0;
	      local_xyzm[1+dim*i] = 0.0;
	      local_xyzm[2+dim*i] = 0.0;
	      local_xyzm[3+dim*i] = 0.0;
	      local_xyzm[4+dim*i] = 0.0;
	      local_xyzm[5+dim*i] = 0.0;
	    }
	}

      for (i=0;i<n;i++)
	if ( (ABS( local_xyzm[5+dim*i] - 1.0 ) > 0.1) &&
	     (ABS( local_xyzm[5+dim*i] - 0.0 ) > 0.1) )
	  {
	    printf("error gather node before sum %d %f\n",
		   first+i, local_xyzm[5+dim*i]);
	  }

      RSS( ref_mpi_sum( local_xyzm, xyzm, dim*n, REF_DBL_TYPE ), "sum" );

      if ( ref_mpi_master )
	for ( i=0; i<n; i++ )
	  {
	    if ( ABS( xyzm[5+dim*i] - 1.0 ) > 0.1 )
	      {
		printf("error gather node %d %f\n",first+i, xyzm[5+dim*i]);
	      }
	    fprintf(file,"%.15e %.15e %.15e %.0f %.0f\n",
		    xyzm[0+dim*i], xyzm[1+dim*i], xyzm[2+dim*i], 
		    xyzm[3+dim*i], xyzm[4+dim*i]);
	  }
    }

  ref_free( xyzm );
  ref_free( local_xyzm );

  return REF_SUCCESS;
}

REF_STATUS ref_gather_node_metric( REF_NODE ref_node, FILE *file )
{
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_INT nnode_written, first, n, i, im;
  REF_INT global, local;
  REF_STATUS status;

  chunk = ref_node_n_global(ref_node)/ref_mpi_n + 1;

  ref_malloc( local_xyzm, 7*chunk, REF_DBL );
  ref_malloc( xyzm, 7*chunk, REF_DBL );

  nnode_written = 0;
  while ( nnode_written < ref_node_n_global(ref_node) )
    {
      first = nnode_written;
      n = MIN( chunk, ref_node_n_global(ref_node)-nnode_written );

      nnode_written += n;

      for (i=0;i<7*chunk;i++)
	local_xyzm[i] = 0.0;

      for (i=0;i<n;i++)
	{
	  global = first + i;
	  status = ref_node_local( ref_node, global, &local );
	  RXS( status, REF_NOT_FOUND, "node local failed" );
	  if ( REF_SUCCESS == status &&
	       ref_mpi_id == ref_node_part(ref_node,local) )
	    {
	      for (im=0;im<6;im++)
		local_xyzm[im+7*i] = ref_node_metric(ref_node,im,local);
	      local_xyzm[6+7*i] = 1.0;
	    }
	  else
	    {
	      for (im=0;im<7;im++)
		local_xyzm[im+7*i] = 0.0;
	    }
	}

      RSS( ref_mpi_sum( local_xyzm, xyzm, 7*n, REF_DBL_TYPE ), "sum" );

      if ( ref_mpi_master )
	for ( i=0; i<n; i++ )
	  {
	    if ( ABS( xyzm[6+7*i] - 1.0 ) > 0.1 )
	      {
		printf("error gather node %d %f\n",first+i, xyzm[6+7*i]);
	      }
	    fprintf(file,"%.15e %.15e %.15e %.15e %.15e %.15e \n",
		    xyzm[0+7*i], xyzm[1+7*i], xyzm[2+7*i], 
		    xyzm[3+7*i], xyzm[4+7*i], xyzm[5+7*i]);
	  }
    }

  ref_free( xyzm );
  ref_free( local_xyzm );

  return REF_SUCCESS;
}

REF_STATUS ref_gather_cell( REF_NODE ref_node, REF_CELL ref_cell, 
			    REF_BOOL faceid_insted_of_c2n, REF_INT faceid,
			    FILE *file )
{
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per = ref_cell_node_per(ref_cell);
  REF_INT size_per = ref_cell_size_per(ref_cell);
  REF_INT ncell;
  REF_INT *c2n;
  REF_INT proc;

  if ( ref_mpi_master )
    {
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
	if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) &&
	     ( ! ref_cell_last_node_is_an_id(ref_cell) ||
	       nodes[ref_cell_node_per(ref_cell)] == faceid ) )
	  {
	    if ( faceid_insted_of_c2n )
	      {
		for ( node = node_per; node < size_per; node++ )
		  {
		    SWAP_INT(nodes[node]);
		    REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),
			 "cel node");
		  }
	      }
	    else
	      {
		for ( node = 0; node < node_per; node++ )
		  {
		    nodes[node] = ref_node_global(ref_node,nodes[node]);
		    nodes[node]++;
		    SWAP_INT(nodes[node]);
		    REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),
			 "cel node");
		  }
	      }
	  }
    }

  if ( ref_mpi_master )
    {
      for (proc=1;proc<ref_mpi_n;proc++)
	{
	  RSS( ref_mpi_recv( &ncell, 1, REF_INT_TYPE, proc ), "recv ncell");
	  ref_malloc(c2n, ncell*size_per, REF_INT);
	  RSS( ref_mpi_recv( c2n, ncell*size_per, 
			     REF_INT_TYPE, proc ), "recv c2n");
	  for ( cell = 0; cell < ncell; cell++ )
	    if ( faceid_insted_of_c2n )
	      {
		for ( node = node_per; node < size_per; node++ )
		  {
		    SWAP_INT(c2n[node+size_per*cell]);
		    REIS(1, fwrite(&(c2n[node+size_per*cell]),
				   sizeof(REF_INT),1,file),"cell");
		  }
	      }
	    else
	      {
		for ( node = 0; node < node_per; node++ )
		  {
		    c2n[node+size_per*cell]++;
		    SWAP_INT(c2n[node+size_per*cell]);
		    REIS(1, fwrite(&(c2n[node+size_per*cell]),
				   sizeof(REF_INT),1,file),"cell");
		  }
	      }
	  ref_free(c2n);
	}
    }
  else
    {
      ncell = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
	if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) &&
	     ( ! ref_cell_last_node_is_an_id(ref_cell) ||
	       nodes[ref_cell_node_per(ref_cell)] == faceid ) )
	  ncell++;
      RSS( ref_mpi_send( &ncell, 1, REF_INT_TYPE, 0 ), "send ncell");
      ref_malloc(c2n, ncell*size_per, REF_INT);
      ncell = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
	if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) &&
	     ( ! ref_cell_last_node_is_an_id(ref_cell) ||
	       nodes[ref_cell_node_per(ref_cell)] == faceid ) )
	  {
	    for ( node = 0; node < node_per; node++ )
	      c2n[node+size_per*ncell] = ref_node_global(ref_node,nodes[node]);
	    for ( node = node_per; node < size_per; node++ )
	      c2n[node+size_per*ncell] = nodes[node];
	    ncell++;
	  }
      RSS( ref_mpi_send( c2n, ncell*size_per, 
			 REF_INT_TYPE, 0 ), "send c2n");

      ref_free(c2n);
    }

  return REF_SUCCESS;
}

REF_STATUS ref_gather_cell_tec( REF_NODE ref_node, REF_CELL ref_cell, 
				FILE *file )
{
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per = ref_cell_node_per(ref_cell);
  REF_INT size_per = ref_cell_size_per(ref_cell);
  REF_INT ncell;
  REF_INT *c2n;
  REF_INT proc;

  if ( ref_mpi_master )
    {
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
	if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) )
	  {
	    for ( node = 0; node < node_per; node++ )
	      {
		nodes[node] = ref_node_global(ref_node,nodes[node]);
		nodes[node]++;
		fprintf(file," %d",nodes[node]);
	      }
	    fprintf(file,"\n");
	  }
    }

  if ( ref_mpi_master )
    {
      for (proc=1;proc<ref_mpi_n;proc++)
	{
	  RSS( ref_mpi_recv( &ncell, 1, REF_INT_TYPE, proc ), "recv ncell");
	  ref_malloc(c2n, ncell*size_per, REF_INT);
	  RSS( ref_mpi_recv( c2n, ncell*size_per, 
			     REF_INT_TYPE, proc ), "recv c2n");
	  for ( cell = 0; cell < ncell; cell++ )
	    {
	      for ( node = 0; node < node_per; node++ )
		{
		  c2n[node+size_per*cell]++;
		  fprintf(file," %d",c2n[node+size_per*cell]);
		}
	      fprintf(file,"\n");
	    }
	  ref_free(c2n);
	}
    }
  else
    {
      ncell = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
	if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) )
	  ncell++;
      RSS( ref_mpi_send( &ncell, 1, REF_INT_TYPE, 0 ), "send ncell");
      ref_malloc(c2n, ncell*size_per, REF_INT);
      ncell = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
	if ( ref_mpi_id == ref_node_part(ref_node,nodes[0]) )
	  {
	    for ( node = 0; node < node_per; node++ )
	      c2n[node+size_per*ncell] = ref_node_global(ref_node,nodes[node]);
	    for ( node = node_per; node < size_per; node++ )
	      c2n[node+size_per*ncell] = nodes[node];
	    ncell++;
	  }
      RSS( ref_mpi_send( c2n, ncell*size_per, 
			 REF_INT_TYPE, 0 ), "send c2n");

      ref_free(c2n);
    }

  return REF_SUCCESS;
}
