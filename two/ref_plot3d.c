
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "ref_plot3d.h"
#include "ref_malloc.h"

#include "ref_endian.h"

/* http://cfl3d.larc.nasa.gov/Cfl3dv6/V5Manual/FileForm.pdf */

REF_STATUS ref_plot3d_from_file( REF_PLOT3D *ref_plot3d_ptr, char *filename )
{
  REF_PLOT3D ref_plot3d;
  FILE *file;
  REF_INT fortran_record_size;
  REF_INT i, j, k, n, ixyz;
  REF_DBL dummy;
  REF_PATCH ref_patch;

  ref_malloc( *ref_plot3d_ptr, 1, REF_PLOT3D_STRUCT );
  ref_plot3d = *ref_plot3d_ptr;

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "srt rec" );
  SWAP_INT(fortran_record_size);
  REIS( 1*4, fortran_record_size, "header start record size" );

  RES( 1, fread( &ref_plot3d_ngrid(ref_plot3d), sizeof(REF_INT), 1,file),"ng" );
  SWAP_INT(ref_plot3d_ngrid(ref_plot3d));

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "end rec" );
  SWAP_INT(fortran_record_size);
  REIS( 1*4, fortran_record_size, "header end record size" );

  ref_malloc( ref_plot3d->patch, ref_plot3d_ngrid(ref_plot3d), 
	      REF_PATCH );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "srt rec" );
  SWAP_INT(fortran_record_size);
  REIS( ref_plot3d_ngrid(ref_plot3d)*3*4, fortran_record_size, 
	"dims start record size" );

  for (n=0;n<ref_plot3d_ngrid(ref_plot3d);n++)
    {
      ref_malloc( ref_plot3d->patch[n], 1, REF_PATCH_STRUCT );
      ref_patch = ref_plot3d->patch[n];
      RES( 1, fread( &(ref_patch->idim), sizeof(REF_INT), 1,file),
	   "id" );
      SWAP_INT(ref_patch->idim);
      RES( 1, fread( &(ref_patch->jdim), sizeof(REF_INT), 1,file),
	   "jd" );
      SWAP_INT(ref_patch->jdim);
      RES( 1, fread( &(ref_patch->kdim), sizeof(REF_INT), 1,file),
	   "kd" );
      SWAP_INT(ref_patch->kdim);
    }

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "end rec" );
  SWAP_INT(fortran_record_size);
  REIS( ref_plot3d_ngrid(ref_plot3d)*3*4, fortran_record_size, 
	"dims end record size" );

  for (n=0;n<ref_plot3d_ngrid(ref_plot3d);n++)
    {
      ref_patch = ref_plot3d->patch[n];
      ref_malloc( ref_patch->xyz, 3*ref_patch->idim*ref_patch->jdim, REF_DBL );

      RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "sxyz" );
      SWAP_INT(fortran_record_size);
      REIS( 8 * 3 * 
	    ref_patch->idim *
	    ref_patch->jdim *
	    ref_patch->kdim, fortran_record_size, "dims end record size" );

      for (ixyz=0;ixyz<3;ixyz++)
	{
	  for (j=0;j<ref_patch->jdim;j++)
	    for (i=0;i<ref_patch->idim;i++)
	      {
		RES( 1, fread( &ref_patch_xyz(ref_patch,ixyz,i,j), 
			       sizeof(REF_DBL), 1, file ), "xyz" );
		SWAP_DBL(ref_patch_xyz(ref_patch,ixyz,i,j));
	      }
	  for (k=1;k<ref_patch->kdim;k++)
	    for (j=0;j<ref_patch->jdim;j++)
	      for (i=0;i<ref_patch->idim;i++)
		RES( 1, fread( &dummy, sizeof(REF_DBL), 1, file ), "dum xyz" );
	}

      RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "sxyz" );
      SWAP_INT(fortran_record_size);
      REIS( 8 * 3 * 
	    ref_patch->idim *
	    ref_patch->jdim *
	    ref_patch->kdim, fortran_record_size, "dims end record size" );

    }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_plot3d_free( REF_PLOT3D ref_plot3d )
{
  REF_INT n;

  for (n=0;n<ref_plot3d_ngrid(ref_plot3d);n++)
    RSS( ref_patch_free( ref_plot3d->patch[n] ), "patch free" );
  
  ref_free( ref_plot3d->patch );
  ref_free( ref_plot3d );

  return REF_SUCCESS;
}

REF_STATUS ref_patch_free( REF_PATCH ref_patch )
{
  ref_free( ref_patch->xyz );
  ref_free( ref_patch );

  return REF_SUCCESS;
}

REF_STATUS ref_plot3d_tec( REF_PLOT3D ref_plot3d, char *filename )
{
  FILE *file;
  REF_INT i, j, n;
  REF_PATCH ref_patch;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"tecplot refine geometry file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  for (n=0;n<ref_plot3d_ngrid(ref_plot3d);n++)
    {
      ref_patch = ref_plot3d->patch[n];

      fprintf(file,
	      "zone t=patch%d, I=%d, J=%d, datapacking=%s\n",
	      n+1, ref_patch->idim, ref_patch->jdim, "point" );

      for (j=0;j<ref_patch->jdim;j++)
	for (i=0;i<ref_patch->idim;i++)
	  fprintf(file, " %.16e %.16e %.16e\n",
		  ref_patch_xyz(ref_patch,0,i,j),
		  ref_patch_xyz(ref_patch,1,i,j),
		  ref_patch_xyz(ref_patch,2,i,j) ) ;
    }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_plot3d_mate( REF_PLOT3D ref_plot3d, REF_GRID ref_grid )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT cell;
  REF_INT node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz[3], uv[2];
  REF_INT n;
  REF_PATCH ref_patch;

  each_ref_cell_valid_cell_with_nodes(ref_cell,cell,nodes)
    {
      switch ( nodes[3] )
	{
	case  9: n = 0; break;
	case 10: n = 1; break;
	case 11: n = 2; break;
	case 12: n = 3; break;
	case 13: n = 4; break;
	default: continue; break;    
	}
      ref_patch = ref_plot3d->patch[n];
      for ( node=0; node<3; node++ )
	{
	  xyz[0] = ref_node_xyz(ref_node, 0, nodes[node] );
	  xyz[1] = ref_node_xyz(ref_node, 1, nodes[node] );
	  xyz[2] = ref_node_xyz(ref_node, 2, nodes[node] );
	  RSS( ref_patch_locate( ref_patch, xyz, uv), "locate");
	}
    }

  return REF_SUCCESS;
}


REF_STATUS ref_patch_locate( REF_PATCH ref_patch, REF_DBL *xyz, REF_DBL *uv )
{
  uv[0]=0.0;
  uv[1]=0.0;

  RSS( ref_patch_xyz_at( ref_patch, uv, xyz ), "at" );

  return REF_SUCCESS;
}

REF_STATUS ref_patch_xyz_at( REF_PATCH ref_patch, REF_DBL *uv, REF_DBL *xyz )
{
  REF_INT i, j;

  i = (REF_INT)(uv[0]*(REF_DBL)(ref_patch->idim));
  j = (REF_INT)(uv[1]*(REF_DBL)(ref_patch->jdim));

  /*
  printf("i %d\n",i);
  printf("j %d\n",j);
  */

  xyz[0] = ref_patch_xyz(ref_patch,0,i,j);
  xyz[1] = ref_patch_xyz(ref_patch,1,i,j);
  xyz[2] = ref_patch_xyz(ref_patch,2,i,j);

  return REF_SUCCESS;
}
