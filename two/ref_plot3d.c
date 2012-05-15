
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

      for (j=0;j<ref_patch->jdim;j++)
	for (i=0;i<ref_patch->idim;i++)
	  for (ixyz=0;ixyz<3;ixyz++)
	    {
	      RES( 1, fread( &ref_patch_xyz(ref_patch,ixyz,i,j), 
			     sizeof(REF_DBL), 1, file ), "xyz" );
	      SWAP_DBL(ref_patch_xyz(ref_patch,ixyz,i,j));
	    }

      for (k=1;k<ref_patch->kdim;k++)
	for (j=0;j<ref_patch->jdim;j++)
	  for (i=0;i<ref_patch->idim;i++)
	    for (ixyz=0;ixyz<3;ixyz++)
	      RES( 1, fread( &dummy, sizeof(REF_DBL), 1, file ), "dummy xyz" );

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
