
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include "grid.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridinsert.h"
#include "gridcad.h"
#include "gridshape.h"
#include "gridmove.h"
#include "gridmpi.h"
#include "gridfiller.h"
#include "CADGeom/CADGeom.h"

Grid *gridDumpBentEdgesForPX(Grid *grid, int order, char *filename)
{
  int conn, total;
  int parent;
  int midnode;
  int nodes[3];
  double xyz[3];
  double uv[2], uv0[2], uv1[2], uv2[2];
  double t, t0, t1;
  double s, b0, b1, b2;
  int face;
  int tri_points;
  FILE *file;

  gridCreateConn(grid);
  total = 0;
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    parent = gridParentGeometry(grid, nodes[0], nodes[1]);
    if (0 != parent ) {
      total++;
    }
  }

  printf("%d edges bent of %d total edges.\n",total,gridNConn(grid));
  file = fopen(filename,"w");
  fprintf(file,"%d geometry order\n",order);
  fprintf(file,"%10d edges in 1-base numbering\n",total);
  
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    parent = gridParentGeometry(grid, nodes[0], nodes[1]);
    if ( 0 < parent ) { /* bend edge in interior of face */
      fprintf(file,"%10d%10d",nodes[0]+1,nodes[1]+1);
      gridNodeUV(grid, nodes[0], parent, uv0);
      gridNodeUV(grid, nodes[1], parent, uv1);
      for (midnode=1;midnode<order;midnode++) {
	s = ((double)midnode) / ((double)order);
	uv[0] = s*uv1[0] + (1.0-s)*uv0[0];
	uv[1] = s*uv1[1] + (1.0-s)*uv0[1];
	gridEvaluateOnFace(grid, parent, uv, xyz );
	fprintf(file,"%24.15e%24.15e%24.15e",xyz[0],xyz[1],xyz[2]);
      }
      fprintf(file,"\n");
    }
    if ( 0 > parent ) { /* bend edge on CAD edge */
      fprintf(file,"%10d%10d",nodes[0]+1,nodes[1]+1);
      gridNodeT(grid, nodes[0], -parent, &t0);
      gridNodeT(grid, nodes[1], -parent, &t1);
      for (midnode=1;midnode<order;midnode++) {
	s = ((double)midnode) / ((double)order);
	t = s*t1 + (1.0-s)*t0;
	gridEvaluateOnEdge(grid, -parent, t, xyz );
	fprintf(file,"%24.15e%24.15e%24.15e",xyz[0],xyz[1],xyz[2]);
      }
      fprintf(file,"\n");
    }
  }
  gridEraseConn(grid);

  if (order>2) {
    switch (order) {
    case 3: tri_points = 1; break;
    case 4: tri_points = 3; break;
    case 5: tri_points = 6; break;
    default:
      printf("ERROR: gridDumpBentEdgesForPX: %s: %d: order %d %s\n",
	     __FILE__, __LINE__, order, "not implemented yet" );
      return NULL;
      break;
    }
    fprintf(file,"%10d face points in 1-base numbering\n",
	    tri_points*gridNFace(grid));
    for (face = 0; face < gridMaxFace(grid) ; face++) {
      if (grid == gridFace(grid, face, nodes, &parent)) {
	gridNodeUV(grid, nodes[0], parent, uv0);
	gridNodeUV(grid, nodes[1], parent, uv1);
	gridNodeUV(grid, nodes[2], parent, uv2);
	switch (order) {
	case 3: 
	  b0 = b1 = b2 = 1.0/3.0;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  break;
	case 4:
	  b0 = 0.5;
	  b1 = b2 = 0.25;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  b1 = 0.5;
	  b0 = b2 = 0.25;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  b2 = 0.5;
	  b0 = b1 = 0.25;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  break;
	case 5:
	  b0 = 0.6;
	  b1 = b2 = 0.2;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  b1 = 0.6;
	  b0 = b2 = 0.2;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  b2 = 0.6;
	  b0 = b1 = 0.2;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");

	  b0 = 0.2;
	  b1 = b2 = 0.4;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  b1 = 0.2;
	  b0 = b2 = 0.4;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  b2 = 0.2;
	  b0 = b1 = 0.4;
	  uv[0] = b0*uv0[0] + b1*uv1[0] + b2*uv2[0];
	  uv[1] = b0*uv0[1] + b1*uv1[1] + b2*uv2[1];
	  gridEvaluateOnFace(grid, parent, uv, xyz );
	  fprintf(file,"%10d%10d%10d",nodes[0]+1,nodes[1]+1,nodes[2]+1);
	  fprintf(file,"%18.15f%18.15f%18.15f%24.15e%24.15e%24.15e",
		  b0, b1, b2, xyz[0], xyz[1], xyz[2]);
	  fprintf(file,"\n");
	  break;
	default:
	  printf("ERROR: gridDumpBentEdgesForPX: %s: %d: order %d %s\n",
		 __FILE__, __LINE__, order, "not implemented yet" );
	  return NULL;
	  break;
	}
      }
    }
  }
  fclose(file);

  return grid;
}

Grid *gridHistogram( Grid *grid, char *filename ) 
{
  int nodes[4];
  int edge, edgeId;
  int face, faceId;
  int cell;
  double length, volume;

  int bins[200];
  int bin;

  char fn[1024];
  FILE *file;
  
  if (NULL == filename) {
    hist_index++;
    sprintf( fn, "histogram%04d.m", hist_index );
    printf("dump %s\n",fn);
    file = fopen( fn, "w" );
  }else{
    printf("dump %s\n",filename);
    file = fopen( filename, "w" );
  }

  for (bin=0;bin<200;bin++) bins[bin] = 0;
  gridCreateConn(grid);
  for ( edge = 0 ; edge < gridNConn(grid) ; edge++ ) {
    gridConn2Node( grid, edge, nodes );
    length = gridEdgeRatio(grid, nodes[0], nodes[1]);
    bin = (int)(log10(length)*100.0)+100;
    bin = MAX(0,bin);
    bin = MIN(bin,199);
    bins[bin]++;
  }
  gridEraseConn(grid);

  fprintf( file, "edge_length = [\n" );
  for (bin=0;bin<200;bin++) 
    fprintf( file, "%f %d\n", (((double)bin)-100.0)/100.0, bins[bin] );
  fprintf( file, "];\n" );

  for (bin=0;bin<200;bin++) bins[bin] = 0;
  for ( cell = 0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if ( grid == gridCell( grid, cell, nodes ) ) {
      volume = gridVolume(grid, nodes);
      bin = (int)(-log10(volume)*10.0);
      bin = MAX(0,bin);
      bin = MIN(bin,199);
      bins[bin]++;
    }
  }
  fprintf( file, "tet_volume = [\n" );
  for (bin=199;bin>=0;bin--) 
    fprintf( file, "%f %d\n", -(((double)bin))/10.0, bins[bin] );
  fprintf( file, "];\n" );

  fprintf( file, "[x,y]=bar(edge_length(:,1), edge_length(:,2));plot(x,y,';;');\n" );
  fprintf( file, "xlabel('log_{10}(L)')\n" );
  fprintf( file, "print -deps histogram_edge.eps\n" );
  fprintf( file, "closeplot\n" );
  fprintf( file, "[x,y]=bar(tet_volume(:,1), tet_volume(:,2));plot(x,y,';Tetrahedral Volume;');\n" );
  fprintf( file, "print -deps histogram_vol.eps\n" );
  fprintf( file, "closeplot\n" );

  fclose(file);

  return grid;
}
