
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#ifndef GRIDFORTRAN_H
#define GRIDFORTRAN_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

void FC_FUNC(gridapiversion,GRIDAPIVERSION)( int *refine_api_version );

void FC_FUNC(gridcreate,GRIDCREATE)( int *partId, int *nnode, double *x, double *y, double *z );
void FC_FUNC(gridfree,GRIDFREE)( void );

void FC_FUNC(gridinsertcells,GRIDINSERTCELLS)( int *nodes_per_cell, int *ncell, int *c2n );
void FC_FUNC(gridinsertbc,GRIDINSERTBC)( int *faceId, int *nodes_per_face, int *nface, int *f2n );
void FC_FUNC(gridsetmap,GRIDSETMAP)( int *nnode, double* map );
void FC_FUNC(gridsetimesh,GRIDSETIMESH)( int *nnode, int* imesh );
void FC_FUNC(gridsetnodelocal2global,GRIDSETNODELOCAL2GLOBAL)( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global );
void FC_FUNC(gridsetnodepart,GRIDSETNODEPART)( int *nnode, int *part );
void FC_FUNC(gridfreezenode,GRIDFREEZENODE)( int *node );
void FC_FUNC(gridparallelloadcapri,GRIDPARALLELLOADCAPRI)( char *url, char *modeler, char *capriProject,
                             int *status );
void FC_FUNC(gridprojectallfaces,GRIDPROJECTALLFACES)( void );
void FC_FUNC(gridtestcadparameters,GRIDTESTCADPARAMETERS)( void );
void FC_FUNC(gridminar,GRIDMINAR)( double *aspectratio );
void FC_FUNC(gridwritetecplotsurfacezone,GRIDWRITETECPLOTSURFACEZONE)( void );
void FC_FUNC(gridexportfast,GRIDEXPORTFAST)( void );

void FC_FUNC(gridsetcostconstraint,GRIDSETCOSTCONSTRAINT)( int *cost_constraint );
void FC_FUNC(gridconstrainsurfacenode,GRIDCONSTRAINSURFACENODE)( void );

void FC_FUNC(gridparallelswap,GRIDPARALLELSWAP)( int *processor, double *ARlimit );
void FC_FUNC(gridparallelsmooth,GRIDPARALLELSMOOTH)( int *processor,
			  double *optimizationLimit, double *laplacianLimit,
                          int *geometryAllowed );
void FC_FUNC(gridparallelrelaxneg,GRIDPARALLELRELAXNEG)( int *processor, int *geometryAllowed );
void FC_FUNC(gridparallelrelaxsurf,GRIDPARALLELRELAXSURF)( int *processor );
void FC_FUNC(gridparalleladapt,GRIDPARALLELADAPT)( int *processor, 
			 double *minLength, double *maxLength );
void FC_FUNC(gridparallelpreproject,GRIDPARALLELPREPROJECT)( int *processor );

void FC_FUNC(queuedumpsize,QUEUEDUMPSIZE)( int *nInt, int *nDouble );
void FC_FUNC(queuedump,QUEUEDUMP)( int *nInt, int *nDouble, int *ints, double *doubles );
void FC_FUNC(gridapplyqueue,GRIDAPPLYQUEUE)( int *nInt, int *nDouble, int *ints, double *doubles );

void FC_FUNC(gridglobalnnode,GRIDGLOBALNNODE)( int *nnodeg );
void FC_FUNC(gridglobalshift,GRIDGLOBALSHIFT)( int *oldnnodeg, int *newnnodeg, int *nodeoffset );
void FC_FUNC(gridrenumberglobalnodes,GRIDRENUMBERGLOBALNODES)( int *nnode, int *new2old );

void FC_FUNC(gridnunusednodeglobal,GRIDNUNUSEDNODEGLOBAL)( int *nunused );
void FC_FUNC(gridgetunusednodeglobal,GRIDGETUNUSEDNODEGLOBAL)( int *nunused, int *unused );
void FC_FUNC(gridjoinunusednodeglobal,GRIDJOINUNUSEDNODEGLOBAL)( int *nunused, int *unused );
void FC_FUNC(gridcopyunusednodeglobal,GRIDCOPYUNUSEDNODEGLOBAL)( int *nunused, int *unused );
void FC_FUNC(grideliminateunusednodeglobal,GRIDELIMINATEUNUSEDNODEGLOBAL)( void );

void FC_FUNC(gridsortfun3d,GRIDSORTFUN3D)( int *nnodes0, int *nnodes01, int *nnodesg );

void FC_FUNC(gridgetnodes,GRIDGETNODES)( int *nnode, int *l2g, double *x, double *y, double *z);
void FC_FUNC(gridgetmap,GRIDGETMAP)( int *nnode, double* map );
void FC_FUNC(gridgetfreezestate,GRIDGETFREEZESTATE)( int *node, int *state );
void FC_FUNC(gridgetimesh,GRIDGETIMESH)( int *nnode, int *imesh);
void FC_FUNC(gridgetncell,GRIDGETNCELL)( int *nodes_per_cell, int *ncell );
void FC_FUNC(gridgetcell,GRIDGETCELL)( int *nodes_per_cell, int *cell, int *nodes );
void FC_FUNC(gridgetbcsize,GRIDGETBCSIZE)( int *ibound, int *nodes_per_face, int *nface );
void FC_FUNC(gridgetbc,GRIDGETBC)( int *ibound, int *nodes_per_face, int *face, int *f2n );

void FC_FUNC(gridsetnaux,GRIDSETNAUX)( int *naux );
void FC_FUNC(gridgetnaux,GRIDGETNAUX)( int *naux );
void FC_FUNC(gridsetauxvector,GRIDSETAUXVECTOR)( int *nnode, int *offset, double *x );
void FC_FUNC(gridsetauxmatrix,GRIDSETAUXMATRIX)( int *ndim, int *nnode, int *offset, double *x );
void FC_FUNC(gridsetauxmatrix3,GRIDSETAUXMATRIX3)( int *ndim, int *nnode, int *offset, double *x );
void FC_FUNC(gridgetauxvector,GRIDGETAUXVECTOR)( int *nnode, int *offset, double *x );
void FC_FUNC(gridgetauxmatrix,GRIDGETAUXMATRIX)( int *ndim, int *nnode, int *offset, double *x );
void FC_FUNC(gridgetauxmatrix3,GRIDGETAUXMATRIX3)( int *ndim, int *nnode, int *offset, double *x );

void FC_FUNC(gridghostcount,GRIDGHOSTCOUNT)( int *nproc, int *count );

void FC_FUNC(gridloadghostnodes,GRIDLOADGHOSTNODES)( int *nproc, int *clientindex,
			 int *clientsize, int *localnode, int *globalnode );
void FC_FUNC(gridloadglobalnodedata,GRIDLOADGLOBALNODEDATA)( int *ndim, int *nnode, int *nodes, double *data );
void FC_FUNC(gridloadlocalnodes,GRIDLOADLOCALNODES)( int *nnode, int *global, int *local );
void FC_FUNC(gridsetlocalnodedata,GRIDSETLOCALNODEDATA)( int *ndim, int *nnode, int *nodes, double *data );

void FC_FUNC(gridmovesetprojectiondisp,GRIDMOVESETPROJECTIONDISP)( void );
void FC_FUNC(gridmoverelaxstartup,GRIDMOVERELAXSTARTUP)( int *relaxationScheme );
void FC_FUNC(gridmoverelaxstartstep,GRIDMOVERELAXSTARTSTEP)( double *position); 
void FC_FUNC(gridmoverelaxsubiter,GRIDMOVERELAXSUBITER)( double *residual);
void FC_FUNC(gridmoverelaxshutdown,GRIDMOVERELAXSHUTDOWN)( void );
void FC_FUNC(gridmoveapplydisplacements,GRIDMOVEAPPLYDISPLACEMENTS)( void );

void FC_FUNC(gridmovedataleadingdim,GRIDMOVEDATALEADINGDIM)( int *ndim );
void FC_FUNC(gridmoveinitializempitest,GRIDMOVEINITIALIZEMPITEST)( void );
void FC_FUNC(gridmovecompletempitest,GRIDMOVECOMPLETEMPITEST)( void );
void FC_FUNC(gridmoveloadlocalnodedata,GRIDMOVELOADLOCALNODEDATA)( int *ndim, int *nnode, 
				 int *nodes, double *data );
void FC_FUNC(gridmovesetlocalnodedata,GRIDMOVESETLOCALNODEDATA)( int *ndim, int *nnode, 
				int *nodes, double *data );

void FC_FUNC(gridmovefree,GRIDMOVEFREE)( void );

void FC_FUNC(gridgeomsize,GRIDGEOMSIZE)( int *nGeomNode, int *nGeomEdge, int *nGeomFace );
void FC_FUNC(gridlocalboundnode,GRIDLOCALBOUNDNODE)( int *nBoundNode );
void FC_FUNC(gridgeomedgeendpoints,GRIDGEOMEDGEENDPOINTS)( int *edgeId, int *endPoints );
void FC_FUNC(gridmaxedge,GRIDMAXEDGE)( int *maxedge );
void FC_FUNC(gridedge,GRIDEDGE)( int *edge, int *edgeId, 
		int *globalnodes, int *nodeparts, 
		double *t, double *xyz);
void FC_FUNC(gridupdateedgegrid,GRIDUPDATEEDGEGRID)(int *edgeId, int *nCurveNode, double *xyz, double *t);
void FC_FUNC(gridmaxface,GRIDMAXFACE)( int *maxface );
void FC_FUNC(gridface,GRIDFACE)( int *face, int *faceId, 
		int *globalnodes, int *nodeparts, 
		double *uv, double *xyz);
void FC_FUNC(gridfaceedgecount,GRIDFACEEDGECOUNT)( int *faceId, int *faceEdgeCount );
void FC_FUNC(gridfaceedgel2g,GRIDFACEEDGEL2G)( int *faceId, int *faceEdgeCount, int *local2global );

void FC_FUNC(gridupdategeometryface,GRIDUPDATEGEOMETRYFACE)( int *faceId, int *nnode, double *xyz, double *uv,
			      int *nface, int *f2n );

void FC_FUNC(gridcreateshellfromfaces,GRIDCREATESHELLFROMFACES)( void );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
