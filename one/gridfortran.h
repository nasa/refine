
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

#ifndef GRIDFORTRAN_H
#define GRIDFORTRAN_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

void REF_FORT(gridapiversion,GRIDAPIVERSION)( int *refine_api_version );

void REF_FORT(gridcreate,GRIDCREATE)( int *partId, int *nnode, double *x, double *y, double *z );
void REF_FORT(gridfree,GRIDFREE)( void );

void REF_FORT(gridinsertcells,GRIDINSERTCELLS)( int *nodes_per_cell, int *ncell, int *c2n );
void REF_FORT(gridinsertbc,GRIDINSERTBC)( int *faceId, int *nodes_per_face, int *nface, int *f2n );
void REF_FORT(gridsetmap,GRIDSETMAP)( int *nnode, double* map );
void REF_FORT(gridsetimesh,GRIDSETIMESH)( int *nnode, int* imesh );
void REF_FORT(gridsetnodelocal2global,GRIDSETNODELOCAL2GLOBAL)( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global );
void REF_FORT(gridsetnodepart,GRIDSETNODEPART)( int *nnode, int *part );
void REF_FORT(gridfreezenode,GRIDFREEZENODE)( int *node );
void REF_FORT(gridparallelloadcapri,GRIDPARALLELLOADCAPRI)( char *url, char *modeler, char *capriProject,
                             int *status );
void REF_FORT(gridprojectallfaces,GRIDPROJECTALLFACES)( void );
void REF_FORT(gridtestcadparameters,GRIDTESTCADPARAMETERS)( void );
void REF_FORT(gridminar,GRIDMINAR)( double *aspectratio );
void REF_FORT(gridwritetecplotsurfacezone,GRIDWRITETECPLOTSURFACEZONE)( void );
void REF_FORT(gridexportfast,GRIDEXPORTFAST)( void );

void REF_FORT(gridsetcostconstraint,GRIDSETCOSTCONSTRAINT)( int *cost_constraint );
void REF_FORT(gridconstrainsurfacenode,GRIDCONSTRAINSURFACENODE)( void );

void REF_FORT(gridparallelswap,GRIDPARALLELSWAP)( int *processor, double *ARlimit );
void REF_FORT(gridparallelsmooth,GRIDPARALLELSMOOTH)( int *processor,
			  double *optimizationLimit, double *laplacianLimit,
                          int *geometryAllowed );
void REF_FORT(gridparallelrelaxneg,GRIDPARALLELRELAXNEG)( int *processor, int *geometryAllowed );
void REF_FORT(gridparallelrelaxsurf,GRIDPARALLELRELAXSURF)( int *processor );
void REF_FORT(gridparalleladapt,GRIDPARALLELADAPT)( int *processor, 
			 double *minLength, double *maxLength );
void REF_FORT(gridparallelpreproject,GRIDPARALLELPREPROJECT)( int *processor );

void REF_FORT(queuedumpsize,QUEUEDUMPSIZE)( int *nInt, int *nDouble );
void REF_FORT(queuedump,QUEUEDUMP)( int *nInt, int *nDouble, int *ints, double *doubles );
void REF_FORT(gridapplyqueue,GRIDAPPLYQUEUE)( int *nInt, int *nDouble, int *ints, double *doubles );

void REF_FORT(gridglobalnnode,GRIDGLOBALNNODE)( int *nnodeg );
void REF_FORT(gridglobalshift,GRIDGLOBALSHIFT)( int *oldnnodeg, int *newnnodeg, int *nodeoffset );
void REF_FORT(gridrenumberglobalnodes,GRIDRENUMBERGLOBALNODES)( int *nnode, int *new2old );

void REF_FORT(gridnunusednodeglobal,GRIDNUNUSEDNODEGLOBAL)( int *nunused );
void REF_FORT(gridgetunusednodeglobal,GRIDGETUNUSEDNODEGLOBAL)( int *nunused, int *unused );
void REF_FORT(gridjoinunusednodeglobal,GRIDJOINUNUSEDNODEGLOBAL)( int *nunused, int *unused );
void REF_FORT(gridcopyunusednodeglobal,GRIDCOPYUNUSEDNODEGLOBAL)( int *nunused, int *unused );
void REF_FORT(grideliminateunusednodeglobal,GRIDELIMINATEUNUSEDNODEGLOBAL)( void );

void REF_FORT(gridsortfun3d,GRIDSORTFUN3D)( int *nnodes0, int *nnodes01, int *nnodesg );

void REF_FORT(gridgetnodes,GRIDGETNODES)( int *nnode, int *l2g, double *x, double *y, double *z);
void REF_FORT(gridgetmap,GRIDGETMAP)( int *nnode, double* map );
void REF_FORT(gridgetfreezestate,GRIDGETFREEZESTATE)( int *node, int *state );
void REF_FORT(gridgetimesh,GRIDGETIMESH)( int *nnode, int *imesh);
void REF_FORT(gridgetncell,GRIDGETNCELL)( int *nodes_per_cell, int *ncell );
void REF_FORT(gridgetcell,GRIDGETCELL)( int *nodes_per_cell, int *cell, int *nodes );
void REF_FORT(gridgetbcsize,GRIDGETBCSIZE)( int *ibound, int *nodes_per_face, int *nface );
void REF_FORT(gridgetbc,GRIDGETBC)( int *ibound, int *nodes_per_face, int *face, int *f2n );

void REF_FORT(gridsetnaux,GRIDSETNAUX)( int *naux );
void REF_FORT(gridgetnaux,GRIDGETNAUX)( int *naux );
void REF_FORT(gridsetauxvector,GRIDSETAUXVECTOR)( int *nnode, int *offset, double *x );
void REF_FORT(gridsetauxmatrix,GRIDSETAUXMATRIX)( int *ndim, int *nnode, int *offset, double *x );
void REF_FORT(gridsetauxmatrix3,GRIDSETAUXMATRIX3)( int *ndim, int *nnode, int *offset, double *x );
void REF_FORT(gridgetauxvector,GRIDGETAUXVECTOR)( int *nnode, int *offset, double *x );
void REF_FORT(gridgetauxmatrix,GRIDGETAUXMATRIX)( int *ndim, int *nnode, int *offset, double *x );
void REF_FORT(gridgetauxmatrix3,GRIDGETAUXMATRIX3)( int *ndim, int *nnode, int *offset, double *x );

void REF_FORT(gridghostcount,GRIDGHOSTCOUNT)( int *nproc, int *count );

void REF_FORT(gridloadghostnodes,GRIDLOADGHOSTNODES)( int *nproc, int *clientindex,
			 int *clientsize, int *localnode, int *globalnode );
void REF_FORT(gridloadglobalnodedata,GRIDLOADGLOBALNODEDATA)( int *ndim, int *nnode, int *nodes, double *data );
void REF_FORT(gridloadlocalnodes,GRIDLOADLOCALNODES)( int *nnode, int *global, int *local );
void REF_FORT(gridsetlocalnodedata,GRIDSETLOCALNODEDATA)( int *ndim, int *nnode, int *nodes, double *data );

void REF_FORT(gridmovesetprojectiondisp,GRIDMOVESETPROJECTIONDISP)( void );
void REF_FORT(gridmoverelaxstartup,GRIDMOVERELAXSTARTUP)( int *relaxationScheme );
void REF_FORT(gridmoverelaxstartstep,GRIDMOVERELAXSTARTSTEP)( double *position); 
void REF_FORT(gridmoverelaxsubiter,GRIDMOVERELAXSUBITER)( double *residual);
void REF_FORT(gridmoverelaxshutdown,GRIDMOVERELAXSHUTDOWN)( void );
void REF_FORT(gridmoveapplydisplacements,GRIDMOVEAPPLYDISPLACEMENTS)( void );

void REF_FORT(gridmovedataleadingdim,GRIDMOVEDATALEADINGDIM)( int *ndim );
void REF_FORT(gridmoveinitializempitest,GRIDMOVEINITIALIZEMPITEST)( void );
void REF_FORT(gridmovecompletempitest,GRIDMOVECOMPLETEMPITEST)( void );
void REF_FORT(gridmoveloadlocalnodedata,GRIDMOVELOADLOCALNODEDATA)( int *ndim, int *nnode, 
				 int *nodes, double *data );
void REF_FORT(gridmovesetlocalnodedata,GRIDMOVESETLOCALNODEDATA)( int *ndim, int *nnode, 
				int *nodes, double *data );

void REF_FORT(gridmovefree,GRIDMOVEFREE)( void );

void REF_FORT(gridgeomsize,GRIDGEOMSIZE)( int *nGeomNode, int *nGeomEdge, int *nGeomFace );
void REF_FORT(gridlocalboundnode,GRIDLOCALBOUNDNODE)( int *nBoundNode );
void REF_FORT(gridgeomedgeendpoints,GRIDGEOMEDGEENDPOINTS)( int *edgeId, int *endPoints );
void REF_FORT(gridmaxedge,GRIDMAXEDGE)( int *maxedge );
void REF_FORT(gridedge,GRIDEDGE)( int *edge, int *edgeId, 
		int *globalnodes, int *nodeparts, 
		double *t, double *xyz);
void REF_FORT(gridupdateedgegrid,GRIDUPDATEEDGEGRID)(int *edgeId, int *nCurveNode, double *xyz, double *t);
void REF_FORT(gridmaxface,GRIDMAXFACE)( int *maxface );
void REF_FORT(gridface,GRIDFACE)( int *face, int *faceId, 
		int *globalnodes, int *nodeparts, 
		double *uv, double *xyz);
void REF_FORT(gridfaceedgecount,GRIDFACEEDGECOUNT)( int *faceId, int *faceEdgeCount );
void REF_FORT(gridfaceedgel2g,GRIDFACEEDGEL2G)( int *faceId, int *faceEdgeCount, int *local2global );

void REF_FORT(gridupdategeometryface,GRIDUPDATEGEOMETRYFACE)( int *faceId, int *nnode, double *xyz, double *uv,
			      int *nface, int *f2n );

void REF_FORT(gridcreateshellfromfaces,GRIDCREATESHELLFROMFACES)( void );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
