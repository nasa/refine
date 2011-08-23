
#ifndef REF_CELL_H
#define REF_CELL_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_CELL_STRUCT REF_CELL_STRUCT;
typedef REF_CELL_STRUCT * REF_CELL;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_CELL_STRUCT {
  REF_INT node_per, edge_per;
  REF_INT n, max;
  REF_INT blank;
  REF_INT *c2n;
  REF_INT *c2e;
};

#define ref_cell_node_per(ref_cell) ((ref_cell)->node_per)
#define ref_cell_edge_per(ref_cell) ((ref_cell)->edge_per)

#define ref_cell_n(ref_cell) ((ref_cell)->n)
#define ref_cell_max(ref_cell) ((ref_cell)->max)
#define ref_cell_blank(ref_cell) ((ref_cell)->blank)
#define ref_cell_c2n(ref_cell,node,cell) \
  ((ref_cell)->c2n[(node)+ref_cell_node_per(ref_cell)*(cell)])
#define ref_cell_valid(ref_cell,cell) \
  ( (cell) >=0 && (cell) < ((ref_cell)->max) && \
    REF_EMPTY != (ref_cell)->c2n[ref_cell_node_per(ref_cell)*(cell)] )

#define ref_cell_for( ref_cell, cell )					\
  for ( (cell) = 0 ;							\
	(cell) < ref_cell_max(ref_cell);				\
	(cell)++ )							\
    if ( ref_cell_valid( ref_cell, cell ) )

#define ref_cell_for_with_nodes( ref_cell, cell, nodes)			\
  for ( (cell) = 0 ;							\
	(cell) < ref_cell_max(ref_cell);				\
	(cell)++ )							\
    if ( REF_SUCCESS == ref_cell_nodes( ref_cell, cell, nodes ) )

REF_STATUS ref_cell_create( REF_INT node_per, REF_CELL *ref_cell );
REF_STATUS ref_cell_free( REF_CELL ref_cell );

REF_STATUS ref_cell_inspect( REF_CELL ref_cell );

REF_STATUS ref_cell_add( REF_CELL ref_cell, REF_INT *nodes, REF_INT *cell );
REF_STATUS ref_cell_nodes( REF_CELL ref_cell, REF_INT cell, REF_INT *nodes );

END_C_DECLORATION

#endif /* REF_CELL_H */

/* http://www.simcenter.msstate.edu/docs/solidmesh/ugridconnectivity.html
in c numbering

                              inode3------5------inode2
                                 / \              . /
                                /   \          .   /
                               /     \      .     /
                              /       \  .       /
                             /        .\        /
                            2      1    4      3
                           /    .        \    /
                          /  .            \  /
                         /.                \/
                      inode0------0------inode1

                     inode3-------7----inode4
                         |    .            | \
                         |       .         |  \
                         |          .      |   \
                         |             5   |    6
                         |                .|     \
                         |                 | .    \
                         |                 |    .  \
                         2                 4       inode2
                         |                 |      . /
                         |                 |   .   /
                         |                 |.     /
                         |               . |     /
                         |            .    |    3
                         |         1       |   /
                         |      .          |  /
                         |   .             | /
                         |.                |/
                       inode0------0-----inode1

                                                  inode5
                                                  . /|
                                               .   / |
                                            .     /  |
                                         .       /   |
                                      .         /    |
                                   7           8     5
                                .             /      |
                             .               /       |
                          .                 /        |
                       inode3-----6------inode4    inode2
                         |                 |      . /
                         |                 |   .   /
                         |                 |.     /
                         |               . |     /
                         2            .    |    3
                         |         1       4   /
                         |      .          |  /
                         |   .             | /
                         |.                |/
                       inode0-----0------inode1

                               inode7-----11-----inode6
                                 /.                /|
                                / .               / |
                               /  .              /  |
                              /   .             /   |
                             9    .           10    6
                            /     7           /     |
                           /      .          /      |
                          /       .         /       |
                       inode4-8----------inode5     |
                         |      inode3.....|...5..inode2
                         |       .         |       /
                         |      .          |      /
                         |     .           |     /
                         2    1            4    3
                         |   .             |   /
                         |  .              |  /
                         | .               | /
                         |.                |/
                       inode0------0-----inode1

*/
