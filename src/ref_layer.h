
#ifndef REF_LAYER_H
#define REF_LAYER_H

#include "ref_defs.h"
#include "ref_grid.h"
#include "ref_dict.h"

BEGIN_C_DECLORATION
typedef struct REF_LAYER_STRUCT REF_LAYER_STRUCT;
typedef REF_LAYER_STRUCT * REF_LAYER;
END_C_DECLORATION

#include "ref_list.h"
#include "ref_node.h"
#include "ref_cell.h"

BEGIN_C_DECLORATION

struct REF_LAYER_STRUCT {
  REF_LIST ref_list;
  REF_NODE ref_node;
  REF_CELL ref_cell;
};

REF_STATUS ref_layer_create( REF_LAYER *ref_layer );
REF_STATUS ref_layer_free( REF_LAYER ref_layer );

#define ref_layer_list(ref_layer) ((ref_layer)->ref_list)
#define ref_layer_node(ref_layer) ((ref_layer)->ref_node)
#define ref_layer_cell(ref_layer) ((ref_layer)->ref_cell)
#define ref_layer_n(ref_layer) (ref_list_n(ref_layer_list(ref_layer)))

REF_STATUS ref_layer_attach( REF_LAYER ref_layer,
			     REF_GRID ref_grid, REF_INT faceid );
REF_STATUS ref_layer_puff( REF_LAYER ref_layer, REF_GRID ref_grid );

REF_STATUS ref_layer_tec( REF_LAYER ref_layer, const char *filename );

END_C_DECLORATION

#endif /* REF_LAYER_H */
