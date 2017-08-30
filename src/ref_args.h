
#ifndef REF_ARGS_H
#define REF_ARGS_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

REF_STATUS ref_args_inspect( REF_INT n, char **args );
REF_STATUS ref_args_find( REF_INT n, char **args,
			  const char *target, REF_INT *pos );

END_C_DECLORATION

#endif /* REF_ARGS_H */
