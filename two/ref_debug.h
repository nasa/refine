
#ifndef REF_DEBUG_H
#define REF_DEBUG_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_print_int( i ) printf( " " #i " = %d\n", i); 
#define ref_print_dbl( d ) printf( " " #d " = %e\n", i); 
#define ref_print_str( s ) printf( " " #s "\n"); 

END_C_DECLORATION

#endif /* REF_DEBUG_H */
