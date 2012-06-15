
#ifndef REF_DEBUG_H
#define REF_DEBUG_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_debug_print_int( i ) printf( " " #i " = %d\n", i); 
#define ref_debug_print_dbl( d ) printf( " " #d " = %e\n", i); 
#define ref_debug_print_str( s ) printf( " " #s "\n"); 

#define ref_debug_print_int_array( a, n )		\
  { REF_INT ref_debug_print_int_array_i;		\
    for( ref_debug_print_int_array_i = 0 ;		\
	 ref_debug_print_int_array_i<(n);		\
	 ref_debug_print_int_array_i++ )		\
      printf( " " #a "[%d] = %d\n",			\
	      ref_debug_print_int_array_i,		\
	      (a)[ref_debug_print_int_array_i]);	\
  }

END_C_DECLORATION

#endif /* REF_DEBUG_H */
