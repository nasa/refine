

#ifndef REF_TEST_H
#define REF_TEST_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define TEST_ECH0(msg) if (0) printf("PASS: %s\n",(msg))

#define TFS(fcn,msg)							\
  {									\
    REF_STATUS ref_private_testing_code;				\
    ref_private_testing_code = (fcn);					\
    if (REF_SUCCESS == ref_private_testing_code){			\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg));	\
      return ref_private_testing_code;					\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

END_C_DECLORATION

#endif /* REF_TEST_H */

