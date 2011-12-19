

#ifndef REF_TEST_H
#define REF_TEST_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define TEST_ECH0(msg) if (0) printf("PASS: %s\n",(msg))

#define TSS(fcn,msg)							\
  {									\
    REF_STATUS ref_private_testing_code;				\
    ref_private_testing_code = (fcn);					\
    if (REF_SUCCESS != ref_private_testing_code){			\
      printf("%s: %d: %s: %d %s\n",__FILE__,__LINE__,__func__,		\
	     ref_private_testing_code,(msg));				\
      return ref_private_testing_code;					\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define TNS(ptr,msg)							\
  {									\
    if (NULL == (ptr)){							\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg));	\
      return REF_NULL;							\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define TES(a,b,msg)							\
  {									\
    if ((a)!=(b)){							\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg));	\
      return REF_FAILURE;						\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define TEIS(a,b,msg)							\
  {									\
    REF_INT ref_private_testing_ai,ref_private_testing_bi;		\
    ref_private_testing_ai = (a);					\
    ref_private_testing_bi = (b);					\
    if (ref_private_testing_ai!=ref_private_testing_bi){		\
      printf("%s: %d: %s: %s\nexpected %d was %d\n",			\
	     __FILE__,__LINE__,__func__,(msg),				\
	     ref_private_testing_ai,ref_private_testing_bi);		\
      return REF_FAILURE;						\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define TWDS(a,b,tol,msg)						\
  {									\
    REF_DBL ref_private_testing_ad,ref_private_testing_bd;		\
    REF_DBL ref_private_testing_del,ref_private_testing_allowed;	\
    ref_private_testing_ad = (a);					\
    ref_private_testing_bd = (b);					\
    ref_private_testing_del =						\
      ABS(ref_private_testing_ad-ref_private_testing_bd);		\
    ref_private_testing_allowed = (tol);				\
    if ( tol < 0.0 ) ref_private_testing_allowed = 1.0E-12;		\
    if (ref_private_testing_del>ref_private_testing_allowed){		\
      printf("%s: %d: %s: %s\nexpected %e, was %e, %e outside of %e\n",	\
	     __FILE__,__LINE__,__func__,(msg),				\
	     ref_private_testing_ad,ref_private_testing_bd,		\
	     ref_private_testing_del,ref_private_testing_allowed);	\
      return REF_FAILURE;						\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define TAS(a,msg)							\
  {									\
    if (!(a)){								\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg));	\
      return REF_FAILURE;						\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

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

#define SKIP_TEST(why)\
  printf(" *** %s *** at %s:%d\n",(why),__FILE__,__LINE__); if (REF_FALSE)

END_C_DECLORATION

#endif /* REF_TEST_H */

