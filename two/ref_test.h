

#ifndef REF_TEST_H
#define REF_TEST_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define TEST_ECH0(msg) if (0) printf("PASS: %s\n",(msg))

#define TSS(fcn,msg)							\
  {									\
    REF_STATUS code;							\
    code = (fcn);							\
    if (REF_SUCCESS != code){						\
      printf("%s: %d: %s: %d %s\n",__FILE__,__LINE__,__func__,code,(msg)); \
      return code;							\
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
    REF_INT ai,bi;							\
    ai = (a);								\
    bi = (b);								\
    if (ai!=bi){							\
      printf("%s: %d: %s: %s\nexpected %d was %d\n",			\
	     __FILE__,__LINE__,__func__,(msg),ai,bi);			\
      return REF_FAILURE;						\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define TWDS(a,b,tol,msg)						\
  {									\
    REF_DBL ad,bd,del,allowed;						\
    ad = (a);								\
    bd = (b);								\
    del = ABS(ad-bd);							\
    allowed = (tol);							\
    if ( tol < 0.0 ) allowed = 1.0E-12;					\
    if (del>allowed){							\
      printf("%s: %d: %s: %s\nexpected %e, was %e, %e outside of %e\n",	\
	     __FILE__,__LINE__,__func__,(msg),ad,bd,del,allowed);	\
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
    REF_STATUS code;							\
    code = (fcn);							\
    if (REF_SUCCESS == code){						\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg));	\
      return code;							\
    }else{								\
      TEST_ECH0(msg);							\
    }									\
  }

#define SKIP_TEST(why)\
  printf(" *** %s *** at %s:%d\n",(why),__FILE__,__LINE__); if (REF_FALSE)

END_C_DECLORATION

#endif /* REF_TEST_H */

