

#ifndef REF_TEST_H
#define REF_TEST_H

BEGIN_C_DECLORATION

#define TSS(fcn,msg)							\
  {									\
    REF_STATUS code;							\
    code = (fcn);							\
    if (REF_SUCCESS != code){						\
      printf("%s: %d: %s: %d %s\n",__FILE__,__LINE__,__func__,code,(msg)); \
      return code;							\
    }else{								\
      printf("PASS: %s\n",(msg));					\
    }									\
  }

#define TNS(ptr,msg)							\
  {									\
    if (NULL == (ptr)){							\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg)); \
      return REF_NULL;							\
    }else{								\
      printf("PASS: %s\n",(msg));					\
    }									\
  }

#define TES(a,b,msg)							\
  {									\
    if ((a)!=(b)){							\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg)); \
      return REF_FAILURE;						\
    }else{								\
      printf("PASS: %s\n",(msg));					\
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
      printf("PASS: %s\n",(msg));					\
    }									\
  }

#define TAS(a,msg)							\
  {									\
    if (!(a)){								\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg)); \
      return REF_FAILURE;						\
    }else{								\
      printf("PASS: %s\n",(msg));					\
    }									\
  }

#define TFS(fcn,msg)							\
  {									\
    REF_STATUS code;							\
    code = (fcn);							\
    if (REF_SUCCESS == code){						\
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,(msg)); \
      return code;							\
    }else{								\
      printf("PASS: %s\n",(msg));					\
    }									\
  }

END_C_DECLORATION

#endif /* REF_TEST_H */

