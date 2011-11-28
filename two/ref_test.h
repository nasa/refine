

#ifndef REF_TEST_H
#define REF_TEST_H

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

