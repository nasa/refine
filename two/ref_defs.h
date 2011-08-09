

#ifndef REF_HEADER_H
#define REF_HEADER_H

#ifdef __cplusplus
#  define BEGIN_C_DECLORATION extern "C" {
#  define END_C_DECLORATION }
#else
#  define BEGIN_C_DECLORATION
#  define END_C_DECLORATION
#endif

BEGIN_C_DECLORATION

#define EMPTY (-1)

#if !defined(ABS)
#define ABS(a)   ((a)>0?(a):-(a))
#endif

#if !defined(MIN)
#define MIN(a,b) ((a)<(b)?(a):(b)) 
#endif

#if !defined(MAX)
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

typedef int REF_STATUS;

#define REF_SUCCESS       (0)
#define REF_FAILURE       (1)

#define TSS(fcn,msg)							\
  {									\
    REF_STATUS code;							\
    code = (fcn);							\
    if (REF_SUCCESS != code){						\
      printf("%s: %d: %s: %d %s\n",__FILE__,__LINE__,__func__,code,(msg)); \
      return code;							\
    }									\
  }

#define SUPRESS_UNUSED_COMPILER_WARNING(ptr)                    \
  if (NULL == &(ptr)) printf("unused macro failed\n");

END_C_DECLORATION

#endif /* REF_HEADER_H */

