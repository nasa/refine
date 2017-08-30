
#ifndef REF_ENDIAN_H
#define REF_ENDIAN_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define SWAP_INT(x) { \
    int y; \
    char *xp = (char *)&(x); \
    char *yp = (char *)&(y); \
    *(yp+3) = *(xp+0); \
    *(yp+2) = *(xp+1); \
    *(yp+1) = *(xp+2); \
    *(yp+0) = *(xp+3); \
    (x) = y; \
  }

#define SWAP_LONG(x) { \
    double y; \
    char *xp = (char *)&(x); \
    char *yp = (char *)&(y); \
    *(yp+7) = *(xp+0); \
    *(yp+6) = *(xp+1); \
    *(yp+5) = *(xp+2); \
    *(yp+4) = *(xp+3); \
    *(yp+3) = *(xp+4); \
    *(yp+2) = *(xp+5); \
    *(yp+1) = *(xp+6); \
    *(yp+0) = *(xp+7); \
    (x) = y; \
  }

#define SWAP_DBL(x) { \
    double y; \
    char *xp = (char *)&(x); \
    char *yp = (char *)&(y); \
    *(yp+7) = *(xp+0); \
    *(yp+6) = *(xp+1); \
    *(yp+5) = *(xp+2); \
    *(yp+4) = *(xp+3); \
    *(yp+3) = *(xp+4); \
    *(yp+2) = *(xp+5); \
    *(yp+1) = *(xp+6); \
    *(yp+0) = *(xp+7); \
    (x) = y; \
  }

END_C_DECLORATION

#endif /* REF_ENDIAN_H */
