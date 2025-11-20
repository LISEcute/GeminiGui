#if !defined(liseStrcpyOS_h)
#define liseStrcpyOS_h

#include <string.h>
size_t strcpyL(char *dst, size_t dstsize, const char *src);
size_t strcatL(char *dst, size_t dstsize, const char *src);
size_t strncpyL(char *dst, size_t dstsize, const char *src, size_t count);

#if defined(linux)
#define strcpy_s(a,b,c) strcpy(a,c)
#define strcat_s(a,b,c) strcat(a,c)
#endif



#endif

