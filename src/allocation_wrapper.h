#ifndef ALLOCATION_WRAPPER_H
#define ALLOCATION_WRAPPER_H


#include <stdlib.h>  /* Avoid breaking the usual definitions */
#include <string.h>

void *calloc_wrapper(size_t nmemb, size_t size, const char *file, unsigned int line, const char *func);
void *malloc_wrapper(size_t size, const char *file, unsigned int line, const char *func);
void *realloc_wrapper(void *ptr, size_t size, const char *file, unsigned int line, const char *func);
char *strdup_wrapper(const char *s, const char *file, unsigned int line, const char *func);

#ifdef calloc
#  undef calloc
#endif
#define calloc(n, s)  calloc_wrapper( (n), (s), __FILE__, __LINE__, __func__)

#ifdef malloc
#  undef malloc
#endif
#define malloc(s)     malloc_wrapper( (s),      __FILE__, __LINE__, __func__)

#ifdef realloc
#  undef realloc
#endif
#define realloc(p, s) realloc_wrapper((p), (s), __FILE__, __LINE__, __func__)

#ifdef strdup
#  undef strdup
#endif
#define strdup(s)     strdup_wrapper( (s),      __FILE__, __LINE__, __func__)

#endif /* ALLOCATION_WRAPPER_H */
