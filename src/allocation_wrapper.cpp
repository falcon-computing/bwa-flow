#include <stdlib.h>  /* Avoid breaking the usual definitions */
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <glog/logging.h>

#include "allocation_wrapper.h"

#ifdef calloc
#  undef calloc
#endif
void *calloc_wrapper(size_t nmemb, size_t size, const char *file, unsigned int line, const char *func) {
  void *p = calloc(nmemb, size);
  if (NULL == p) {
    LOG(ERROR) << strerror(errno) << " due to "
               << ((errno==12) ? "out-of-memory" : "internal failure") ;
    exit(EXIT_FAILURE);
  }
  return p;
}

#ifdef malloc
#  undef malloc
#endif
void *malloc_wrapper(size_t size, const char *file, unsigned int line, const char *func) {
  void *p = malloc(size);
  if (NULL == p) {
    LOG(ERROR) << strerror(errno) << " due to "
               << ((errno==12) ? "out-of-memory" : "internal failure") ;
    exit(EXIT_FAILURE);
  }
  return p;
}

#ifdef realloc
#  undef realloc
#endif
void *realloc_wrapper(void *ptr, size_t size, const char *file, unsigned int line, const char *func) {
  void *p = realloc(ptr, size);
  if (NULL == p) {
    LOG(ERROR) << strerror(errno) << " due to "
               << ((errno==12) ? "out-of-memory" : "internal failure") ;
    exit(EXIT_FAILURE);
  }
  return p;
}

#ifdef strdup
#  undef strdup
#endif
char *strdup_wrapper(const char *s, const char *file, unsigned int line, const char *func) {
  char *p = strdup(s);
  if (NULL == p) {
    LOG(ERROR) << strerror(errno) << " due to "
               << ((errno==12) ? "out-of-memory" : "internal failure") ;
    exit(EXIT_FAILURE);
  }
  return p;
}
