#include <stdlib.h>  /* Avoid breaking the usual definitions */
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <errno.h>
#include <glog/logging.h>

#include "allocation_wrapper.h"

#ifdef calloc
#  undef calloc
#endif
void *calloc_wrapper(size_t nmemb, size_t size, const char *file, unsigned int line, const char *func) {
  void *p = calloc(nmemb, size);
  if (NULL == p) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  return p;
}

#ifdef malloc
#  undef malloc
#endif
void *malloc_wrapper(size_t size, const char *file, unsigned int line, const char *func) {
  void *p = malloc(size);
  if (NULL == p) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  return p;
}

#ifdef realloc
#  undef realloc
#endif
void *realloc_wrapper(void *ptr, size_t size, const char *file, unsigned int line, const char *func) {
  void *p = realloc(ptr, size);
  if (NULL == p) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  return p;
}

#ifdef strdup
#  undef strdup
#endif
char *strdup_wrapper(const char *s, const char *file, unsigned int line, const char *func) {
  char *p = strdup(s);
  if (NULL == p) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  return p;
}
