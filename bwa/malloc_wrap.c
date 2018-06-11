#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#ifdef USE_MALLOC_WRAPPERS
/* Don't wrap ourselves */
#  undef USE_MALLOC_WRAPPERS
#endif
#include "malloc_wrap.h"

void *wrap_calloc(size_t nmemb, size_t size,
				  const char *file, unsigned int line, const char *func) {
	void *p = calloc(nmemb, size);
	if (NULL == p) {
                char time_buff[25];
                struct timeval cur_tv;
                time_t cur_time;
                gettimeofday(&cur_tv, NULL);
                cur_time = (time_t)cur_tv.tv_sec;
                strftime(time_buff, 24, "%F %T", localtime(&cur_time));
                fprintf(stderr, "%s.%3d [%5d] ERROR bwa/%s:%u - %s due to out-of-memory\n",
                        time_buff, (int)(cur_tv.tv_usec/1000.0), getpid(), file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

void *wrap_malloc(size_t size,
				  const char *file, unsigned int line, const char *func) {
	void *p = malloc(size);
	if (NULL == p) {
                char time_buff[25];
                struct timeval cur_tv;
                time_t cur_time;
                gettimeofday(&cur_tv, NULL);
                cur_time = (time_t)cur_tv.tv_sec;
                strftime(time_buff, 24, "%F %T", localtime(&cur_time));
                fprintf(stderr, "%s.%3d [%5d] ERROR bwa/%s:%u - %s due to out-of-memory\n",
                        time_buff, (int)(cur_tv.tv_usec/1000.0), getpid(), file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

void *wrap_realloc(void *ptr, size_t size,
				   const char *file, unsigned int line, const char *func) {
	void *p = realloc(ptr, size);
	if (NULL == p) {
                char time_buff[25];
                struct timeval cur_tv;
                time_t cur_time;
                gettimeofday(&cur_tv, NULL);
                cur_time = (time_t)cur_tv.tv_sec;
                strftime(time_buff, 24, "%F %T", localtime(&cur_time));
                fprintf(stderr, "%s.%3d [%5d] ERROR bwa/%s:%u - %s due to out-of-memory\n",
                        time_buff, (int)(cur_tv.tv_usec/1000.0), getpid(), file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}

char *wrap_strdup(const char *s,
				  const char *file, unsigned int line, const char *func) {
	char *p = strdup(s);
	if (NULL == p) {
                char time_buff[25];
                struct timeval cur_tv;
                time_t cur_time;
                gettimeofday(&cur_tv, NULL);
                cur_time = (time_t)cur_tv.tv_sec;
                strftime(time_buff, 24, "%F %T", localtime(&cur_time));
                fprintf(stderr, "%s.%3d [%5d] ERROR bwa/%s:%u - %s due to out-of-memory\n",
                        time_buff, (int)(cur_tv.tv_usec/1000.0), getpid(), file, line, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return p;
}
