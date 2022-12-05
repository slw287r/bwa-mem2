#ifndef iZLIB_H
#define iZLIB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/stat.h>
#include "igzip_lib.h"

#ifndef UNIX
#define UNIX 3
#endif

#ifndef BUF_SIZE
#define BUF_SIZE (1<<22)
//#define BUF_SIZE (1e7)
#endif

#ifndef HDR_SIZE
#define HDR_SIZE (1<<16)
#endif

#ifndef MIN_COM_LVL
#define MIN_COM_LVL 0
#endif

#ifndef MAX_COM_LVL
#define MAX_COM_LVL 3
#endif

#ifndef COM_LVL_DEFAULT
#define COM_LVL_DEFAULT 2
#endif

#define Z_OK 0
#define Z_ERRNO (-1)
#define Z_NEED_DICT 2

const char * const z_errmsg[10] = {
	(const char *)"need dictionary",     /* Z_NEED_DICT       2  */
	(const char *)"stream end",          /* Z_STREAM_END      1  */
	(const char *)"",                    /* Z_OK              0  */
	(const char *)"file error",          /* Z_ERRNO         (-1) */
	(const char *)"stream error",        /* Z_STREAM_ERROR  (-2) */
	(const char *)"data error",          /* Z_DATA_ERROR    (-3) */
	(const char *)"insufficient memory", /* Z_MEM_ERROR     (-4) */
	(const char *)"buffer error",        /* Z_BUF_ERROR     (-5) */
	(const char *)"incompatible version",/* Z_VERSION_ERROR (-6) */
	(const char *)""
};

const int com_lvls[4] = {
	ISAL_DEF_LVL0_DEFAULT,
	ISAL_DEF_LVL1_DEFAULT,
	ISAL_DEF_LVL2_DEFAULT,
	ISAL_DEF_LVL3_DEFAULT
};

typedef struct
{
	FILE *fp;
	char *mode;
	int is_plain;
	struct isal_gzip_header *gzip_header;
	struct inflate_state *state;
	struct isal_zstream *zstream;
	uint8_t *buf_in;
	size_t buf_in_size;
	uint8_t *buf_out;
	size_t buf_out_size;
} gzFile_t;

typedef gzFile_t* gzFile;

#ifdef __cplusplus
extern "C" {
#endif
int is_gz(FILE* fp);
uint32_t get_posix_filetime(FILE* fp);
gzFile gzopen(const char *in, const char *mode);
gzFile gzdopen(int fd, const char *mode);
int gzread(gzFile fp, void *buf, size_t len);
int gzwrite(gzFile fp, void *buf, size_t len);
const char *zError(int err);
const char *gzerror(gzFile fp, int *errnum);
int set_compress_level(gzFile fp, int level);
int gzclose(gzFile fp);
#ifdef __cplusplus
}
#endif
#endif
