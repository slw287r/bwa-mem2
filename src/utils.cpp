/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "izlib.h"
#include <errno.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <sys/resource.h>
#include <sys/time.h>
#include "utils.h"

#include "ksort.h"
#define pair64_lt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y < (b).y))
KSORT_INIT(128, pair64_t, pair64_lt)
KSORT_INIT(64,  uint64_t, ks_lt_generic)

#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

/********************
 * System utilities *
 ********************/

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
    FILE *fp = 0;
    if (strcmp(fn, "-") == 0)
        return (strstr(mode, "r"))? stdin : stdout;
    if ((fp = fopen(fn, mode)) == 0) {
        err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
    }
    return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
    if (freopen(fn, mode, fp) == 0) {
        err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
    }
    return fp;
}

gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
    gzFile fp;
    if (strcmp(fn, "-") == 0) {
        fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
        /* According to zlib.h, this is the only reason gzdopen can fail */
        if (!fp) err_fatal(func, "Out of memory");
        return fp;
    }
    if ((fp = gzopen(fn, mode)) == 0) {
        err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
    }
    return fp;
}

void err_fatal(const char *header, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[%s] ", header);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[%s] ", header);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, " Abort!\n");
    va_end(args);
    abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
    fprintf(stderr, "[%s] %s\n", func, msg);
    exit(EXIT_FAILURE);
}

void _err_fatal_simple_core(const char *func, const char *msg)
{
    fprintf(stderr, "[%s] %s Abort!\n", func, msg);
    abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t ret = fwrite(ptr, size, nmemb, stream);
    if (ret != nmemb) 
        _err_fatal_simple("fwrite", strerror(errno));
    return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t ret = fread(ptr, size, nmemb, stream);
    if (ret != nmemb)
    {
        _err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
    }
    return ret;
}

int err_fseek(FILE *stream, long offset, int whence)
{
    int ret = fseek(stream, offset, whence);
    if (0 != ret)
    {
        _err_fatal_simple("fseek", strerror(errno));
    }
    return ret;
}

long err_ftell(FILE *stream)
{
    long ret = ftell(stream);
    if (-1 == ret)
    {
        _err_fatal_simple("ftell", strerror(errno));
    }
    return ret;
}

int err_printf(const char *format, ...) 
{
    va_list arg;
    int done;
    va_start(arg, format);
    done = vfprintf(stdout, format, arg);
    int saveErrno = errno;
    va_end(arg);
    if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
    return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
    va_list arg;
    int done;
    va_start(arg, format);
    done = vfprintf(stream, format, arg);
    int saveErrno = errno;
    va_end(arg);
    if (done < 0) _err_fatal_simple("vfprintf", strerror(saveErrno));
    return done;
}

int err_fputc(int c, FILE *stream)
{
    int ret = putc(c, stream);
    if (EOF == ret)
    {
        _err_fatal_simple("fputc", strerror(errno));
    }

    return ret;
}

char* err_fgets(char *str, int size, FILE *stream)
{
    char* ret = fgets(str, size, stream );
    if (ret == NULL)
    {
        _err_fatal_simple("fgets", strerror(errno));
    }

    return ret;
}

int err_fputs(const char *s, FILE *stream)
{
    int ret = fputs(s, stream);
    if (EOF == ret)
    {
        _err_fatal_simple("fputs", strerror(errno));
    }

    return ret;
}

int err_puts(const char *s)
{
    int ret = puts(s);
    if (EOF == ret)
    {
        _err_fatal_simple("puts", strerror(errno));
    }

    return ret;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) _err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
    /* Calling fflush() ensures that all the data has made it to the
       kernel buffers, but this may not be sufficient for remote filesystems
       (e.g. NFS, lustre) as an error may still occur while the kernel
       is copying the buffered data to the file server.  To be sure of
       catching these errors, we need to call fsync() on the file
       descriptor, but only if it is a regular file.  */
    {
        struct stat sbuf;
        if (0 != fstat(fileno(stream), &sbuf))
            _err_fatal_simple("fstat", strerror(errno));
        
        if (S_ISREG(sbuf.st_mode))
        {
            if (0 != fsync(fileno(stream)))
                _err_fatal_simple("fsync", strerror(errno));
        }
    }
#endif
    return ret;
}

int err_fclose(FILE *stream) 
{
    int ret = fclose(stream);
    if (ret != 0) _err_fatal_simple("fclose", strerror(errno));
    return ret;
}

/*********
 * Timer *
 *********/

double cputime()
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

/**********
 * mmap *
 *********/
const char *get_username()
{
    uid_t uid = geteuid();
    struct passwd *pw = getpwuid(uid);
    if (pw) return pw->pw_name;
    return "*";
}

void file_size(const char *fn, int64_t *size)
{
    int fd = open(fn, O_RDONLY);
    xassert(fd > -1, "Cannot open file");
    struct stat buf;
    int s = fstat(fd, &buf);
    xassert(s > -1, "cannot stat file");
    *size = buf.st_size;
}

int64_t get_memory()
{
    int64_t mem = 0;
#ifdef __linux__
	mem = sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGESIZE);
#elif defined __APPLE__
    size_t len = sizeof(mem);
    sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
#endif
    return mem;
}

int64_t max_locked_mem()
{
	int64_t mlm = 0;
	char line[LINE_MAX];
	// max locked memory (kbytes, -l)
	char cmd[] = "ulimit -l";
	FILE *fp = popen(cmd, "r");
	fgets(line, LINE_MAX, fp);
	line[strchr(line, '\n') - line] = '\0';
	if (!strncmp(line, "unlimited", strlen("unlimited")))
		mlm = get_memory();
	else
		mlm = strtoul(line, 0, 10) * 1024;
	pclose(fp);
	return mlm;
}

void *mmap_file(const char *fn, int64_t size)
{
    int fd = open(fn, O_RDONLY);
    xassert(fd > -1, "Cannot open file");

    struct stat buf;
    int s = fstat(fd, &buf);
    xassert(s > -1, "cannot stat file");

    off_t st_size = buf.st_size;
    if (size > 0) {
        xassert(st_size >= size, "bad file size");
        st_size = size;
    }

    // mmap flags:
    // MAP_PRIVATE: copy-on-write mapping. Writes not propagated to file.
    // MAP_POPULATE: prefault page tables for mapping.  Use read-ahead. Only supported for MAP_PRIVATE
    // MAP_HUGETLB: use huge pages.  Manual says it's only supported since kernel ver. 2.6.32
    //              and requires special system configuration.
    // MAP_NORESERVE: don't reserve swap space
    // MAP_LOCKED:  Lock the pages of the mapped region into memory in the manner of mlock(2)
    //              Because we try to lock the pages in memory, this call will fail if the system
    //              doesn't have sufficient physical memory.  However, without locking, if the
    //              system can't quite fit the reference the call to mmap will succeed but aligning
    //              will take forever as parts of the reference are evicted and/or reloaded from disk.
    int map_flags = MAP_FLAGS;
    //fprintf(stderr, "* mmapping file %s (%0.1fMB)\n", fn, ((double)st_size) / (1024*1024));
    void* m = mmap(0, st_size, PROT_READ, map_flags, fd, 0);
    if (m == MAP_FAILED) {
        perror(__func__);
        err_fatal("Failed to map %s file to memory\n", fn);
    }
    //fprintf(stderr, "* File %s locked in memory\n", fn);
    close(fd);
    // MADV_WILLNEED:  Expect access in the near future
    madvise(m, st_size, MADV_WILLNEED);
    return m;
}

void unmap_file(void *map, size_t map_size)
{
    if (munmap(map, map_size) < 0)
        perror(__func__);
}
