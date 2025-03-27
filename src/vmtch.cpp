#include "vmtch.h"
/*
 * To find out if the stat results from a single file correspond to a file we
 * have already seen, we need to compare both the device and the inode
 */
struct dev_and_inode
{
	dev_t dev;
	ino_t ino;
};

int64_t offset = 0;
int curr_crawl_depth = 0;
ino_t crawl_inodes[PATH_MAX];

// remember all inodes (for files with inode count > 1) to find duplicates
void *seen_inodes = NULL;
dev_t orig_device = 0;
int orig_device_inited = 0;

static void fatal(const char *fmt, ...)
{
	va_list ap;
	char buf[4096];
	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);
	fprintf(stderr, "FATAL: %s\n", buf);
	exit(1);
}

static void warning(const char *fmt, ...) {
	va_list ap;
	char buf[4096];
	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);
	fprintf(stderr, "WARNING: %s\n", buf);
}

int aligned_p(void *p)
{
	return 0 == ((long)p & (sysconf(_SC_PAGESIZE) - 1));
}

void increment_nofile_rlimit()
{
	struct rlimit r;
	if (getrlimit(RLIMIT_NOFILE, &r))
		fatal("increment_nofile_rlimit: getrlimit (%s)", strerror(errno));
	r.rlim_cur = r.rlim_max + 1;
	r.rlim_max = r.rlim_max + 1;
	if (setrlimit(RLIMIT_NOFILE, &r))
	{
		if (errno == EPERM)
		{
			if (getuid() == 0 || geteuid() == 0) fatal("system open file limit reached");
			fatal("open file limit reached and unable to increase limit. retry as root");
		}
		fatal("increment_nofile_rlimit: setrlimit (%s)", strerror(errno));
	}
}

static void vmtouch_core(char *path)
{
	int fd = -1;
	void *mem = 0;
	struct stat sb;
	int64_t len_of_file=0;
	int64_t len_of_range=0;
	int res, open_flags;
retry_open:
	open_flags = O_RDONLY;
#if defined(O_NOATIME)
	open_flags |= O_NOATIME;
#endif
	fd = open(path, open_flags, 0);
#if defined(O_NOATIME)
	if (fd == -1 && errno == EPERM)
	{
		open_flags &= ~O_NOATIME;
		fd = open(path, open_flags, 0);
	}
#endif
	if (fd == -1)
	{
		if (errno == ENFILE || errno == EMFILE)
		{
			increment_nofile_rlimit();
			goto retry_open;
		}
		warning("unable to open %s (%s), skipping", path, strerror(errno));
		goto bail;
	}
	res = fstat(fd, &sb);
	if (res)
	{
		warning("unable to fstat %s (%s), skipping", path, strerror(errno));
		goto bail;
	}
	if (S_ISBLK(sb.st_mode))
	{
#if defined(__linux__)
		if (ioctl(fd, BLKGETSIZE64, &len_of_file))
		{
			warning("unable to ioctl %s (%s), skipping", path, strerror(errno));
			goto bail;
		}
#else
		fatal("discovering size of block devices not (yet?) supported on this platform");
#endif
	}
	else
		len_of_file = sb.st_size;
	if (len_of_file == 0) goto bail;
	if (offset >= len_of_file)
	{
		warning("file %s smaller than offset, skipping", path);
		goto bail;
	}
	else
		len_of_range = len_of_file - offset;
	mem = mmap(NULL, len_of_range, PROT_READ, MAP_SHARED, fd, offset);
	if (mem == MAP_FAILED)
	{
		warning("unable to mmap file %s (%s), skipping", path, strerror(errno));
		goto bail;
	}
	if (!aligned_p(mem)) fatal("mmap(%s) wasn't page aligned", path);
bail:
	if (mem && munmap(mem, len_of_range))
		warning("unable to munmap file %s (%s)", path, strerror(errno));
	if (fd != -1) close(fd);
}

// compare device and inode information
int compare_func(const void *p1, const void *p2)
{
	const struct dev_and_inode *kp1 = (struct dev_and_inode *)p1, \
							   *kp2 = (struct dev_and_inode *)p2;
	int cmp1 = (kp1->ino > kp2->ino) - (kp1->ino < kp2->ino);
	if (cmp1 != 0) return cmp1;
	return (kp1->dev > kp2->dev) - (kp1->dev < kp2->dev);
}

// add device and inode information to the tree of known inodes
static inline void add_object (struct stat *st)
{
	struct dev_and_inode *newp = (struct dev_and_inode *)malloc(sizeof(struct dev_and_inode));
	if (newp == NULL) fatal("malloc: out of memory");
	newp->dev = st->st_dev;
	newp->ino = st->st_ino;
	if (tsearch(newp, &seen_inodes, compare_func) == NULL)
		fatal("tsearch: out of memory");
}

// return true only if the device and inode information has not been added before
static inline int find_object(struct stat *st)
{
	struct dev_and_inode obj;
	void *res;
	obj.dev = st->st_dev;
	obj.ino = st->st_ino;
	res = (void *) tfind(&obj, &seen_inodes, compare_func);
	return res != (void *) NULL;
}

void vmtouch(char *path)
{
	struct stat sb;
	DIR *dirp;
	struct dirent *de;
	char npath[PATH_MAX];
	int i, res;
	int tp_path_len = strlen(path);
	if (path[tp_path_len-1] == '/' && tp_path_len > 1)
		path[tp_path_len - 1] = '\0'; // prevent ugly double slashes when printing path names
	res = lstat(path, &sb);

	if (res)
	{
		warning("unable to stat %s (%s)", path, strerror(errno));
		return;
	}
	else
	{
		if (S_ISLNK(sb.st_mode))
		{
			warning("not following symbolic link %s", path);
			return;
		}
		if (sb.st_nlink > 1)
		{
			/*
			 * For files with more than one link to it, ignore it if we already know
			 * inode.  Without this check files copied as hardlinks (cp -al) are
			 * counted twice (which may lead to a cache usage of more than 100% of
			 * RAM).
			 */
			if (find_object(&sb))
				// we already saw the device and inode referenced by this file
				return;
			else
				add_object(&sb);
		}

		if (S_ISDIR(sb.st_mode))
		{
			for (i=0; i<curr_crawl_depth; i++)
			{
				if (crawl_inodes[i] == sb.st_ino)
				{
					warning("symbolic link loop detected: %s", path);
					return;
				}
			}
			if (curr_crawl_depth == PATH_MAX)
				fatal("maximum directory crawl depth reached: %s", path);
			crawl_inodes[curr_crawl_depth] = sb.st_ino;
			retry_opendir:
			dirp = opendir(path);
			if (dirp == NULL)
			{
				if (errno == ENFILE || errno == EMFILE)
				{
					increment_nofile_rlimit();
					goto retry_opendir;
				}
				warning("unable to opendir %s (%s), skipping", path, strerror(errno));
				return;
			}
			while((de = readdir(dirp)) != NULL)
			{
				if (strcmp(de->d_name, ".") == 0 || strcmp(de->d_name, "..") == 0) continue;
				if (snprintf(npath, sizeof(npath), "%s/%s", path, de->d_name) >= sizeof(npath))
				{
					warning("path too long %s", path);
					goto bail;
				}
				curr_crawl_depth++;
				vmtouch(npath);
				curr_crawl_depth--;
			}
bail:
			if (closedir(dirp))
			{
				warning("unable to closedir %s (%s)", path, strerror(errno));
				return;
			}
		}
		else if (S_ISLNK(sb.st_mode))
		{
			warning("not following symbolic link %s", path);
			return;
		}
		else if (S_ISREG(sb.st_mode) || S_ISBLK(sb.st_mode))
			vmtouch_core(path);
		else
			warning("skipping non-regular file: %s", path);
	}
}
