#include "izlib.h"

int is_gz(FILE* fp)
{
	if(!fp) return 0;
	char buf[2];
	int gzip = 0;
	if(fread(buf, 1, 2, fp) == 2){
		if(((int)buf[0] == 0x1f) && ((int)(buf[1]&0xFF) == 0x8b)) gzip = 1;
	}
	fseek(fp, 12, SEEK_SET);
	if(fread(buf, 1, 2, fp) == 2){
		if((int)buf[0] == 0x42 && (int)(buf[1]&0xFF) == 0x43) gzip = 2;
	}
	fseek(fp, 0, SEEK_SET);
	return gzip;
}

uint32_t get_posix_filetime(FILE* fp)
{
	struct stat file_stats;
	fstat(fileno(fp), &file_stats);
	return file_stats.st_mtime;
}

gzFile gzopen(const char *in, const char *mode)
{
	gzFile fp = (gzFile_t *)calloc(1, sizeof(gzFile_t));
	fp->fp = fopen(in, mode);
	if(!fp->fp)
	{
		gzclose(fp);
		return NULL;
	}
	fp->mode = strdup(mode);
	// plain file
	if(*mode == 'r')
	{
		fp->is_plain = !is_gz(fp->fp);
		if (fp->is_plain) return fp;
	}
	// gz file
	fp->gzip_header = (struct isal_gzip_header *)calloc(1, sizeof(struct isal_gzip_header));
	isal_gzip_header_init(fp->gzip_header);
	if (*mode == 'r') // read
	{
		fp->state = (struct inflate_state *)calloc(1, sizeof(struct inflate_state));
		fp->buf_in_size = BUF_SIZE;
		fp->buf_in = malloc(fp->buf_in_size * sizeof(uint8_t));
		isal_inflate_init(fp->state);
		fp->state->crc_flag = ISAL_GZIP_NO_HDR_VER;
		fp->state->next_in = fp->buf_in;
		fp->state->avail_in = fread(fp->state->next_in, 1, fp->buf_in_size, fp->fp);
		int ret = isal_read_gzip_header(fp->state, fp->gzip_header);
		if(ret != ISAL_DECOMP_OK)
		{
			gzclose(fp);
			return NULL;
		}
	}
	else if (*mode == 'w') // write
	{
		fp->gzip_header->os = UNIX; // FIXME auto parse OS
		fp->gzip_header->time = get_posix_filetime(fp->fp);
		fp->gzip_header->name = strdup(in); 
		fp->gzip_header->name_buf_len = strlen(fp->gzip_header->name) + 1;
		fp->buf_out_size = BUF_SIZE;
		fp->buf_out = (uint8_t *)calloc(fp->buf_out_size, sizeof(uint8_t));
		fp->zstream = (struct isal_zstream *)calloc(1, sizeof(struct isal_zstream));
		isal_deflate_init(fp->zstream);
		fp->zstream->avail_in = 0;
		fp->zstream->flush = NO_FLUSH;
		fp->zstream->level = COM_LVL_DEFAULT;
		fp->zstream->level_buf_size = com_lvls[fp->zstream->level];
		fp->zstream->level_buf = (uint8_t *)calloc(fp->zstream->level_buf_size, sizeof(uint8_t));
		fp->zstream->gzip_flag = IGZIP_GZIP_NO_HDR;
		fp->zstream->avail_out = fp->buf_out_size;
		fp->zstream->next_out = fp->buf_out;
		int ret = isal_write_gzip_header(fp->zstream, fp->gzip_header);
		if(ret != ISAL_DECOMP_OK)
		{
			gzclose(fp);
			return NULL;
		}
	}
	return fp;
}

gzFile gzdopen(int fd, const char *mode)
{
	char path[10];         /* identifier for error messages */
	if (fd == -1)
		return NULL;
	sprintf(path, "<fd:%d>", fd);   /* for debugging */
	gzFile fp = calloc(1, sizeof(gzFile_t));
	fp->fp = fdopen(fd, mode);
	if(!fp->fp)
	{
		gzclose(fp);
		return NULL;
	}
	fp->mode = strdup(mode);
	// plain file
	if(*mode == 'r')
	{
		fp->is_plain = !is_gz(fp->fp);
		if (fp->is_plain) return fp;
	}
	// gz file
	fp->gzip_header = calloc(1, sizeof(struct isal_gzip_header));
	isal_gzip_header_init(fp->gzip_header);
	if (*mode == 'r') // read
	{
		fp->state = calloc(1, sizeof(struct inflate_state));
		fp->buf_in_size = BUF_SIZE;
		fp->buf_in = malloc(fp->buf_in_size * sizeof(uint8_t));
		isal_inflate_init(fp->state);
		fp->state->crc_flag = ISAL_GZIP_NO_HDR_VER;
		fp->state->next_in = fp->buf_in;
		fp->state->avail_in = fread(fp->state->next_in, 1, fp->buf_in_size, fp->fp);
		int ret = isal_read_gzip_header(fp->state, fp->gzip_header);
		if(ret != ISAL_DECOMP_OK)
		{
			gzclose(fp);
			return NULL;
		}
	}
	else if (*mode == 'w') // write
	{
		fp->gzip_header->os = UNIX; // FIXME auto parse OS
		fp->gzip_header->time = get_posix_filetime(fp->fp);
		fp->gzip_header->name = strdup(path); 
		fp->gzip_header->name_buf_len = strlen(fp->gzip_header->name) + 1;
		fp->buf_out_size = BUF_SIZE;
		fp->buf_out = calloc(fp->buf_out_size, sizeof(uint8_t));
		fp->zstream = calloc(1, sizeof(struct isal_zstream));
		isal_deflate_init(fp->zstream);
		fp->zstream->avail_in = 0;
		fp->zstream->flush = NO_FLUSH;
		fp->zstream->level = COM_LVL_DEFAULT;
		fp->zstream->level_buf_size = com_lvls[fp->zstream->level];
		fp->zstream->level_buf = calloc(fp->zstream->level_buf_size, sizeof(uint8_t));
		fp->zstream->gzip_flag = IGZIP_GZIP_NO_HDR;
		fp->zstream->avail_out = fp->buf_out_size;
		fp->zstream->next_out = fp->buf_out;
		int ret = isal_write_gzip_header(fp->zstream, fp->gzip_header);
		if(ret != ISAL_DECOMP_OK)
		{
			gzclose(fp);
			return NULL;
		}
	}
	return fp;
}


const char *zError(int err)
{
	return z_errmsg[Z_NEED_DICT-(err)];
}

const char *gzerror(gzFile fp, int *errnum)
{
	/* get internal structure and check integrity */
	if (fp == NULL)
		return NULL;
	if (*fp->mode != 'r' && *fp->mode != 'w')
		return NULL;
}

int gzclose(gzFile fp)
{
	if(!fp) return -1;
	if(fp->mode) free(fp->mode);
	if(fp->zstream && fp->fp) gzwrite(fp, NULL, 0);
	if(fp->gzip_header)
	{
		if(fp->gzip_header->name) free(fp->gzip_header->name);
		free(fp->gzip_header);
	}
	if(fp->state) free(fp->state);
	if(fp->buf_in) free(fp->buf_in);
	if(fp->buf_out) free(fp->buf_out);
	if(fp->zstream){
		if(fp->zstream->level_buf) free(fp->zstream->level_buf);
		free(fp->zstream);
	}
	if(fp->fp) fclose(fp->fp);
	free(fp);
	return 0;
}

int gzread(gzFile fp, void *buf, size_t len)
{
	int buf_data_len = 0;
	if (fp->is_plain)
	{
		if(!feof(fp->fp))
			buf_data_len = fread((uint8_t *)buf, 1, len, fp->fp);
		return buf_data_len;
	}
	do // Start reading in compressed data and decompress
	{
		if (!fp->state->avail_in)
		{
			fp->state->next_in = fp->buf_in;
			fp->state->avail_in = fread(fp->state->next_in, 1, fp->buf_in_size, fp->fp);
		}
		fp->state->next_out = (uint8_t *)buf;
		fp->state->avail_out = len;
		if (isal_inflate(fp->state) != ISAL_DECOMP_OK)
			return -3;
		if ((buf_data_len = fp->state->next_out - (uint8_t *)buf))
			return buf_data_len;
	} while (fp->state->block_state != ISAL_BLOCK_FINISH // while not done
		&& (!feof(fp->fp) || !fp->state->avail_out)); // and work to do
	// Add the following to look for and decode additional concatenated files
	if (!feof(fp->fp) && !fp->state->avail_in)
	{
		fp->state->next_in = fp->buf_in;
		fp->state->avail_in = fread(fp->state->next_in, 1, fp->buf_in_size, fp->fp);
	}
	while (fp->state->avail_in && fp->state->next_in[0] == 31) // 0x1f
	{
		// Look for magic numbers for gzip header. Follows the gzread() decision
		// whether to treat as trailing junk
		if (fp->state->avail_in > 1 && fp->state->next_in[1] != 139) // 0x8b
			break;
		isal_inflate_reset(fp->state);
		fp->state->crc_flag = ISAL_GZIP; // Let isal_inflate() process extra headers
		do
		{
			if (!feof(fp->fp) && !fp->state->avail_in)
			{
				fp->state->next_in = fp->buf_in;
				fp->state->avail_in = fread(fp->state->next_in, 1, fp->buf_in_size, fp->fp);
			}
			fp->state->next_out = (uint8_t *)buf;
			fp->state->avail_out = len;
			if (isal_inflate(fp->state) != ISAL_DECOMP_OK)
				return -3;
			if((buf_data_len = fp->state->next_out - (uint8_t *)buf))
				return buf_data_len;
		} while (fp->state->block_state != ISAL_BLOCK_FINISH
				&& (!feof(fp->fp) || !fp->state->avail_out));
	}
	return buf_data_len;
}

int set_compress_level(gzFile fp, int level)
{
	if (!fp || !fp->mode || *fp->mode != 'w') return -1;
	if (level < MIN_COM_LVL || level > MAX_COM_LVL) return -1;
	if (fp->zstream->level != level)
	{
		fp->zstream->level = level;
		fp->zstream->level_buf_size = com_lvls[fp->zstream->level];
		fp->zstream->level_buf = realloc(fp->zstream->level_buf,
				fp->zstream->level_buf_size * sizeof(uint8_t));
	}
	return 0;
}

int gzwrite(gzFile fp, void *buf, size_t _len)
{
	fp->zstream->next_in = (uint8_t *)buf;
	fp->zstream->avail_in = _len;
	fp->zstream->end_of_stream = !buf;
	size_t len = 0;
	do
	{
		if(!fp->zstream->next_out)
		{
			fp->zstream->next_out = fp->buf_out;
			fp->zstream->avail_out = fp->buf_out_size;
		}
		int ret = isal_deflate(fp->zstream);
		if (ret != ISAL_DECOMP_OK) return -3;
		len += fwrite(fp->buf_out, 1, fp->zstream->next_out - fp->buf_out, fp->fp);
		fp->zstream->next_out = NULL;
	} while (!fp->zstream->avail_out);
	return len;
}
