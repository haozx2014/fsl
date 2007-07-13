#include "flags.h"
#include "dbh.h"
#include "xfm.h"

#define TOO_LARGE_PPE 10000

int read_varian_header(char *infile1,struct dsr *header,struct xfm *transform);
int flags_from_procpar(char *infile,struct flags *proc);
int read_varian_chunks(char *infile1,unsigned long chunksize,unsigned int chunks,unsigned int offset,unsigned int step,float* buffer);
int read_procpar(char *filename,char *parameter,char *value);
int read_procpar_strings(char *filename,char *parameter,char *value);
long str_to_date(char *datestr);
float sw_readout(char *filename);
unsigned long sizeof_fidfile(char *infile1);
int flash_converted(char *filename);
int tab_converted(char *filename);
int volumes_in_varian_header(char *infile1);
