#include "data.h"
#include "data_byte.h"

int find_byte_order(int *reverse,int *sshort,int *slong);
void convert_filehead(struct datafilehead_byte *in, struct datafilehead *fh,
		      int rev, int sshort, int slong);
void get_dc_from_blockhead(struct datablockhead_byte bhb,float *dcr,float *dci);
void two_byte_array_to_float(char *in,float *out,long num,
				 int rev, int sshort,float dcr,float dci);
void four_byte_array_to_float(char *in,float *out,long num,
				 int rev, int slong,float dcr,float dci);

long four_byte_to_four_byte_long(char *in,int rev);
long four_byte_to_eight_byte_long(char *in,int rev);
short two_byte_to_two_byte_short(char *in,int rev);
short two_byte_to_four_byte_short(char *in,int rev);
void eight_byte_long_to_four_bytes(long in,char out[4],int rev);
void eight_byte_long_to_two_bytes(long in,char out[2],int rev);
void four_byte_long_to_four_bytes(long in,char out[4],int rev);
void four_byte_long_to_two_bytes(long in,char out[2],int rev);
