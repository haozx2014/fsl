/*
  byteorder.c: Performs byte reordering from Sun/SGI to Intel/Alpha
  
  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "byteorder.h"
#include <cstdio>

int find_byte_order(int *reverse,int *sshort,int *slong)
{
  union u {short s;char c;} t;
  t.s = 1;
  if(t.c) *reverse=1;
  else *reverse=0;
  *sshort=sizeof(short);
  *slong=sizeof(long);

  if((sizeof(long)!=4)&&(sizeof(long)!=8))return(1);
  if((sizeof(short)!=2)&&(sizeof(short)!=4))return(1);

  return(0);
}

void convert_filehead(struct datafilehead_byte *in, struct datafilehead *fh,
		      int rev, int sshort, int slong)
{
  if(slong==4){
    fh->nblocks=four_byte_to_four_byte_long(in->nblocks,rev);
    fh->ntraces=four_byte_to_four_byte_long(in->ntraces,rev);
    fh->np=four_byte_to_four_byte_long(in->np,rev);
    fh->ebytes=four_byte_to_four_byte_long(in->ebytes,rev);
    fh->tbytes=four_byte_to_four_byte_long(in->tbytes,rev);
    fh->bbytes=four_byte_to_four_byte_long(in->bbytes,rev);
    fh->nbheaders=four_byte_to_four_byte_long(in->nbheaders,rev);
  }
  else {
    fh->nblocks=four_byte_to_eight_byte_long(in->nblocks,rev);
    fh->ntraces=four_byte_to_eight_byte_long(in->ntraces,rev);
    fh->np=four_byte_to_eight_byte_long(in->np,rev);
    fh->ebytes=four_byte_to_eight_byte_long(in->ebytes,rev);
    fh->tbytes=four_byte_to_eight_byte_long(in->tbytes,rev);
    fh->bbytes=four_byte_to_eight_byte_long(in->bbytes,rev);
    fh->nbheaders=four_byte_to_eight_byte_long(in->nbheaders,rev);
  }
  if(sshort==2){
    fh->vers_id=two_byte_to_two_byte_short(in->vers_id,rev);
    fh->status=two_byte_to_two_byte_short(in->status,rev);
  }
}
void get_dc_from_blockhead(struct datablockhead_byte bhb,float *dcr,float *dci)
{
  int reverse;  
  union u {short s;char c;} p;
  union v {float f;char c[4];} t;
  p.s = 1;
  if(p.c) reverse=1;
  else reverse=0;

  /* Can't deal with change in data type size yet */
  if(sizeof(float)!=4){
    *dcr=*dci=0.0;
    return;
  }
  /* Real */
  if(reverse){
    t.c[0]=bhb.lvl[3];
    t.c[1]=bhb.lvl[2];
    t.c[2]=bhb.lvl[1];
    t.c[3]=bhb.lvl[0];
  }
  else {
    t.c[0]=bhb.lvl[0];
    t.c[1]=bhb.lvl[1];
    t.c[2]=bhb.lvl[2];
    t.c[3]=bhb.lvl[3];
  }
  *dcr=t.f;

  /* Imaginary */
  if(reverse){
    t.c[0]=bhb.tlt[3];
    t.c[1]=bhb.tlt[2];
    t.c[2]=bhb.tlt[1];
    t.c[3]=bhb.tlt[0];
  }
  else {
    t.c[0]=bhb.tlt[0];
    t.c[1]=bhb.tlt[1];
    t.c[2]=bhb.tlt[2];
    t.c[3]=bhb.tlt[3];
  }
  *dci=t.f;
}
void two_byte_array_to_float(char *in,float *out,long num,
				 int rev, int sshort,float dcr,float dci)
{
  long l;
  float dc;

  for(l=0;l<num;l++){
  if(l%2)dc=dci;
  else dc=dcr;
    if(sshort==2){
      out[l] = (float)two_byte_to_two_byte_short(&in[l*2],rev) - dc;
    }
    else {
      out[l] = (float)two_byte_to_four_byte_short(&in[l*2],rev) - dc;
    }
  }
}
void four_byte_array_to_float(char *in,float *out,long num,
				 int rev, int slong,float dcr,float dci)
{
  long l;
  float dc;
  
  for(l=0;l<num;l++){
    if(l%2)dc=dci;
    else dc=dcr;
    if(slong==4){
      out[l] = (float)four_byte_to_four_byte_long(&in[l*4],rev) - dc;
    }
    else {
      out[l] = (float)four_byte_to_eight_byte_long(&in[l*4],rev) - dc;
    }
  }
}

/* Sun to native conversions */

long four_byte_to_four_byte_long(char *in,int rev)
{
  union u {long l;char c[4];} t;
  if(rev){
    t.c[0]=in[3];
    t.c[1]=in[2];
    t.c[2]=in[1];
    t.c[3]=in[0];
  }
  else {
    t.c[0]=in[0];
    t.c[1]=in[1];
    t.c[2]=in[2];
    t.c[3]=in[3];
  }
  return(t.l);
}
long four_byte_to_eight_byte_long(char *in,int rev)
{
  union u {long l;char c[8];} t;
  t.l=0;
  if(rev){
    if(in[0]&128){
      t.c[0]=~in[3];
      t.c[1]=~in[2];
      t.c[2]=~in[1];
      t.c[3]=~in[0];
      t.l=~t.l;
    }
    else {
      t.c[0]=in[3];
      t.c[1]=in[2];
      t.c[2]=in[1];
      t.c[3]=in[0];
    }
  }
  else {
    if(in[0]&128){
      t.c[4]=~in[0];
      t.c[5]=~in[1];
      t.c[6]=~in[2];
      t.c[7]=~in[3];
      t.l=~t.l;
    }
    else {
      t.c[4]=in[0];
      t.c[5]=in[1];
      t.c[6]=in[2];
      t.c[7]=in[3];
    }
  }
  return(t.l);
}
short two_byte_to_two_byte_short(char *in,int rev)
{
  union u {short s;char c[2];} t;
  if(rev){
    t.c[0]=in[1];
    t.c[1]=in[0];
  }
  else {
    t.c[0]=in[0];
    t.c[1]=in[1];
  }
  return(t.s);
}
short two_byte_to_four_byte_short(char *in,int rev)
{
  union u {short s;char c[4];} t;
  t.s=0;
  if(rev){
    if(t.c[0]&128){
      t.c[0]=~in[1];
      t.c[1]=~in[0];
      t.s=~t.s;
    }
    else {
      t.c[0]=in[1];
      t.c[1]=in[0];
    }
  }
  else {
    if(t.c[0]&128){
      t.c[2]=~in[0];
      t.c[3]=~in[1];
      t.s=~t.s;
    }
    else {
      t.c[2]=in[0];
      t.c[3]=in[1];
    }
  }
  return(t.s);
}

/* Native to Sun conversions */

void eight_byte_long_to_four_bytes(long in,char out[4],int rev)
{
  union u {long l;char c[8];} t;
  t.l=in;
  if(rev){
    if(t.c[0]&128){
      t.l=~t.l;
      out[0]=~t.c[3];
      out[1]=~t.c[2];
      out[2]=~t.c[1];
      out[3]=~t.c[0];
    } else {
      out[0]=t.c[3];
      out[1]=t.c[2];
      out[2]=t.c[1];
      out[3]=t.c[0];     
    }
  } else {
    /* Not checked! */
    if(t.c[4]&128){
      t.l=~t.l;
      out[0]=~t.c[4];
      out[1]=~t.c[5];
      out[2]=~t.c[6];
      out[3]=~t.c[7];
    } else {
      out[0]=t.c[4];
      out[1]=t.c[5];
      out[2]=t.c[6];
      out[3]=t.c[7];     
    }
  }
}
void eight_byte_long_to_two_bytes(long in,char out[2],int rev)
{
  union u {long l;char c[8];} t;
  t.l=in;
  if(rev){
    if(t.c[0]&128){
      t.l=~t.l;
      out[0]=~t.c[1];
      out[1]=~t.c[0];
    } else {
      out[0]=t.c[1];
      out[1]=t.c[0];
    }
  } else {
    /* Not checked! */
    if(t.c[6]&128){
      t.l=~t.l;
      out[0]=~t.c[6];
      out[1]=~t.c[7];
    } else {
      out[0]=t.c[6];
      out[1]=t.c[7];
    }
  }  
}
void four_byte_long_to_four_bytes(long in,char out[4],int rev)
{
  union u {long l;char c[4];} t;
  t.l=in;
  if(rev){
    out[0]=t.c[3];
    out[1]=t.c[2];
    out[2]=t.c[1];
    out[3]=t.c[0];
  } else {
    out[0]=t.c[0];
    out[1]=t.c[1];
    out[2]=t.c[2];
    out[3]=t.c[3];
  }
}
void four_byte_long_to_two_bytes(long in,char out[2],int rev)
{
  union u {long l;char c[4];} t;
  t.l=in;
  if(rev){
    if(t.c[1]&128){
      t.l=~t.l;
      out[0]=~t.c[1];
      out[1]=~t.c[0];
    } else {
      out[0]=t.c[1];
      out[1]=t.c[0];
    }
  } else {
    if(t.c[2]&128){
      t.l=~t.l;
      out[0]=~t.c[2];
      out[1]=~t.c[3];
    } else {
      out[0]=t.c[2];
      out[1]=t.c[3];
    }
  }
}
