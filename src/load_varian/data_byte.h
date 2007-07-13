#if !defined(__DATA_MOD_H)
#define __DATA_MOD_H

struct datafilehead_byte
{
  char nblocks[4];
  char ntraces[4];
  char np[4];
  char ebytes[4];
  char tbytes[4];
  char bbytes[4];
  char vers_id[2];
  char status[2];
  char nbheaders[4];
};

struct datablockhead_byte
{
   char scale[2];
   char status[2];
   char index[2];
   char mode[2];
   char ctcount[4];
   char lpval[4];
   char rpval[4];
   char lvl[4];
   char tlt[4];
};

#endif
