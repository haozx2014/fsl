//     fslcorrecthd.cc - check and correct a nifti file for bad vox-offset
//     Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2007 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include <iostream>
using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslcorrecthd <input> <output>" << endl;
}

int main(int argc,char *argv[])
{
  if (argc < 3) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  FSLIO* fslio=NULL;
  fslio = FslOpen(FslMakeBaseName(argv[1]),"rb");
  FslClose(fslio);
  int ft;
  struct dsr *hdr;
  ft = FslGetFileType(fslio);
  hdr = (struct dsr *)calloc(1,sizeof(struct dsr));
  FslReadRawHeader(hdr,fslio->niftiptr->fname);
  if (fslio->niftiptr->byteorder != nifti_short_order()) 
  {
    cout << "Byte swapping" << endl;
    AvwSwapHeader(hdr);
  } 
  //check nifti-libs output versus raw header info
  ft = fslio->niftiptr->iname_offset - hdr->dime.vox_offset;
  int minft=MIN(fslio->niftiptr->iname_offset,hdr->dime.vox_offset);
  cerr << "number of bytes wrong: " << ft << endl;
  cerr << "start at byte location: " << minft << endl;

 ifstream input_file;
 ofstream output_file;
 char byte[1];
 input_file.open(argv[1],ios::in | ios :: binary);
 output_file.open(argv[2],ofstream::out | ofstream::binary);
 for(int i=1;i<=minft;i++)
 {
   input_file.read(byte,1);
   if (input_file.eof()) break;
   output_file.write(byte,1);
 }

 for(int i=1;i<=abs(ft) && ft>0;i++)
 {
   byte[0]=0;
   output_file.write(byte,1);
   }

 for(int i=1;i<=abs(ft) && ft<0;i++)
 {
   input_file.read(byte,1);
 }

 while(true)
 {
   input_file.read(byte,1);
   if (input_file.eof()) break;
   output_file.write(byte,1);
 }  
 output_file.close();
 input_file.close();
  return 0;
}


