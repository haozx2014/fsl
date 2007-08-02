#include <iostream>
#include <fstream>

using namespace std;

int main(int argc,char *argv[])
{
ifstream input_file;
ofstream output_file;
char byte[1];
 input_file.open(argv[1],ios::in | ios :: binary);
 output_file.open(argv[2],ofstream::out | ofstream::binary);
 for(int i=1;i<=348;i++)
 {
   input_file.read(byte,1);
   if (input_file.eof()) break;
   output_file.write(byte,1);
 }
 for(int i=1;i<=4;i++)
 {
   byte[0]=0;
   output_file.write(byte,1);
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
