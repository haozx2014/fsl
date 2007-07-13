/*  test.cc
    
    Copyright (C) 1999-2004 University of Oxford */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 3.3 (c) 2006, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */

#include <iostream>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include <string>
#include "utils/log.h"
#include "melgmix.h"
#include "meloptions.h"
#include "melhlprfns.h"
#include <time.h>
#include "miscmaths/miscprob.h"


using namespace std;
using namespace Utilities;
using namespace NEWIMAGE;
using namespace Melodic;


void usage(void)
{
  cout << "Usage: test infile ICfile [mixfile]" << endl;
  exit(1);
}

int main(int argc, char *argv[])
{
  Matrix Mats;

  Matrix A,C;
  RowVector B;

  Mats = normrnd(5,10);
  std_pca(Mats,C,A,B);

  cout << " STD PCA : " << endl << A << endl << endl << B << endl << endl;
  em_pca(Mats, A, B, 4);
  cout << "  EM PCA : " << endl << A << endl << endl << B << endl << endl;

  string RAWfname;
  RAWfname = string(argv[1]);

  volume4D<float> RawData;
  volumeinfo VolInfo;
  cout << " Reading orig. data " << RAWfname << " ... " << endl<< endl;
  read_volume4D(RawData,RAWfname,VolInfo);
  Matrix mat;
  
  volume<float> RawData2;
  read_volume(RawData2,RAWfname);

  print_size(RawData2);

  A = normrnd(7,3);
  C = normrnd(5,3);

  cout << " here " << endl;
  Mats = krprod(A,C);

  cout << " A : " << A << endl << endl;
  cout << " B : " << C << endl << endl;
  
  cout << " C : " << Mats << endl << endl;
  
  write_ascii_matrix(string("A"),A);
  write_ascii_matrix(string("B"),C);
  write_ascii_matrix(string("C"),Mats);

  Matrix tmpA, tmpB;
  tmpA = zeros(A.Nrows(),A.Ncols());
  tmpB = zeros(C.Nrows(),C.Ncols());

  krfact(Mats, tmpA, tmpB);

  cout << " A2 : " << tmpA << endl << endl; 
  cout << " B2 : " << tmpB << endl << endl; 
  write_ascii_matrix(string("A2"),tmpA);
  write_ascii_matrix(string("B2"),tmpB);

  tmpA = krapprox(Mats,tmpA.Nrows());
  write_ascii_matrix(string("C2"),tmpA);
  
  return 0;
}
