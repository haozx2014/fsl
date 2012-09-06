// Declarations of virtual base class for
// making predictions about DWI data.
//
// DWIPredictionMaker.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef DWIPredictionMaker_h
#define DWIPredictionMaker_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

namespace EDDY {

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class DWIPredictionMaker
//
// Virtual base class for classes used to make 
// predictions about diffusion data.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class DWIPredictionMaker
{
public:
  DWIPredictionMaker() {}
  virtual ~DWIPredictionMaker() {}
  virtual NEWIMAGE::volume<float> Predict(unsigned int indx) const = 0;
  virtual NEWIMAGE::volume<float> Predict(const DiffPara& dpar) const = 0;
  virtual NEWIMAGE::volume4D<float> PredictAll() const = 0;
  virtual bool IsValid() const = 0;
  virtual void SetNoOfScans(unsigned int n) = 0;
  virtual void AddScan(const NEWIMAGE::volume<float>& scan,  // NOT thread safe
                       const DiffPara&                dp) = 0;
  virtual void SetScan(const NEWIMAGE::volume<float>& scan,  // May be thread safe if used "sensibly"
                       const DiffPara&                dp,
                       unsigned int                   indx) = 0;
  virtual void EvaluateModel(const NEWIMAGE::volume<float>& mask) = 0;
};

// }}} End of fold.

} // End namespace EDDY

#endif // End #ifndef DWIPredictionMaker_h
