// Declarations of classes that implements a scan
// or a collection of scans within the EC project.
// 
// EddyScanClasses.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef ECScanClasses_h
#define ECScanClasses_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "topup/topup_file_io.h"
#include "EddyHelperClasses.h"
#include "ECModels.h"

namespace EDDY {

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class ECScan
//
// This class manages one diffusion weighted scan
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class ECScan
{
public:
  ECScan(const NEWIMAGE::volume<float>    ima,
         const AcqPara&                   acqp,
         const DiffPara&                  diffp,
         boost::shared_ptr<ScanECModel>   ecp,
         unsigned int                     sess) : _ima(ima), _acqp(acqp), _diffp(diffp), _ecp(ecp->Clone()), _fwhm(0.0), _sess(sess){}
  ECScan(const ECScan& inp) : _ima(inp._ima), _sima(inp._sima), _acqp(inp._acqp), _diffp(inp._diffp), _ecp(inp._ecp->Clone()), _fwhm(inp._fwhm), _sess(inp._sess) {}
  virtual ~ECScan() {}
  ECScan& operator=(const ECScan& rhs) {
    if (this == &rhs) return(*this);
    _ima=rhs._ima; _sima=rhs._sima; _acqp=rhs._acqp; _diffp=rhs._diffp; _ecp=rhs._ecp->Clone(); _fwhm=rhs._fwhm; _sess=rhs._sess;
    return(*this);
  }
  unsigned int Session() const { return(_sess); }
  bool HasFieldOffset() const { return(_ecp->HasFieldOffset()); }
  double GetFieldOffset() const { return(_ecp->GetFieldOffset()); }
  void SetFieldOffset(double ofst) { _ecp->SetFieldOffset(ofst); }
  virtual NEWMAT::RowVector GetLinearParameters() const { return(_ecp->GetLinearParameters()); }
  const NEWIMAGE::volume<float>& GetOriginalIma() const { return(_ima); }
  NEWIMAGE::volume<float> GetMotionCorrectedOriginalIma(NEWIMAGE::volume<float>& omask) const { return(motion_correct(GetOriginalIma(),&omask)); }
  NEWIMAGE::volume<float> GetMotionCorrectedOriginalIma() const { return(motion_correct(GetOriginalIma(),NULL)); }
  NEWIMAGE::volume<float> GetUnwarpedOriginalIma(// Input
						 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
						 // Output
						 NEWIMAGE::volume<float>&                            omask) const;
  NEWIMAGE::volume<float> GetUnwarpedOriginalIma(// Input
						 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc) const;
  const NEWIMAGE::volume<float>& GetIma() const { return( (_fwhm) ? _sima : _ima); } 
  NEWIMAGE::volume<float> GetUnwarpedIma(// Input
					 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
					 // Output
					 NEWIMAGE::volume<float>&                            omask) const;
  NEWIMAGE::volume<float> GetUnwarpedIma(// Input
					 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc) const;
  AcqPara GetAcqPara() const { return(_acqp); }
  DiffPara GetDiffPara() const { return(_diffp); }
  NEWMAT::ColumnVector GetHz2mmVector() const;
  unsigned int NParam(EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->NParam(whichp)); }
  NEWMAT::ColumnVector GetParams(EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->GetParams(whichp)); }
  unsigned int NDerivs(EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->NDerivs(whichp)); }
  double GetDerivParam(unsigned int indx, EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->GetDerivParam(indx,whichp)); }
  double GetDerivScale(unsigned int indx, EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->GetDerivScale(indx,whichp)); }
  void SetDerivParam(unsigned int indx, double p, EDDY::Parameters whichp=EDDY::ALL) { _ecp->SetDerivParam(indx,p,whichp); }
  double GetFWHM() const { return(_fwhm); }
  // Returns matrix denoted \mathbf{R} in paper.
  NEWMAT::Matrix ForwardMovementMatrix() const { return(_ecp->ForwardMovementMatrix(_ima)); }
  // Returns matrix denoted \mathbf{R}^{-1} in paper.
  NEWMAT::Matrix InverseMovementMatrix() const { return(_ecp->InverseMovementMatrix(_ima)); }
  // Returns the EC-field relevant for the Observation->Model transform
  NEWIMAGE::volume<float> ECField() const { return(_ecp->ECField(_ima)); }

  // Returns the total field relevant for the Observation->Model transform
  NEWIMAGE::volume4D<float> FieldForScanToModelTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask,
							 NEWIMAGE::volume<float>&                          jac) const {
    return(field_for_scan_to_model_transform(susc,&omask,&jac));
  }
  NEWIMAGE::volume4D<float> FieldForScanToModelTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask) const {
    return(field_for_scan_to_model_transform(susc,&omask,NULL));
  }
  NEWIMAGE::volume4D<float> FieldForScanToModelTransform(boost::shared_ptr<const NEWIMAGE::volume<float> > susc) const {
    return(field_for_scan_to_model_transform(susc,NULL,NULL));
  }

  // Returns the total field relevant for the Model->Observation transform
  NEWIMAGE::volume4D<float> FieldForModelToScanTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask,
							 NEWIMAGE::volume<float>&                          jac) const {
    return(field_for_model_to_scan_transform(susc,&omask,&jac));
  }
  NEWIMAGE::volume4D<float> FieldForModelToScanTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask) const {
    return(field_for_model_to_scan_transform(susc,&omask,NULL));
  }
  NEWIMAGE::volume4D<float> FieldForModelToScanTransform(boost::shared_ptr<const NEWIMAGE::volume<float> > susc) const {
    return(field_for_model_to_scan_transform(susc,NULL,NULL));
  }
  NEWIMAGE::volume4D<float> FieldForModelToScanTransformWithJac(// Input
								boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
								// Output
								NEWIMAGE::volume<float>&                          jac) const {
    return(field_for_model_to_scan_transform(susc,NULL,&jac));
  }
  void SetParams(const NEWMAT::ColumnVector& mpep, EDDY::Parameters whichp=EDDY::ALL) { _ecp->SetParams(mpep,whichp); }
  void SetFWHM(double fwhm);
  void ReplaceSlices(const NEWIMAGE::volume<float>&                     rep,
		     boost::shared_ptr<const NEWIMAGE::volume<float> >  susc,
		     const NEWIMAGE::volume<float>&                     inmask, 
		     const std::vector<unsigned int>&                   ol);
private:
  NEWIMAGE::volume<float>         _ima;
  NEWIMAGE::volume<float>         _sima;
  AcqPara                         _acqp;
  DiffPara                        _diffp;
  boost::shared_ptr<ScanECModel>  _ecp;
  double                          _fwhm;
  unsigned int                    _sess;

  NEWIMAGE::volume<float> motion_correct(// Input
					 const NEWIMAGE::volume<float>&  inima,
					 // Output (optional)
					 NEWIMAGE::volume<float>         *omask) const;
  NEWIMAGE::volume<float> transform_to_model_space(// Input
						   const NEWIMAGE::volume<float>&                    inima,
						   boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
						   // Output
						   NEWIMAGE::volume<float>&                          omask) const;
  NEWIMAGE::volume4D<float> field_for_scan_to_model_transform(// Input
							      boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							      // Output
							      NEWIMAGE::volume<float>                           *omask,
							      NEWIMAGE::volume<float>                           *jac) const;
  NEWIMAGE::volume4D<float> field_for_model_to_scan_transform(// Input
							      boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							      // Output
							      NEWIMAGE::volume<float>                           *omask,
							      NEWIMAGE::volume<float>                           *jac) const;
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class ECScanManager
//
// This class manages all diffusion weighted scans.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class ECScanManager
{
public:
  ECScanManager(const std::string&               simafname,
                const std::string&               maskfname,
                const std::string&               acqfname,
                const std::string&               topupfname,
                const std::string&               bvecfname,
                const std::string&               bvalsfname,
                EDDY::ECModel                    ecmodel,
                const std::vector<unsigned int>& indicies,
		const std::vector<unsigned int>& session_indicies,
                unsigned int                     no_of_sessions);
  ~ECScanManager() {}
  // Number of scans of type st
  unsigned int NScans(ScanType st=ANY) const;
  // Number of LSR pairs of type st
  unsigned int NLSRPairs(ScanType st=ANY) const;
  // Return number of sessions over which data was acquired
  unsigned int NoOfSessions() const { return(_nsess); }
  bool IsDWI(unsigned int indx) const { if (indx>_fi.size()-1) throw EddyException("ECScanManager::IsDWI: index out of range"); else return(!_fi[indx].first); }
  bool IsB0(unsigned int indx) const { return(!IsDWI(indx)); }
  double ScaleFactor() const { return(_sf); }
  // Returns vector of Global (file) indicies for all dwi indicies.
  std::vector<unsigned int> GetDwi2GlobalIndexMapping() const;
  // Returns true if a susceptibilty induced off-resonance field has been set
  bool HasSuscHzOffResField() const { return(_has_topup); }
  // Returns true if the EC model includes a field offset
  bool HasFieldOffset() const { return(Scan(0,DWI).HasFieldOffset()); }
  // Returns true if data allows for LSR resampling
  bool CanDoLSRResampling() const;
  // Returns the individual indicies for scans that constitute the i'th LSR pair
  std::pair<unsigned int,unsigned int> GetLSRPair(unsigned int i, ScanType st) const;
  // Sets parameters for all scans
  void SetParameters(const NEWMAT::Matrix& pM, ScanType st=ANY);
  void SetParameters(const std::string& fname, ScanType st=ANY) { NEWMAT::Matrix pM = MISCMATHS::read_ascii_matrix(fname); SetParameters(pM,st); }
  // Returns a read-only reference to a scan given by indx
  const ECScan& Scan(unsigned int indx, ScanType st=ANY) const;
  // Returns a read-write reference to a scan given by indx
  ECScan& Scan(unsigned int indx, ScanType st=ANY);
  // Returns an "original" (no smoothing) "unwarped" into model space
  NEWIMAGE::volume<float> GetUnwarpedOrigScan(unsigned int              indx,
					      NEWIMAGE::volume<float>&  omask,
					      ScanType                  st=ANY) const; 
  NEWIMAGE::volume<float> GetUnwarpedOrigScan(unsigned int indx,
					      ScanType     st=ANY) const { 
    NEWIMAGE::volume<float> mask=_scans[0].GetIma(); mask=1.0; 
    return(GetUnwarpedOrigScan(indx,mask,st));
  }
  // Returns an image "unwarped" into model space
  NEWIMAGE::volume<float> GetUnwarpedScan(unsigned int              indx,
					  NEWIMAGE::volume<float>&  omask,
					  ScanType                  st=ANY) const; 
  NEWIMAGE::volume<float> GetUnwarpedScan(unsigned int indx,
					  ScanType     st=ANY) const { 
    NEWIMAGE::volume<float> mask=_scans[0].GetIma(); mask=1.0; 
    return(GetUnwarpedScan(indx,mask,st));
  }
  // Resamples a matching pair of images using least-squares resampling.
  NEWIMAGE::volume<float> LSRResamplePair(// Input
					  unsigned int              i, 
					  unsigned int              j, 
					  ScanType                  st,
					  // Output
					  NEWIMAGE::volume<float>&  omask) const;
  // Sets movement and EC parameters for a scan given by indx
  void SetScanParameters(unsigned int indx, const NEWMAT::ColumnVector& p, ScanType st=ANY) { Scan(indx,st).SetParams(p,ALL); }

  // Smooths all scans to the requested FWHM
  void SetFWHM(double fwhm) {
    # pragma omp parallel for shared(fwhm)
    for (int i=0; i<int(_scans.size()); i++) {
      _scans[i].SetFWHM(fwhm);
    } 
    # pragma omp parallel for shared(fwhm)
    for (int i=0; i<int(_b0scans.size()); i++) {
      _b0scans[i].SetFWHM(fwhm);
    } 
  }

  // Returns the current smoothness of (all) the scans
  double GetFWHM() const { if (_scans.size()) return(_scans[0].GetFWHM());  else if (_b0scans.size()) return(_b0scans[0].GetFWHM()); else return(0.0); }

  // Returns the user defined mask
  const NEWIMAGE::volume<float>& Mask() const { return(_mask); }

  // Returns a pointer to susceptibility induced off-resonance field (in Hz). Returns NULL when no field set.
  const boost::shared_ptr<NEWIMAGE::volume<float> > GetSuscHzOffResField() const { return(_topup_field); }

  // Returns off-resonance field pertaining to EC only. Movements have not been taken into
  // account and its main use is for visualiation/demonstration.
  NEWIMAGE::volume<float> GetScanHzECOffResField(unsigned int indx, ScanType st=ANY) const { return(Scan(indx,st).ECField()); }

  // Attempt to distinguish field offset (which may be due to image FOV centre not
  // coinciding with scanner iso-centre) for actual subject movement in PE direction.
  void SeparateFieldOffsetFromMovement();

  // Set specified scan as reference for location.
  void SetDWIReference(unsigned int ref=0) { set_reference(ref,DWI); }
  void Setb0Reference(unsigned int ref=0) { set_reference(ref,B0); }

  // Various methods for writing results to disk
  void WriteRegisteredImages(const std::string& fname, FinalResampling  resmethod, ScanType st=ANY)
  {
    if (resmethod==EDDY::JAC) write_jac_registered_images(fname,st);
    else if (resmethod==EDDY::LSR) write_lsr_registered_images(fname,st);
    else throw EddyException("ECScanManager::WriteRegisteredImages: Unknown resampling method");
  }
  void WriteParameterFile(const std::string& fname, ScanType st=ANY) const;
  void WriteECFields(const std::string& fname, ScanType st=ANY) const;
private:
  bool                                          _has_topup;
  boost::shared_ptr<NEWIMAGE::volume<float> >   _topup_field;
  NEWIMAGE::volume<float>                       _mask;
  double                                        _sf;           // Scale factor applied to scans
  std::vector<pair<int,int> >                   _fi;           // Used to keep track of index into file
  std::vector<ECScan>                           _scans;
  std::vector<ECScan>                           _b0scans;
  unsigned int                                  _nsess;

  NEWMAT::ColumnVector hz_vector_with_everything() const;
  NEWMAT::Matrix offset_design_matrix() const;
  void set_reference(unsigned int ref, ScanType st);
  double mean_of_first_b0(const NEWIMAGE::volume4D<float>&   vols,
                          const NEWIMAGE::volume<float>&     mask,
                          const NEWMAT::Matrix&              bvecs,
                          const NEWMAT::Matrix&              bvals) const;
  bool index_kosher(unsigned int indx, ScanType st) const
  {
    if (st==DWI) return(indx<_scans.size());
    else if (st==B0) return(indx<_b0scans.size());
    else return(indx<_fi.size());
  }
  void write_jac_registered_images(const std::string& fname, ScanType st);
  void write_lsr_registered_images(const std::string& fname, ScanType st);
};

} // End namespace EDDY

#endif // End #ifndef ECScanClasses_h
