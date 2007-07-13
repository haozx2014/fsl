//
// Definitions for class BFMatrix
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newmatio.h"
#include "bfmatrix.h"

namespace MISCMATHS {

//
// Member functions for BFMatrix
//

void BFMatrix::print(const NEWMAT::Matrix&  m,
                     const std::string&     fname) const
{
  if (!fname.length()) {
    cout << endl << m << endl;
  }
  else {
    try {
      ofstream  fout(fname.c_str());
      fout << setprecision(10) << m;
    }
    catch(...) {
      std::string  errmsg("BFMatrix::print: Failed to write to file " + fname);
      throw BFMatrixException(errmsg);
    }
  }
}

//
// Member functions for SparseBFMatrix
//

//
// Concatenation of two matrices returning a third
//
void SparseBFMatrix::HorConcat(const BFMatrix& B, BFMatrix& AB) const
{
  try {
    const SparseBFMatrix& lB = dynamic_cast<const SparseBFMatrix&>(B);
    SparseBFMatrix& lAB = dynamic_cast<SparseBFMatrix&>(AB);
    if (Nrows() != lB.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat: Matrices must have same # of rows");}
    *(lAB.mp) = *mp | *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::HorConcat: dynamic cast error"); 
  }
}

void SparseBFMatrix::HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  try {
    SparseBFMatrix& lAB = dynamic_cast<SparseBFMatrix&>(AB);
    if (int(Nrows()) != B.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat: Matrices must have same # of rows");}
    *(lAB.mp) = *mp | B;
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::HorConcat: dynamic cast error"); 
  }
}

void SparseBFMatrix::VertConcat(const BFMatrix& B, BFMatrix& AB) const
{
  try {
    const SparseBFMatrix& lB = dynamic_cast<const SparseBFMatrix&>(B);
    SparseBFMatrix& lAB = dynamic_cast<SparseBFMatrix&>(AB);
    if (Ncols() != lB.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcat: Matrices must have same # of columns");}
    *(lAB.mp) = *mp & *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::VertConcat: dynamic cast error"); 
  }
}

void SparseBFMatrix::VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  try {
    SparseBFMatrix& lAB = dynamic_cast<SparseBFMatrix&>(AB);
    if (int(Ncols()) != B.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcat: Matrices must have same # of columns");}
    *(lAB.mp) = *mp & B;
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::VertConcat: dynamic cast error"); 
  }
}

//
// Concatenate another matrix to *this
//
void SparseBFMatrix::HorConcat2MyRight(const BFMatrix& B)
{
  try {
    const SparseBFMatrix& lB = dynamic_cast<const SparseBFMatrix&>(B);
    if (Nrows() != lB.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}
    *mp |= *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::VertConcat: dynamic cast error"); 
  }
}

void SparseBFMatrix::HorConcat2MyRight(const NEWMAT::Matrix& B)
{
  if (int(Nrows()) != B.Nrows()) {throw BFMatrixException("SparseBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}
  *mp |= B;
}

void SparseBFMatrix::VertConcatBelowMe(const BFMatrix& B)
{
  try {
    const SparseBFMatrix& lB = dynamic_cast<const SparseBFMatrix&>(B);
    if (Ncols() != lB.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}
    *mp &= *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::VertConcat: dynamic cast error"); 
  }
}

void SparseBFMatrix::VertConcatBelowMe(const NEWMAT::Matrix& B)
{
  if (int(Ncols()) != B.Ncols()) {throw BFMatrixException("SparseBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}
  *mp &= B;
}

// Multiply by vector
NEWMAT::ReturnMatrix SparseBFMatrix::MulByVec(const NEWMAT::ColumnVector& invec) const
{
  if (invec.Nrows() != int(Ncols())) {throw BFMatrixException("Matrix-vector size mismatch");}
  NEWMAT::ColumnVector   outvec = *mp * invec;
  outvec.Release();
  return(outvec);
}

// Add another matrix to this one
void SparseBFMatrix::AddToMe(const BFMatrix& M, double s)
{
  try {
    const SparseBFMatrix& lM = dynamic_cast<const SparseBFMatrix&>(M);
    if (Ncols() != lM.Ncols() || Nrows() != lM.Nrows()) {
      throw BFMatrixException("SparseBFMatrix::AddToMe: Matrix size mismatch");
    }
    if (s == 1.0) *mp += *(lM.mp);
    else *mp += s * *(lM.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("SparseBFMatrix::AddToMe: dynamic cast error"); 
  }
}

// Given A*x=b, solve for x
NEWMAT::ReturnMatrix SparseBFMatrix::SolveForx(const NEWMAT::ColumnVector& b,
					       MISCMATHS::MatrixType       type,
					       double                      tol,
                                               int                         miter) const
{
  if (b.Nrows() != int(Nrows())) {
    throw BFMatrixException("SparseBFMatrix::SolveForx: Matrix-vector size mismatch");
  }
  NEWMAT::ColumnVector  x = mp->SolveForx(b,type,tol,miter);
  x.Release();
  return(x);
}

//
// Member functions for FullBFMatrix
//

boost::shared_ptr<BFMatrix> FullBFMatrix::Transpose(boost::shared_ptr<BFMatrix>& pA) const
{
  boost::shared_ptr<FullBFMatrix>  tmp(new FullBFMatrix);
  *(tmp->mp) = mp->t();
  return(tmp);
}

//
// Concatenate two matrices yielding a third
//

void FullBFMatrix::HorConcat(const BFMatrix& B, BFMatrix& AB) const
{
  try {
    const FullBFMatrix& lB = dynamic_cast<const FullBFMatrix&>(B);
    FullBFMatrix& lAB = dynamic_cast<FullBFMatrix&>(AB);
    if (Nrows() != lB.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat: Matrices must have same # of rows");}
    *(lAB.mp) = *mp | *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::HorConcat: dynamic cast error"); 
  }
}

void FullBFMatrix::HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  try {
    FullBFMatrix& lAB = dynamic_cast<FullBFMatrix&>(AB);
    if (int(Nrows()) != B.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat: Matrices must have same # of rows");}
    *(lAB.mp) = *mp | B;
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::HorConcat: dynamic cast error"); 
  }
}

void FullBFMatrix::VertConcat(const BFMatrix& B, BFMatrix& AB) const
{
  try {
    const FullBFMatrix& lB = dynamic_cast<const FullBFMatrix&>(B);
    FullBFMatrix& lAB = dynamic_cast<FullBFMatrix&>(AB);
    if (Ncols() != lB.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcat: Matrices must have same # of columns");}
    *(lAB.mp) = *mp & *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::VertConcat: dynamic cast error"); 
  }
}

void FullBFMatrix::VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const
{
  try {
    FullBFMatrix& lAB = dynamic_cast<FullBFMatrix&>(AB);
    if (int(Ncols()) != B.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcat: Matrices must have same # of columns");}
    *(lAB.mp) = *mp & B;
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::VertConcat: dynamic cast error"); 
  }
}

//  
// Concatenation of another matrix to *this
//
void FullBFMatrix::HorConcat2MyRight(const BFMatrix& B)
{
  try {
    const FullBFMatrix& lB = dynamic_cast<const FullBFMatrix&>(B);
    if (Nrows() != lB.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}
    *mp |= *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::HorConcat2MyRight: dynamic cast error"); 
  }
}

void FullBFMatrix::HorConcat2MyRight(const NEWMAT::Matrix& B)
{
  if (int(Nrows()) != B.Nrows()) {throw BFMatrixException("FullBFMatrix::HorConcat2MyRight: Matrices must have same # of rows");}
  *mp |= B;
}

void FullBFMatrix::VertConcatBelowMe(const BFMatrix& B)
{
  try {
    const FullBFMatrix& lB = dynamic_cast<const FullBFMatrix&>(B);
    if (Ncols() != lB.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}
    *mp &= *(lB.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::HorConcatBelowMe: dynamic cast error"); 
  }
}

void FullBFMatrix::VertConcatBelowMe(const NEWMAT::Matrix& B)
{
  if (int(Ncols()) != B.Ncols()) {throw BFMatrixException("FullBFMatrix::VertConcatBelowMe: Matrices must have same # of columns");}
  *mp &= B;
}

// Multiply this matrix with scalar

void FullBFMatrix::MulMeByScalar(double s)
{
  *mp = s*(*mp);
}

// Multiply by vector
NEWMAT::ReturnMatrix FullBFMatrix::MulByVec(const NEWMAT::ColumnVector& invec) const
{
  if (invec.Nrows() != int(Ncols())) {throw BFMatrixException("FullBFMatrix::MulByVec: Matrix-vector size mismatch");}
  NEWMAT::ColumnVector  ret;
  ret = (*mp)*invec;
  ret.Release();
  return(ret);
}

// Add another matrix to this one
void FullBFMatrix::AddToMe(const BFMatrix& m, double s)
{
  try {
    const FullBFMatrix& lm = dynamic_cast<const FullBFMatrix&>(m);
    if (Ncols() != lm.Ncols() || Nrows() != lm.Nrows()) {
      throw BFMatrixException("FullBFMatrix::AddToMe: Matrix size mismatch");
    }
    *mp += s*(*lm.mp);
  }
  catch (std::bad_cast) {
    throw BFMatrixException("FullBFMatrix::AddToMe: dynamic cast error"); 
  }
}

// Given A*x=b, solve for x
NEWMAT::ReturnMatrix FullBFMatrix::SolveForx(const NEWMAT::ColumnVector& b,       // Ignoring all parameters except b
					     MISCMATHS::MatrixType       type,
					     double                      tol,
                                             int                         miter) const
{
  if (int(Nrows()) != b.Nrows()) {throw BFMatrixException("FullBFMatrix::SolveForx: Matrix-vector size mismatch");}
  NEWMAT::ColumnVector  ret;
  ret = mp->i()*b;
  ret.Release();
  return(ret);
}
    
} // End namespace MISCMATHS
