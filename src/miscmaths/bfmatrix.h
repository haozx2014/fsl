// Declarations for class BFMatrix.
//
// The purpose of class BFmatrix is to have a class from which
// to derive 2 other classes; FullBFMatrix and SparseBFMatrix.
// The reason for this is that the two classes SplineField and
// DCTField will return Hessian matrices that are either Sparse
// (SplineField) or full (DCTField). By defining a pure virtual
// class BFMatrix with a minimal (only what is needed for non-
// linear reg.) functionality I will be able to write code that
// is independent of type of matrix returned, and hence of type
// field.
//
// The syntax for the (little) functionality is sort of a mixture
// of Newmat and SparseMatrix. Mostly SparseMatrix actually.
// I hope this will not complicate the use of the nonlin package
// for those who are only interested in the full (normal) case.
//
// At one point SparseMatrix was replaced by SpMat as the underlying
// sparse matrix representation in SparseBFMatrix. SpMat was written
// with an API that largely mimicks that of NEWMAT. This is the 
// "historical" reason why a wrapper class was written, rather than
// using templatisation which would have been possible given the
// similarities in API between SpMat and NEWMAT.
//

#ifndef BFMatrix_h
#define BFMatrix_h

#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "SpMat.h"
#include "cg.h"
#include "bicg.h"

namespace MISCMATHS {


class BFMatrixException: public std::exception
{
private:
  std::string m_msg;
public:
  BFMatrixException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("BFMatrix::" + m_msg).c_str();
  }

  ~BFMatrixException() throw() {}
};

class BFMatrix
{
protected:
virtual void print(const NEWMAT::Matrix& m, const std::string& fname) const;
public:
  // Constructors, destructors and stuff
  BFMatrix() {}
  BFMatrix(unsigned int m, unsigned int n) {}
  virtual ~BFMatrix() {}

  // Access as NEWMAT::Matrix
  virtual NEWMAT::ReturnMatrix AsMatrix() const = 0;
  
  // Basic properties
  virtual unsigned int Nrows() const = 0;
  virtual unsigned int Ncols() const = 0;

  // Print matrix (for debugging)
  virtual void Print(const std::string fname=std::string("")) const = 0;

  // Setting, deleting or resizing the whole sparse matrix.
  virtual void SetMatrix(const MISCMATHS::SpMat<double>& M) = 0;
  virtual void SetMatrix(const NEWMAT::Matrix& M) = 0;
  virtual void Clear() = 0;
  virtual void Resize(unsigned int m, unsigned int n) = 0;

  // Accessing
  inline double operator()(unsigned int r, unsigned int c) const {return(Peek(r,c));}
  virtual double Peek(unsigned int r, unsigned int c) const = 0;

  // Assigning
  virtual void Set(unsigned int x, unsigned int y, double val) = 0;
  virtual void Insert(unsigned int x, unsigned int y, double val) = 0;
  virtual void AddTo(unsigned int x, unsigned int y, double val) = 0;  

  // Transpose
  virtual boost::shared_ptr<BFMatrix> Transpose(boost::shared_ptr<BFMatrix>& pA) const = 0;

  // Concatenation. Note that desired polymorphism prevents us from using BFMatrix->NEWMAT::Matrix conversion
  // Concatenate two matrices yielding a third
  // AB = [*this B] in Matlab lingo
  virtual void HorConcat(const BFMatrix& B, BFMatrix& AB) const = 0;
  virtual void HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const = 0;
  // AB = [*this; B] in Matlab lingo
  virtual void VertConcat(const BFMatrix& B, BFMatrix& AB) const = 0;
  virtual void VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const = 0;
  // Concatenate another matrix to *this
  virtual void HorConcat2MyRight(const BFMatrix& B) = 0;
  virtual void HorConcat2MyRight(const NEWMAT::Matrix& B) = 0;
  virtual void VertConcatBelowMe(const BFMatrix& B) = 0;
  virtual void VertConcatBelowMe(const NEWMAT::Matrix& B) = 0;

  // Multiply by scalar
  virtual void MulMeByScalar(double s) = 0;
  // Multiply by vector
  virtual NEWMAT::ReturnMatrix MulByVec(const NEWMAT::ColumnVector& v) const = 0;
  // Add another matrix to this one 
  virtual void AddToMe(const BFMatrix& m, double s=1.0) = 0;
  // Given A*x=b, solve for x.
  virtual NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector& b,
					 MISCMATHS::MatrixType       type=SYM_POSDEF,
					 double                      tol=1e-6,
                                         int                         miter=200) const = 0;  
};

class SparseBFMatrix : public BFMatrix
{
private:
  boost::shared_ptr<MISCMATHS::SpMat<double> >    mp;
  
public:
  // Constructors, destructor and assignment
  SparseBFMatrix() 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>())) {}
  SparseBFMatrix(unsigned int m, unsigned int n) 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(m,n))) {}
  SparseBFMatrix(unsigned int m, unsigned int n, const unsigned int *irp, const unsigned int *jcp, const double *sp)
  : mp(boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(m,n,irp,jcp,sp))) {}
  SparseBFMatrix(const MISCMATHS::SpMat<double>& M) 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(M))) {}
  SparseBFMatrix(const NEWMAT::Matrix& M) 
  : mp(boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(M))) {}
  virtual ~SparseBFMatrix() {}
  virtual const SparseBFMatrix& operator=(const SparseBFMatrix& M) {
    mp = boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(*(M.mp))); return(*this);
  }

  // Access as NEWMAT::Matrix
  virtual NEWMAT::ReturnMatrix AsMatrix() const {NEWMAT::Matrix ret; ret = mp->AsNEWMAT(); ret.Release(); return(ret);}
  
  // Basic properties
  virtual unsigned int Nrows() const {return(mp->Nrows());}
  virtual unsigned int Ncols() const {return(mp->Ncols());}

  // Print matrix (for debugging)
  virtual void Print(const std::string fname=std::string("")) const {print(mp->AsNEWMAT(),fname);}

  // Setting, deleting or resizing the whole sparse matrix.
  virtual void SetMatrix(const MISCMATHS::SpMat<double>& M) {mp = boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(M));}
  virtual void SetMatrix(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(M));}
  virtual void SetMatrixPtr(boost::shared_ptr<MISCMATHS::SpMat<double> >& mptr) {mp = mptr;}
  virtual void Clear() {mp = boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>());}
  virtual void Resize(unsigned int m, unsigned int n) {mp = boost::shared_ptr<MISCMATHS::SpMat<double> >(new MISCMATHS::SpMat<double>(m,n));}

  // Accessing values
  virtual double Peek(unsigned int r, unsigned int c) const {return(mp->Peek(r,c));}

  // Setting and inserting values
  virtual void Set(unsigned int x, unsigned int y, double val) {mp->Set(x,y,val);}
  virtual void Insert(unsigned int x, unsigned int y, double val) {mp->Set(x,y,val);}
  virtual void AddTo(unsigned int x, unsigned int y, double val) {mp->AddTo(x,y,val);}

  // Transpose. 
  virtual boost::shared_ptr<BFMatrix> Transpose(boost::shared_ptr<BFMatrix>& pA) const {throw BFMatrixException("SparseBFMatrix::Transpose: Not yet implemented");}
  
  // Concatenation of two matrices returning a third
  // AB = [*this B] in Matlab lingo
  virtual void HorConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;
  // AB = [*this; B] in Matlab lingo
  virtual void VertConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;

  // Concatenation of another matrix to *this
  virtual void HorConcat2MyRight(const BFMatrix& B);
  virtual void HorConcat2MyRight(const NEWMAT::Matrix& B);
  virtual void VertConcatBelowMe(const BFMatrix& B);
  virtual void VertConcatBelowMe(const NEWMAT::Matrix& B);

  // Multiply by scalar
  virtual void MulMeByScalar(double s) {(*mp)*=s;}

  // Multiply by vector
  virtual NEWMAT::ReturnMatrix MulByVec(const NEWMAT::ColumnVector& invec) const;

  // Add another matrix to this one
  virtual void AddToMe(const BFMatrix& m, double s=1.0);

  // Given A*x=b, solve for x
  virtual NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector& b,
					 MISCMATHS::MatrixType       type,
					 double                      tol,
                                         int                         miter) const;
};

class FullBFMatrix : public BFMatrix
{
private:
  boost::shared_ptr<NEWMAT::Matrix>    mp;
public:
  // Constructors, destructor and assignment
  FullBFMatrix() {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix());}
  FullBFMatrix(unsigned int m, unsigned int n) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(m,n));}
  FullBFMatrix(const MISCMATHS::SpMat<double>& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M.AsNEWMAT()));}
  FullBFMatrix(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M));}
  virtual ~FullBFMatrix() {}
  virtual const FullBFMatrix& operator=(const FullBFMatrix& M) {
    mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(*(M.mp))); return(*this);
  }

  virtual NEWMAT::ReturnMatrix AsMatrix() const {NEWMAT::Matrix ret; ret = *mp; ret.Release(); return(ret);}

  // Basic properties
  virtual unsigned int Nrows() const {return(mp->Nrows());}
  virtual unsigned int Ncols() const {return(mp->Ncols());}

  // Print matrix (for debugging)
  virtual void Print(const std::string fname=std::string("")) const {print(*mp,fname);}

  // Setting, deleting or resizing the whole matrix.
  virtual void SetMatrix(const MISCMATHS::SpMat<double>& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M.AsNEWMAT()));} 
  virtual void SetMatrix(const NEWMAT::Matrix& M) {mp = boost::shared_ptr<NEWMAT::Matrix>(new NEWMAT::Matrix(M));}
  virtual void SetMatrixPtr(boost::shared_ptr<NEWMAT::Matrix>& mptr) {mp = mptr;}
  virtual void Clear() {mp->ReSize(0,0);}
  virtual void Resize(unsigned int m, unsigned int n) {mp->ReSize(m,n);}

  // Accessing values
  virtual double Peek(unsigned int r, unsigned int c) const {return((*mp)(r,c));}

  // Setting and inserting values.
  virtual void Set(unsigned int x, unsigned int y, double val) {(*mp)(x,y)=val;}
  virtual void Insert(unsigned int x, unsigned int y, double val) {(*mp)(x,y)=val;}
  virtual void AddTo(unsigned int x, unsigned int y, double val) {(*mp)(x,y)+=val;}

  // Transpose. 
  virtual boost::shared_ptr<BFMatrix> Transpose(boost::shared_ptr<BFMatrix>& pA) const;
  
  // Concatenation of two matrices returning a third
  virtual void HorConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void HorConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;
  virtual void VertConcat(const BFMatrix& B, BFMatrix& AB) const;
  virtual void VertConcat(const NEWMAT::Matrix& B, BFMatrix& AB) const;

  // Concatenation of another matrix to *this
  virtual void HorConcat2MyRight(const BFMatrix& B);
  virtual void HorConcat2MyRight(const NEWMAT::Matrix& B);
  virtual void VertConcatBelowMe(const BFMatrix& B);
  virtual void VertConcatBelowMe(const NEWMAT::Matrix& B);

  // Multiply by scalar
  virtual void MulMeByScalar(double s);

  // Multiply by vector
  virtual NEWMAT::ReturnMatrix MulByVec(const NEWMAT::ColumnVector& invec) const;

  // Add another matrix to this one
  virtual void AddToMe(const BFMatrix& m, double s);

  // Given A*x=b, solve for x
  virtual NEWMAT::ReturnMatrix SolveForx(const NEWMAT::ColumnVector& b,
					 MISCMATHS::MatrixType       type,
					 double                      tol,
                                         int                         miter) const;
    
};

} // End namespace MISCMATHS

#endif // End #ifndef BFMatrix_h
