#if !defined(__LIBCEPHES_H)
#define __LIBCEPHES_H

namespace CEPHES {

  double zetac ( double x );
  double cephes_j0(double x);
  double cephes_y0(double x);
  double cephes_j1(double x);
  double cephes_y1(double x);
  double cephes_i0(double x);
  double cephes_i0e(double x);
  double cephes_i1(double x);
  double cephes_i1e(double x);
  double cephes_chbevl(double x,double array[],int  n);

}

#endif
