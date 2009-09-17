/* BLAS function prototypes */

#ifndef BLAS_H
#define BLAS_H

#ifdef BLAS__SUNPRO_CC
#define dcopy dcopy_
#define daxpy daxpy_
#define dgemv dgemv_
#define dtrmv dtrmv_
#define dtrsv dtrsv_
#define dgemm dgemm_
#define dtrsm dtrsm_
#define dsyrk dsyrk_
#define dger dger_
#define dpotrf dpotrf_
#endif

#ifdef __cplusplus
extern "C" {
#endif 
  void dcopy( const int*, const double*, const int*, double*, const int*);
  void daxpy( const int*, const double*, const double*, const int*, double*, const int*);


  void dgemv( const char*, const int*, const int*, 
              const double*, const double*, const int*, 
              const double*, const int*, const double*, 
              double*, const int*);
  void dtrmv(const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
  void dtrsv(const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);


  void dtrsm( const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
  void dgemm(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, const double*, const int*) ;void dsyrk( const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, double*, const int*);
  
  void dpotrf( const char*, const int*, const double*, const int*, const int*) ;

#ifdef __cplusplus
}
#endif



#endif