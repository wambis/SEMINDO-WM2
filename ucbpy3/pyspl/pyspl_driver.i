%module pyspl_driver

%{
        #define SWIG_FILE_WITH_INIT
        #include "spl_driver.h"
%}

%include "numpy.i"

%init %{
        import_array();
%}

%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* knots, int nknots ) }
%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* r, int nr ) }
%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* knot_lons, int nknot_lons ) }
%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* knot_lats, int nknot_lats ) }
%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* knot_dists, int nknot_dists ) }
%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* lons, int nlons ) }
%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* lats, int nlats ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* res, int nres ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* resval, int nresval ) }
%apply ( int* INPLACE_ARRAY1, int DIM1 ) { ( int* resind, int nresind ) }
%apply ( int* INPLACE_ARRAY1, int DIM1 ) { ( int* resloc, int nresloc ) }

%include "spl_driver.h"
