%module pyminos_driver

%{
        #define SWIG_FILE_WITH_INIT
        #include "pyminos_driver.h"
%}

%include "numpy.i"

%init %{
        import_array();
%}

%apply ( double* IN_ARRAY1, int DIM1 ) { ( double* model_vec, int ndata ) }
%apply ( int* INPLACE_ARRAY1, int DIM1 ) { ( int* n_, int ndim_ ) }
%apply ( int* INPLACE_ARRAY1, int DIM1 ) { ( int* l_, int ldim_ ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* w_, int wdim_ ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* c_, int cdim_ ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* U_, int Udim_ ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* q_, int qdim_ ) }
%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* eig_, int eigdim_ ) }

%include "pyminos_driver.h"
