%module pyseism_driver

%{
        #define SWIG_FILE_WITH_INIT
        #include "pyseism_driver.h"
%}

%include "numpy.i"

%init %{
        import_array();
%}

%apply ( double* INPLACE_ARRAY1, int DIM1 ) { ( double* u_t_, int ndata ) }

%include "pyseism_driver.h"
