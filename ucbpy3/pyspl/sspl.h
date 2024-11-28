
#ifndef _SSPL_H

#define _SSPL_H

double
getdel(double a0, double o0, double a, double o);

void
getdel_vec(double *a0_vec, double *o0_vec, int N, double a, double o, double *res);

double 
spbsp(double hdel, double ahdel);

double 
gspbsp(double hdel, double ahdel);

#endif
