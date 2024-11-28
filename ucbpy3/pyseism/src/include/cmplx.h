#ifndef CMPLXDEF
#define CMPLXDEF 1
typedef struct {float real,imag;} complex;
#endif

float	c_amp();
complex *c_cos(),*c_div(),*c_mult(),*c_add();
complex *c_sub(),*c_exp(),*c_log(),*c_sin();
complex *c_sqrt();
complex *c_conj();
