/*   single precision complex arthimatic routines.  Routines require
     the type structure "complex" to be initilized in main.
     Most routines return pointers to the answer.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <messageiolib.h>

typedef struct { float real, imag; } complex;
complex dummy;

float c_amp(z)
complex z;
{
float f;
	f=z.real*z.real+z.imag*z.imag;
	f=(float) sqrt( (double)f );
	return f;
}

/*
 *	"@(#)cabs.c	1.1"
 */

double cabs(real, imag)
double real, imag;
{
double temp, sqrt();

if(real < 0)
	real = -real;
if(imag < 0)
	imag = -imag;
if(imag > real){
	temp = real;
	real = imag;
	imag = temp;
}
if((real+imag) == real)
	return(real);

temp = imag/real;
temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
return(temp);
}

/*
 *	"@(#)c_cos.c	1.1"
 */
complex  *c_cos(z)
complex z;
{
double sin(), cos(), sinh(), cosh();

dummy.real = cos(z.real) * cosh(z.imag);
dummy.imag = - sin(z.real) * sinh(z.imag);
return(&dummy);
}

complex *c_div(a, b)
complex a, b;
{
float den;
        den=b.real*b.real+b.imag*b.imag;
        if(!den) 
          stop("c_div: division by zero");
        
	dummy.real=(a.real*b.real+a.imag*b.imag)/den;
	dummy.imag=(a.imag*b.real-a.real*b.imag)/den;
	return &dummy;
}


complex *c_mult(a,b)
complex a,b;
{
dummy.real = a.real * b.real - a.imag * b.imag;
dummy.imag = a.imag * b.real + a.real * b.imag;
return(&dummy);
}

complex *c_add(a,b)
complex a,b;
{
dummy.real = a.real + b.real;
dummy.imag = a.imag + b.imag;
return(&dummy);
}

complex *c_sub(a,b)
complex a,b;
{
dummy.real = a.real - b.real;
dummy.imag = a.imag - b.imag;
return(&dummy);
}


/*
 *	"@(#)c_exp.c	1.1"
 */
complex *c_exp(z)
complex z;
{
double expx;
double exp(), cos(), sin();

expx = exp(z.real);
dummy.real = expx * cos(z.imag);
dummy.imag = expx * sin(z.imag);
return(&dummy);
}


/*
 *	"@(#)c_log.c	1.1"
 */
complex *c_log(z)
complex z;
{
double log(), cabs(), atan2();

dummy.imag = atan2(z.imag, z.real);
dummy.real = log( cabs(z.real, z.imag) );
return(&dummy);
}


/*
 *	"@(#)c_sin.c	1.1"
 */
complex *c_sin(z)
complex z;
{
double sin(), cos(), sinh(), cosh();

dummy.real = sin(z.real) * cosh(z.imag);
dummy.imag = cos(z.real) * sinh(z.imag);
return(&dummy);
}


/*
 *	"@(#)c_sqrt.c	1.1"
 */
complex *c_sqrt(z)
complex z;
{
double mag, sqrt(), cabs();

if( (mag = cabs(z.real, z.imag)) == 0.)
	dummy.real = dummy.imag = 0.;
else if(z.real > 0)
	{
	dummy.real = sqrt(0.5 * (mag + z.real) );
	dummy.imag = z.imag / dummy.real / 2;
	}
else
	{
	dummy.imag = sqrt(0.5 * (mag - z.real) );
	if(z.imag < 0)
		dummy.imag = - dummy.imag;
	dummy.real = z.imag / dummy.imag /2;
	}
return(&dummy);
}



complex *c_conj(z)
complex z;
{
	dummy.real=z.real;
	dummy.imag=-z.imag;
	return(&dummy);
}
