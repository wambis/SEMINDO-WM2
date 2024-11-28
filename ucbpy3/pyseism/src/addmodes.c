#include <math.h>

#include "cmplx.h"

void addmodes(ndata,u,t0,dt,w,alpha,amp)
int ndata;
float *u,t0,dt,w,alpha;
complex *amp;
{
int i;
float t,e,creal,cimag;
double wt;
	for(i=0;i<ndata;i++)
	{ t=t0+dt*(float)i;
	  e=(float)exp(-(double)(alpha*t));
	  wt=w*t;
	  creal=(float) cos(wt);
	  cimag=(float) sin(wt);
	  e=e*(creal*amp->real - cimag*amp->imag);
	  u[i]+=e;
	}
}
