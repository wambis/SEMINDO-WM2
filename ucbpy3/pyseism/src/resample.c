#define TPI 6.2831853
#define PI 3.14159265
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <messageiolib.h>
#include <struct.h>

int min2n();
float *farray1();

void resample(ndata,nfac,dt,u,trace)
int ndata,nfac;
float dt,*u;
syntrace_st *trace;
{
static int prvdim=0;
int n,nn,dim,i;
static float *work;

	n=min2n(ndata);
	nn=n*nfac;
	dim=nn+2;
	
	if(dim>prvdim)
	{ if(prvdim) free_farray1(work,0);
	  prvdim=dim;
	  work=farray1(0,dim);
	}
	
	for(i=0;i<ndata;i++) work[i]=u[i];
	for(i=ndata;i<n;i++) work[i]=0.;
	
	cfftr(work,n);
	
	for(i=0;i<n+2;i++) work[i]*=2./(float)n;
	for(i=n+2;i<nn+2;i++) work[i]=0.0;

	cfftri(work,nn);
	
	for(i=0;i<trace->ndata;i++) u[i]=work[i];

}
