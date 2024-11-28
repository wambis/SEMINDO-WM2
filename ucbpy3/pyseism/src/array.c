#define U (unsigned)
#include <stdio.h>
#include <stdlib.h>
#include <messageiolib.h>
#include <cmplx.h>

/*======================= char version ==============================*/

/*******1-D*******/
char *carray1(n11,n12)
int n11,n12;
{
char *m;
	m=(char *)malloc(U(n12-n11+1));
	if(!m) stop("allocation error in carray1");
	return(m-n11);
}
void free_carray1(a,n11)
char *a;
int n11;
{
	free(&a[n11]);
}
/*******2-D*******/
char **carray2(n11,n12,n21,n22)
int n11,n12,n21,n22;
{
int i;
char **m;
	m=(char **)malloc( U (n12-n11+1)*sizeof(char *) );
	if(!m) stop("allocation error 1 in carray2");
	m-=n11;
	
	for(i=n11;i<=n12;i++)
	{ m[i]=(char *)malloc( U (n22-n21+1)*sizeof(char) );
	  if(!m[i]) stop("allocation error 2 in carray2");
	  m[i]-=n21;
	}
	return(m);
}
void free_carray2(a,n11,n12,n21)
char **a;
int n11,n12,n21;
{
int i;	
	for(i=n11;i<=n12;i++)
	  free(&(a[i][n21]));
	free(&(a[n11]));

}

/*======================= int version ==============================*/

/*******1-D*******/
int *iarray1(n11,n12)
int n11,n12;
{
int *m;
	m=(int *)malloc(U(n12-n11+1)*sizeof(int));
	if(!m) stop("allocation error in iarray1");
	return(m-n11);
}
void free_iarray1(a,n11)
int *a,n11;
{
	free(&a[n11]);
}

/*******2-D*******/
int **iarray2(n11,n12,n21,n22)
int n11,n12,n21,n22;
{
int i, **m;
	m=(int **)malloc( U (n12-n11+1)*sizeof(int*) );
	if(!m) stop("allocation error 1 in iarray2");
	m-=n11;
	
	for(i=n11;i<=n12;i++)
	{ m[i]=(int *)malloc( U (n22-n21+1)*sizeof(int) );
	  if(!m[i]) stop("allocation error 2 in iarray2");
	  m[i]-=n21;
	}
	return(m);
}
void free_iarray2(a,n11,n12,n21)
int **a,n11,n12,n21;
{
int i;	
	for(i=n11;i<=n12;i++)
	  free(&(a[i][n21]));
	free(&(a[n11]));

}
/*======================= float version ==============================*/

/*******1-D*******/
float *farray1(n11,n12)
int n11,n12;
{
float *m;
	m=(float *)malloc( U (n12-n11+1)*sizeof(float) );
	if(!m) stop("allocation error in farray1");
	return(m-n11);
}
void free_farray1(a,n11)
float *a;
int n11;
{
	free(&a[n11]);
}

/*******2-D*******/
float **farray2(n11,n12,n21,n22)
int n11,n12,n21,n22;
{
int i;
float **m;
	m=(float **)malloc( U (n12-n11+1)*sizeof(float*) );
	if(!m) stop("allocation error 1 in farray2");
	m-=n11;
	
	for(i=n11;i<=n12;i++)
	{ m[i]=(float *)malloc( U (n22-n21+1)*sizeof(float) );
	  if(!m[i]) stop("allocation error 2 in farray2");
	  m[i]-=n21;
	}
	return(m);
}
void free_farray2(a,n11,n12,n21)
float **a;
int n11,n12,n21;
{
int i;	
    for(i=n11;i<=n12;i++) {
	/* Free second dimension array (which contains data). */
        free(&(a[i][n21]));
    }
    /* Free first dimension array which pointed to second dimension. */
    free(&(a[n11]));
}

/*******3-D*******/
float ***farray3(n11,n12,n21,n22,n31,n32)
int n11,n12,n21,n22,n31,n32;
{ 
int i,j;
float ***a;
	a=(float ***)malloc( U (n12-n11+1)*sizeof(float**) );
	if(!a) stop("allocation error 1 in farray3");
	a-=n11;
	
	for(i=n11;i<=n12;i++)
	{ a[i]=(float **)malloc( U(n22-n21+1)*sizeof(float*) );
	  if(!a[i]) stop("allocation error 2 in farray3");
	  a[i]-=n21;
	  
	  for(j=n21;j<=n22;j++)
	  { a[i][j]=(float *)malloc(U(n32-n31+1)*sizeof(float));
	    if(!a[i][j]) stop("allocation error 3 in farray3");
	    a[i][j]-=n31;
	  }
	}
	return(a);
}
void free_farray3(a,n11,n12,n21,n22,n31)
float ***a;
int n11,n12,n21,n22,n31;
{
int i,j;
    for(i=n11;i<=n12;i++) {
        for(j=n21;j<=n22;j++) {
            /* Free third dimension (arrays containing data). */
            free(&(a[i][j][n31])); 
        }
        /* Free second dimension array which pointed to third dimension. */
        free(&(a[i][n21]));
    }
    /* Free first dimension array which pointed to second dimension. */
    free(&(a[n11]));
}

/*******4-D*******/
float ****farray4(n11,n12,n21,n22,n31,n32,n41,n42)
int n11,n12,n21,n22,n31,n32,n41,n42;
{ 
int i,j,k;
float ****a;
	a=(float ****)malloc( U (n12-n11+1)*sizeof(float***) );
	if(!a) stop("allocation error 1 in farray4");
	a-=n11;
	
	for(i=n11;i<=n12;i++)
	{ a[i]=(float ***)malloc( U(n22-n21+1)*sizeof(float**) );
	  if(!a[i]) stop("allocation error 2 in farray4");
	  a[i]-=n21;
	  
	  for(j=n21;j<=n22;j++)
	  { a[i][j]=(float **)malloc(U(n32-n31+1)*sizeof(float*));
	    if(!a[i][j]) stop("allocation error 3 in farray4");
	    a[i][j]-=n31;
	    
	    for(k=n31;k<=n32;k++)
	    { a[i][j][k]=(float *)malloc(U(n42-n41+1)*sizeof(float));
	      if(!a[i][j][k]) stop("allocation error 4 in farray4");
	      a[i][j][k]-=n41;
	    }
	  }
	}
	return(a);
}
void free_farray4(a,n11,n12,n21,n22,n31,n32,n41)
float ****a;
int n11,n12,n21,n22,n31,n32,n41;
{
int i,j,k;	
	for(i=n11;i<=n12;i++)
	{  for(j=n21;j<=n22;j++)
	  { for(k=n31;k<=n32;k++) free(&a[i][j][k][n41]);
	 free(&a[i][j][n31]);
	  }
	 free(&a[i][n21]);
	}
	free(a[n11]);

}
/*======================= short version ==============================*/

/*******1-D*******/
short *sarray1(n11,n12)
int n11,n12;
{
short *m;
	m=(short *)malloc( U (n12-n11+1)*sizeof(short) );
	if(!m) stop("allocation error in sarray1");
	return(m-n11);
}
void free_sarray1(a,n11)
short *a;
int n11;
{
	free(&a[n11]);
}
/*======================= complex version ===========================*/

/*******1-D*******/
complex *cmplxarray1(n11,n12)
int n11,n12;
{
complex *m;
	m=(complex *)malloc( U (n12-n11+1)*sizeof(complex) );
	if(!m) stop("allocation error in cmplxarray1");
	return(m-n11);
}
void free_cmplxarray1(a,n11)
complex *a;
int n11;
{
	free(&a[n11]);
}
/*******2-D*******/
complex **cmplxarray2(n11,n12,n21,n22)
int n11,n12,n21,n22;
{
int i;
complex **m;
	m=(complex **)malloc( U (n12-n11+1)*sizeof(complex*) );
	if(!m) stop("allocation error 1 in cmplxarray2");
	m-=n11;
	
	for(i=n11;i<=n12;i++)
	{ m[i]=(complex *)malloc( U (n22-n21+1)*sizeof(complex) );
	  if(!m[i]) stop("allocation error 2 in cmplxarray2");
	  m[i]-=n21;
	}
	return(m);
}
void free_cmplxarray2(a,n11,n12,n21)
complex **a;
int n11,n12,n21;
{
int i;	
	for(i=n11;i<=n12;i++)
	  free(&(a[i][n21]));
	free(&(a[n11]));
}
/*******3-D*******/
complex ***cmplxarray3(n11,n12,n21,n22,n31,n32)
int n11,n12,n21,n22,n31,n32;
{
int i,j;
complex ***m;
	m=(complex ***)malloc( U (n12-n11+1)*sizeof(complex**) );
	if(!m) stop("allocation error 1 in cmplxarray3");
	m-=n11;
	
	for(i=n11;i<=n12;i++)
	{ m[i]=(complex **)malloc( U (n22-n21+1)*sizeof(complex*) );
	  if(!m[i]) stop("allocation error 2 in cmplxarray3");
	  m[i]-=n21;
	  
	  for(j=n21;j<=n22;j++)
	  { m[i][j]=(complex*)malloc( U (n32-n31+1)*sizeof(complex));
	    if(!m[i][j]) stop("allocation error 3 in cmplxarray3");
	    m[i][j]-=n31;
	  }
	}
	return(m);
}
void free_cmplxarray3(a,n11,n12,n21,n22,n31)
complex ***a;
int n11,n12,n21,n22,n31;
{
int i,j;	
	for(i=n11;i<=n12;i++)
	{  for(j=n21;j<=n22;j++) free(&(a[i][j][n31]));
	free(&(a[i][n21]));
	}
	free(&(a[n11]));

}

/*======================= double version ===========================*/

/*******1-D*******/
double *darray1(n11,n12)
int n11,n12;
{
double *m;
	m=(double *)malloc( U (n12-n11+1)*sizeof(double) );
	if(!m) stop("allocation error in darray1");
	return(m-n11);
}
void free_darray1(a,n11)
double *a;
int n11;
{
	free(&a[n11]);
}
