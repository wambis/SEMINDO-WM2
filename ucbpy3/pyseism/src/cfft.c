/*
     (1)    *********** cfour(data,n,isign)  **************     
            (cifour is a version for integer data)  

	input to cfour is an array of complex values either defined as a 
	complex structure or a real array with real and imaginary parts 
	alternating.  cfour makes no assumptions about the data. isign 
	is either + or - 1 depending on whether you are taking a 
	transform or inverse transform. To get the correct amplitudes 
	back from cfour after inverse transforming divide all values by 
	the number of data points n. 	

        discrete Fourier transform functions based on the fast Fourier 
        transform algorithm which computes the quantity :
              n-1
         x(j)=sum  x(k) exp( 2*pi*i*isign*k*j/n )
              k=0
       n must be a integral power of 2. beware of roundoff and overflow 
       error if you use cifour.
       

     (2)    ************cfftr(data,n)************            
            (cifftr is a version for integer data)
            
	    ***********cfftri(data,m)************
	    (cifftri is a version for integer data)
	
	cfftr assumes the input data is pure real. To get the correct 
	amplitudes back from cfftri divide all values by (n/2).
        computing discrete Fourier transform. and the transform has 
        n/2+1 complex values starting with frequency zero and ending at 
        the Nyquist frequency. parameters: n ... the number of samples, 
        an integral power of 2

        cfftri computes the inverse fast Fourier transform of a real 
        time series. m/2 + 1 complex frequency values are input and m 
        real timeseries are returned, where m is a power of two. the 
        complex frequencies are stores in the array x, with real and
        imaginary parts alternating. To get correct amplitude, divide 
        the output by n/2 (where n has been used in cfftr).
	
     (3)  **********discrete_fft(data,n,isign)*********
           (real_dis_fft takes a array of real data)

        subroutine discrete_fft constructs the discrete fast Fourier 
        transform of a time series using the Cooley-Tukey mixed Radix 
        method.  Algorithum is modified from Conte and de Boor, 
        Elementary Numerical Analysis, page 281.

        Input:
 	z1=pointer to complex array of data to be Fourier transformed
                Length of array should be n+1
 	n=number of data points in data
 	isign=sign of transform; 1 for fft, -1 for inverse fft
 
        Output:
 	z1=pointer to complex array Fourier transformed data
 		note: input data array is over written
 		Also if a data set is transformed and then
 		inverse transformed, output amplitudes should
 		be divided by the n
               ----------------------------------------
        real_dis_fft(data,npts,isign) computes the discrete Fourier 
        transform of  a real input data array. 
 
        Inputs:
 		data=real array of data
 			Note: must be dimensioned to 2*(npts+1)
 		npts=number of data points in data
 		isign=sign of Fourier transform 1 for transform -1 for 
 		      inverse
 
        Outputs:
 		data=array of transformed data
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmplx.h>

cfour( data, n, isign )
int n, isign;
float data[];
{
	int ip0, ip1, ip2, ip3, i3rev;
	int i1, i2a, i2b, i3;
	float sinth, wstpr, wstpi, wr, wi, tempr, tempi, theta;
	double sin();
	ip0=2;
	ip3=ip0*n;
	i3rev=1;
	for( i3=1; i3<=ip3; i3+=ip0 ) {
		if( i3 < i3rev ) {
			tempr = data[i3-1];
			tempi = data[i3];
			data[i3-1] = data[i3rev-1];
			data[i3] = data[i3rev];
			data[i3rev-1] = tempr;
			data[i3rev] = tempi;
			}
		ip1 = ip3 / 2;
		do {
			if( i3rev <= ip1 )
				break;
			i3rev -= ip1;
			ip1 /= 2;
			}
			while ( ip1 >= ip0 );
		i3rev += ip1;
		}
		ip1 = ip0;
		while ( ip1 < ip3 ) {
			ip2 = ip1 * 2;
			theta = 6.283185/( (float) (isign*ip2/ip0) );
			sinth = (float) sin( (double) (theta/2.) );
			wstpr = -2.*sinth*sinth;
			wstpi = (float) sin( (double) theta );
			wr = 1.;
			wi = 0.;
			for ( i1=1; i1<=ip1; i1+=ip0 ) {
			for ( i3=i1; i3<ip3; i3+=ip2 ) {
				i2a=i3;
				i2b=i2a+ip1;
				tempr = wr*data[i2b-1] - wi*data[i2b];
				tempi = wr*data[i2b] + wi*data[i2b-1];
				data[i2b-1] = data[i2a-1] - tempr;
				data[i2b] = data[i2a] - tempi;
				data[i2a-1] += tempr;
				data[i2a] += tempi;
				}
				tempr = wr;
				wr = wr*wstpr - wi*wstpi + wr;
				wi = wi*wstpr + tempr*wstpi + wi;
				}
			ip1=ip2;
			}
		return;
		}
cifour( data, n, isign )
int data[], n, isign;
{
	int ip0, ip1, ip2, ip3, i3rev, itempr, itempi;
	int i1, i2a, i2b, i3;
	float sinth, wstpr, wstpi, wr, wi, tempr, tempi, theta;
	double sin();
	ip0=2;
	ip3=ip0*n;
	i3rev=1;
	for( i3=1; i3<=ip3; i3+=ip0 ) {
		if( i3 < i3rev ) {
			itempr = data[i3-1];
			itempi = data[i3];
			data[i3-1] = data[i3rev-1];
			data[i3] = data[i3rev];
			data[i3rev-1] = itempr;
			data[i3rev] = itempi;
			}
		ip1 = ip3 / 2;
		do {
			if( i3rev <= ip1 )
				break;
			i3rev -= ip1;
			ip1 /= 2;
			}
			while ( ip1 >= ip0 );
		i3rev += ip1;
		}
		ip1 = ip0;
		while ( ip1 < ip3 ) {
			ip2 = ip1 * 2;
			theta = 6.283185/( (float) (isign*ip2/ip0) );
			sinth = (float) sin( (double) (theta/2.) );
			wstpr = -2.*sinth*sinth;
			wstpi = (float) sin( (double) theta );
			wr = 1.;
			wi = 0.;
			for ( i1=1; i1<=ip1; i1+=ip0 ) {
			for ( i3=i1; i3<ip3; i3+=ip2 ) {
				i2a=i3;
				i2b=i2a+ip1;
				itempr = ( (int) (wr*( (float) data[i2b-1] ) - wi*( (float) data[i2b] )));
				itempi = ( (int) (wr*( (float) data[i2b] ) + wi*( (float) data[i2b-1] )));
				data[i2b-1] = data[i2a-1] - itempr;
				data[i2b] = data[i2a] - itempi;
				data[i2a-1] += itempr;
				data[i2a] += itempi;
				}
				tempr = wr;
				wr = wr*wstpr - wi*wstpi + wr;
				wi = wi*wstpr + tempr*wstpi + wi;
				}
			ip1=ip2;
			}
		return;
		}

cfftr( x, n )
int n;
float x[];
{
	int nn, is, nm, j, i;
	int k1j, k1i, k2j, k2i;
	float s, fn, ex, wr, wi, wwr, wrr, wwi, a1, a2, b1, b2;
	double sin(), cos();
	nn = n/2;
	is = 1;
	cfour( x, nn, is );
	nm = nn/2;
	s = x[0];
	x[0] += x[1];
	x[n] = s - x[1];
	x[1] = 0.0 ;
	x[n+1] = 0.0;
	x[nn+1] = (-x[nn+1]);
	fn = (float) n;
	ex = 6.2831853 / fn;
	j = nn;
	wr = 1.0;
	wi = 0.0;
	wwr = (float) cos( (double) ex );
	wwi = (float) (-sin( (double) ex ));
	for (i=2; i<=nm; i++) {
		wrr = wr*wwr-wi*wwi;
		wi = wr*wwi+wi*wwr;
		wr = wrr;
		k1j = 2*j-1;
		k1i = 2*i-1;
		k2j = 2*j;
		k2i = 2*i;
		a1 = 0.5*(x[k1i-1]+x[k1j-1]);
		a2 = 0.5*(x[k2i-1]-x[k2j-1]);
		b1 = 0.5*(-x[k1i-1]+x[k1j-1]);
		b2 = 0.5*(-x[k2i-1]-x[k2j-1]);
		s = b1;
		b1 = b1*wr+b2*wi;
		b2 = b2*wr-s*wi;
		x[k1i-1] = a1-b2;
		x[k2i-1] = (-a2-b1);
		x[k1j-1] = a1+b2;
		x[k2j-1] =  a2-b1;
		j -= 1;
	}
	return;
}
cifftr( x, n )
int x[], n;
{
	int nn, is, nm, j, i, sv;
	int k1j, k1i, k2j, k2i;
	float s, fn, ex, wr, wi, wwr, wrr, wwi, a1, a2, b1, b2;
	double sin(), cos();
	nn = n/2;
	is = 1;
	cifour( x, nn, is );
	nm = nn/2;
	sv = x[0];
	x[0] += x[1];
	x[n] = sv - x[1];
	x[1] = 0.0 ;
	x[n+1] = 0.0;
	x[nn+1] = (-x[nn+1]);
	fn = (float) n;
	ex = 6.2831853 / fn;
	j = nn;
	wr = 1.0;
	wi = 0.0;
	wwr = (float) cos( (double) ex );
	wwi = (float) (-sin( (double) ex ));
	for (i=2; i<=nm; i++) {
		wrr = wr*wwr-wi*wwi;
		wi = wr*wwi+wi*wwr;
		wr = wrr;
		k1j = 2*j-1;
		k1i = 2*i-1;
		k2j = 2*j;
		k2i = 2*i;
		a1 = 0.5* (float)(x[k1i-1]+x[k1j-1]);
		a2 = 0.5* (float)(x[k2i-1]-x[k2j-1]);
		b1 = 0.5* (float)(-x[k1i-1]+x[k1j-1]);
		b2 = 0.5* (float)(-x[k2i-1]-x[k2j-1]);
		s = b1;
		b1 = b1*wr+b2*wi;
		b2 = b2*wr-s*wi;
		x[k1i-1] = (int) (a1-b2);
		x[k2i-1] = (int) (-a2-b1);
		x[k1j-1] = (int) (a1+b2);
		x[k2j-1] = (int) (a2-b1);
		j -= 1;
	}
	return;
}
cfftri( x, n )
int n;
float x[];
{
	int nn, is, nm, j, i, k1j, k1i, k2j, k2i;
	float s, fn, ex, wr, wi, wwr, wwi, wrr, a1, a2, b1, b2;
	double sin(), cos();
	nn = n/2;
	s = x[0];
	x[0] = 0.5 * ( x[0] + x[n] );
	x[1] = 0.5 * ( s - x[n] );
	x[nn+1] = (-x[nn+1]);
	is = -1;
	nm = nn/2;
	fn = (float) n ;
	ex = 6.2831853 / fn ;
	j = nn;
	wr = 1.0;
	wi = 0.0;
	wwr = (float) cos ( (double) ex );
	wwi = (float) ( - sin( (double) ex ) );
	for ( i=2; i<=nm; i++ ) {
		wrr = wr*wwr-wi*wwi;
		wi = wr*wwi+wi*wwr;
		wr = wrr;
		k1j = 2*j-1;
		k1i = 2*i-1;
		k2j = 2*j;
		k2i = 2*i;
		a1 = 0.5 * ( x[k1i-1] + x[k1j-1] );
		a2 = 0.5 * ( x[k2i-1] - x[k2j-1] );
		b1 = 0.5 * (-x[k1i-1] + x[k1j-1] );
		b2 = 0.5 * (-x[k2i-1] - x[k2j-1] );
		s = b1;
		b1 = b1*wr+b2*wi;
		b2 = b2*wr-s*wi;
		x[k1i-1] = (a1 - b2);
		x[k2i-1] = (-a2-b1);
		x[k1j-1] = (a1+b2);
		x[k2j-1] = (a2-b1);
		j -= 1;
		}
	cfour(x, nn, is);
	return;
	}
cifftri( x, n )
int x[], n;
{
	int nn, is, nm, j, i, k1j, k1i, k2j, k2i, sv;
	float s, fn, ex, wr, wi, wwr, wwi, wrr, a1, a2, b1, b2;
	double sin(), cos();
	nn = n/2;
	sv = x[0];
	x[0] = 0.5 * ( x[0] + x[n] );
	x[1] = 0.5 * ( sv - x[n] );
	x[nn+1] = (-x[nn+1]);
	is = -1;
	nm = nn/2;
	fn = (float) n ;
	ex = 6.2831853 / fn ;
	j = nn;
	wr = 1.0;
	wi = 0.0;
	wwr = (float) cos ( (double) ex );
	wwi = (float) ( - sin( (double) ex ) );
	for ( i=2; i<=nm; i++ ) {
		wrr = wr*wwr-wi*wwi;
		wi = wr*wwi+wi*wwr;
		wr = wrr;
		k1j = 2*j-1;
		k1i = 2*i-1;
		k2j = 2*j;
		k2i = 2*i;
		a1 = 0.5 * (float) ( x[k1i-1] + x[k1j-1] );
		a2 = 0.5 * (float) ( x[k2i-1] - x[k2j-1] );
		b1 = 0.5 * (float) (-x[k1i-1] + x[k1j-1] );
		b2 = 0.5 * (float) (-x[k2i-1] - x[k2j-1] );
		s = b1;
		b1 = b1*wr+b2*wi;
		b2 = b2*wr-s*wi;
		x[k1i-1] = (int) (a1 - b2);
		x[k2i-1] = (int) (-a2-b1);
		x[k1j-1] = (int) (a1+b2);
		x[k2j-1] = (int) (a2-b1);
		j -= 1;
		}
	cifour(x, nn, is);
	return;
	}

discrete_fft(z1,n,isign)

	int n,isign;
	complex *z1;
{

	complex *z2;   /* temporary work space */
	int after,now,before,next,nextmx,inz,i;	
	static int prime[12]={2,3,5,7,11,13,17,19,23,29,31,37};

	after=1;
	before=n;
	next=0;
	nextmx=12;
	inz=1; /* flip-flop variable to keep track of transform */

/* open up temporary work space */

	if((z2=(complex *)calloc((unsigned)(n*2+1),sizeof(float))) == NULL)
	{
		fprintf(stderr,"Error allocating space for fft\n");
		exit(-2);
	}

	while(before != 1)
	{
nextprime:	if((before/prime[next])*prime[next] < before) /* find smallest prime */
		{
			++next;
			if(next < nextmx) goto nextprime; /* I know - I know its ugly */
			now=before;  /* set up variables for transform */
			before=1;
		}
		else
		{
			now=prime[next]; /* set up variables for transform */
			before=before/prime[next];
		}
		if(inz == 1) /* perform one step of transform */
		{
			fftstp(z1,after,now,before,z2,isign);
		}
		else
		{
			fftstp(z2,after,now,before,z1,isign);
		}
		inz=3-(inz); /* keep track of output */
		after*=now;
	}
	
	if(inz == 2) bcopy((char *)z2,(char *)z1,sizeof(complex)*n);  /* output is in z2 move it to z1 */
	cfree(z2); /* free temporary space */
}

/* subroutine carries out one step of the discrete Fourier transform 
 */

fftstp(zin,after,now,before,zout,isign)

	int after,before,now,isign;
	complex *zin,*zout;

{

	int ia,ib,in,j,inter_fact2,inter_fact3;
	complex arg,omega,value,temp,*pointer1,*pointer2,*array_end2;
	complex *array_end1,*pointer3;
	float angle,twopi=6.283185307;

/* pointer1 & pointer2 -> zin ; pointer3 -> zout */

	angle=twopi/(float)(isign*now*after);
	omega.real=cos(angle);
	omega.imag=(-1.0*sin(angle));
	arg.real=1.0;
	arg.imag=0.0;

	inter_fact2=after*before;  /* array increment for zin */
	inter_fact3=now*after;  /* array increment for zout */

	for(j=1;j<=now;j++)
	{
		for(ia=1;ia<=after;ia++)
		{

			array_end1=zin+(ia-1+after*(before-1+before*(now-1))); /* loop end criteria for zin */
			for(ib=1,pointer3=zout+ia-1+after*(j-1),pointer1=zin+ia-1+after*before*(now-1);pointer1<=array_end1;ib++,pointer3+=inter_fact3,pointer1+=after)
			{
				value.real=pointer1->real;
				value.imag=pointer1->imag;

				array_end2=zin+(ia-1+after*(ib-1)); /* loop end criteria for zin */
				for(pointer2=zin+(ia-1+after*(ib-1+before*(now-2)));pointer2>=array_end2;pointer2-=inter_fact2)
				{
					temp.real=value.real*arg.real-value.imag*arg.imag;
					temp.imag=value.real*arg.imag+value.imag*arg.real;
					value.real=temp.real+pointer2->real;
					value.imag=temp.imag+pointer2->imag;
				}
				pointer3->real=value.real;
				pointer3->imag=value.imag;
			}
			temp.real=arg.real*omega.real-arg.imag*omega.imag;
			arg.imag=arg.real*omega.imag+arg.imag*omega.real;
			arg.real=temp.real;
		}
	}
}
int real_dis_fft(data,npts,isign)
	float *data;
	int npts,isign;
{
	int i;

	if(isign == 1)  /* convert real array to complex array with zero imag */
	{
		for(i=npts;i>=0;--i)
		{
			data[i*2]=data[i];
			data[i*2+1]=0.0;
		}
	}
	else /* complete complex array using inherent symmetry */
	{
		for(i=0;i<npts;i+=2)
		{
			data[2*(npts)-i]=data[i];
			data[2*(npts)-i+1]=(-data[i+1]);
		}
	}
	discrete_fft((complex *)data,npts,isign);

	if(isign == -1) /* collapse complex array to a real array */
	{
		for(i=0;i<=npts;i++)
		{
			data[i]=data[i*2];
		}
	}
}

/******************************************************************************
 *                            FORTRAN PROTOCOLS 
 ******************************************************************************
                          All integers are integer*2 
 */


four_( data, n, isign )
int data[], *n, *isign;
{
	cfour( data, *n, *isign );
	return;
}
ifour_( data, n, isign )
int *n, *isign;
float data[];
{
	cifour( data, *n, *isign );
	return;
}

fftr_( x, n )
int *n;
float x[];
{
	cfftr( x, *n );
	return;
}
ifftr_( x, n )
int x[], *n ;
{
	cifftr( x, *n );
	return;
}

fftri_( x, n )
int *n;
float x[];
{
	cfftri( x, *n );
	return;
}
ifftri_( x, n )
int x[], *n ;
{
	cifftri( x, *n );
	return;
}
