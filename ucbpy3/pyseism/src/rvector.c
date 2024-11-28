/* calculate R_k^m in Woodhouse & Girnius */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cmplx.h"
#include "messageiolib.h"

float **farray2();

complex *rvector(theta, phi, l, rvec, compnt, mtype)
float theta, phi;		/* epicentral coordinates of reciever */
int l;
complex *rvec;
char compnt;			/* 'Z', 'L' or 'T' */
char mtype;
{
	static complex r[3], eip;	/* exp(i*phi) */
	static int prvl = -1;
	static float prvtheta = -1.0, prvphi = 2.1 * 3.1416;
	static float **p;
	int m;
	complex fac, temp;
	float **PNm();

	/*if(l!=prvl || theta!=prvtheta)
	   { if(l!=prvl)
	   { int NN,mm;
	   if(prvl!=-1) free_farray2(p,-1,1,-prvl);
	   p=farray2(-1,1,-l,l);
	   for(NN=-1;NN<=1;NN++)
	   for(mm=-l;mm<=l;mm++) p[NN][mm]=0.0;
	   prvl=l;
	   }
	   prvtheta=theta;
	   if(l) rotmatrix(1,l,theta,p);
	   else p[0][0]=1.0;
	   } */
	p = PNm(l, theta);

	if (phi != prvphi) {
		prvphi = phi;
		eip.real = (float)cos((double)phi);
		eip.imag = (float)sin((double)phi);
	}

	fac.real = 1.0;
	fac.imag = 0.0;

	if (compnt == 'Z') {
		if (mtype == 'T')
			stop("rvector: illegal combination");
		for (m = 0; m < 3; m++) {
			temp.real = p[0][m] * rvec[0].real;
			temp.imag = 0.;
			r[m] = *c_mult(temp, fac);
			fac = *c_mult(fac, eip);
		}
	} else {
		if ((mtype == 'S' && compnt == 'L') || (mtype == 'T' && compnt == 'T')) {
			for (m = 0; m < 3; m++) {
				temp.real = rvec[1].real * (-p[-1][m] + p[1][m]);
				temp.imag = 0.0;
				r[m] = *c_mult(temp, fac);
				fac = *c_mult(fac, eip);
			}
		} else {
			for (m = 0; m < 3; m++) {
				temp.real = 0.0;
				temp.imag = rvec[1].imag * (p[-1][m] + p[1][m]);
				r[m] = *c_mult(temp, fac);
				fac = *c_mult(fac, eip);
			}
		}
	}

	return (r);

}
