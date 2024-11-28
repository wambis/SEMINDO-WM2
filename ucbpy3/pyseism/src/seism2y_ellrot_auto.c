#define PI 3.14159265
#define HPI 1.570796
#define CONST 9.52278e-26	/* rn^{-4}*wn^{-2}*rhobar^{-1} */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ellatt.h"
#include "cmplx.h"
#include "syn_struct.h"
#include "dimensiony_prem1002_10s.h"
#include "dimensionsynth.h"
#include "messageiolib.h"

extern complex **cmplxarray2();
extern complex *rvector();
extern void loadmode4y_ellrot_auto();

void seism2y_ellrot_auto(flag, source, receiver, trace, mode, ifprint, no_ell, u_t)
int flag;			/* flag==1 --> event/band/modetype changed; */
int no_ell;
source_st *source;
recvlctn_st *receiver;
syntrace_st *trace;
int *mode;
int ifprint;
float *u_t;			/* seismogram */
{
	static int first = 1, dim = 0;
	static float ptheta = -1., pphi = 361.;
	static int nmode;
	static short ll[MAXMODE];
	static float ww[MAXMODE], aa[MAXMODE], ee[MAXMODE], grv[MAXMODE];
	static complex **rvect, **svect;
	static complex amp[MAXMODE];
	static float theta, phi;
	static float ellcor;
	complex rs;
	static band_st band;
	int i, nfac, ndata;
	float dt, phv, fac;
	double eta;

	if (trace->com != 'Z' && trace->com != 'L' && trace->com != 'T')
		stop("seism2: illegal component");
	if (first) {
		first = 0;
		rvect = cmplxarray2(0, MAXMODE - 1, 0, 1);
		svect = cmplxarray2(0, MAXMODE - 1, 0, 2, 0, 5);
	}

	if (flag == 1) {
		band.w1 = trace->whghcut;
		band.w2 = trace->whghcorner;
		band.w3 = trace->wlowcorner;
		band.w4 = trace->wlowcut;

		loadmode4y_ellrot_auto(trace->modetype, trace->whichbrnch, source,
				       &band, &nmode, ll, ww, aa, ee, grv, svect, rvect);

		if (nmode > MAXMODE)
			stop("seism2: not enough working space");
	}

	if (receiver->theta != ptheta || receiver->phi != pphi || flag == 1) {
		double thr, ths, dphi, a, temp;
		float coscth;
		complex *rr, sum, cc;
		ifprint = 0;
		if (ifprint)
			printf("\nCalculating amplitudes");

		ptheta = receiver->theta;
		pphi = receiver->phi;

		thr = receiver->theta;
		ths = source->theta;
		dphi = receiver->phi - source->phi;
		a = cos(ths) * cos(thr) + sin(ths) * sin(thr) * cos(dphi);
		if (a < -1.)
			a = -1.;	/* cure round-off error */
		if (a > 1.)
			a = 1.;

/* rotate to epicentral coordinates */
		theta = (float)acos(a);
		temp = sin(thr) * sin(dphi) / sin((double)theta);
		if (temp < -1.)
			temp = -1.;	/* cure round-off error */
		if (temp > 1.)
			temp = 1.;
		phi = (float)asin(temp);
		temp = cos(ths) * sin(thr) * cos(dphi) - sin(ths) * cos(thr);
		if (temp < 0.0)
			phi = PI - phi;

/* pole for ellipticity correction */
		a = a * a;
		temp = sin(ths) * sin(thr) * sin(dphi);
		if (a == 1.0)
			coscth = 0.;
		else
			coscth = temp * temp / (1. - a);
		if (coscth < -0.0001 || coscth > 1.0001)
			stop("seism2: error 1");
		ellcor = 1. - 3. * coscth;

		if (no_ell)
			ellcor = 0.0;

		/* product of source and receiver vector */
		for (i = 0; i < nmode; i++) {
			{	//if(!(i%1000)) if(ifprint) printf("\n  mode # %d",i);
				rr = rvector(theta, phi, (int)ll[i], rvect[i], trace->com, trace->modetype);
			}
			cc = *c_mult(svect[i][2], rr[2]);
			sum = *c_add(cc, *c_conj(cc));
			cc = *c_mult(svect[i][1], rr[1]);
			sum = *c_add(sum, *c_add(cc, *c_conj(cc)));
			cc = *c_mult(svect[i][0], rr[0]);
			rs = *c_add(sum, cc);

			/* amplitude */
			if (ww[i] < band.w2) {
				eta = PI * (ww[i] - band.w1) / (band.w2 - band.w1);
				fac = 0.5 * CONST * (1. - (float)cos(eta));
			} else if (ww[i] <= band.w3)
				fac = CONST;
			else {
				eta = PI * (ww[i] - band.w3) / (band.w4 - band.w3);
				fac = 0.5 * CONST * (1. + (float)cos(eta));
			}

			switch (trace->datatype) {
			case 'a':
				amp[i].real = fac * rs.real;
				amp[i].imag = fac * rs.imag;
				break;
			case 'v':
				amp[i].real = fac * rs.imag / ww[i];
				amp[i].imag = -fac * rs.real / ww[i];
				break;
			case 'd':
				amp[i].real = -fac * rs.real / ww[i] / ww[i];
				amp[i].imag = -fac * rs.imag / ww[i] / ww[i];
				break;
			default:
				stop("seism2: to be fixed");
			}
		}
	}

	dt = trace->dt;
	if (trace->ndata > 1000) {
		while (dt <= HPI / trace->wlowcut) {
			dt *= 2.;
		}
		nfac = (int)(dt / trace->dt);
		ndata = 1 + (trace->ndata - 1) / nfac;
	} else {
		nfac = 1;
		ndata = trace->ndata;
	}
	for (i = 0; i < ndata; i++)
		u_t[i] = 0.0;
	if (ifprint)
		printf("\nSuperposing modes");
	*mode = 0;
	for (i = 0; i < nmode; i++) {

		if (!ATTEN)
			aa[i] = 0.0;
		if (grv[i] >= trace->gv1 && grv[i] <= trace->gv2) {
			phv = ww[i] / (0.5 + (float)ll[i]);
			if (phv > trace->pv1 && phv <= trace->pv2) {
				fac = 1. + ellcor * ee[i];
				addmodes(ndata, u_t, trace->t0, dt, fac * ww[i], aa[i], &amp[i]);
				(*mode)++;
			}
		}
	}
	if (nfac > 1)
		resample(ndata, nfac, dt, u_t, trace);
}
