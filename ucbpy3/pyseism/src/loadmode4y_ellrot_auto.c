/*calculating receiver and source vectors and store their productions*/
#include <stdio.h>
#include <string.h>

#include "messageiolib.h"
#include "stdlib.h"
#include "syn_struct.h"
#include "dimensiony_prem1002_10s.h"
#include "dimensionsynth.h"
#include "premy_st.h"
#include "cmplx.h"

extern int geteigys_ellrot();

void loadmode4y_ellrot_auto(modetype, whichbrnch, src, band, nmode, ll, ww, alpha, ell, grv, svect, rvect)
char modetype;			/* 'S' or 'T' */
char whichbrnch;		/* 'a'= all; 'f'= fundamentals; 'o'= overtunes */
source_st *src;			/* source parameters */
band_st *band;			/* information for band pass */
int *nmode;			/* number of modes loaded */
short *ll;			/* degree vector */
complex **svect;		/* svect[mode][m] */
complex **rvect;		/* rvect[mode][m] */
float *ww;			/* angular frequency vector */
float *alpha;			/* half-q vector */
float *ell;			/* ellipticity correction vector */
float *grv;			/* group velocity vector */
{
	int l, n, l0, minn, maxn, exist;
	sprem_st prem;

	if (modetype == 'S')
		l0 = 0;
	else if (modetype == 'T')
		l0 = 1;
	else
		stop("loadmode4: illegal modetype");

	switch (whichbrnch) {
	case 'a':
		minn = 0;
		maxn = MAXN;
		break;
	case 'f':
		minn = 0;
		maxn = 0;
		break;
	case 'o':
		minn = 1;
		maxn = MAXN;
		break;
	default:
		stop("loadmode4: illegal whichbrnch");
	}

	*nmode = 0;
	for (l = l0; l < MAXL; l++) {

		for (n = minn; n <= maxn; n++) {

			exist = geteigys_ellrot(n, modetype, l, src->dep, &prem);

			if (n == 0 && l == 1)
				exist = 0;
			if (exist == 0) {
				if (n == minn && l != 1)
					goto finished;
				if (n > minn && l != 1)
					break;
			}
			if (exist && prem.w > band->w4) {
				if (n == minn && l != 1)
					goto finished;
				if (n > minn)
					break;
			}

			if (exist && prem.w > band->w1) {
				ll[*nmode] = (short)l;

				ww[*nmode] = prem.w;

				alpha[*nmode] = 0.5 * prem.q * prem.w;

				ell[*nmode] = prem.ell;

				grv[*nmode] = prem.grv / 6371.;

				rvec(&prem, rvect[*nmode]);

				svec(src, &prem, svect[(*nmode)++]);

			}
		}
	}
 finished:;			// printf("\n%d modes loaded.",*nmode);
}
