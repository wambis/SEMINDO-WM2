/* This routine calculates the rotation matrix d_Nm^(l)(beta) 
 */

#define PI 3.14159265
#define D double
#define F float
#define BIG 1.0e34
#define LOGBIG 34.0
#define SMALL 1.0e-34

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "messageiolib.h"

void rotmatrix(Nmax, l, beta, d)
int Nmax, l;
float beta;
float **d;			/* 2*Nmax+1 by 2*l+1  matrix */
{
	int swtch, change;
	int N, m;
	float shbt, chbt, sbt, cbt, logf, logs;
	float fN, t1, sign, part, fac, temp;
	int mult, i;

	if (Nmax > l)
 STOP("rotmatrix:N larger than l");
	if (l == 0) {
		**d = 1.0;
		return;
	}

	change = 0;
	if (beta < 0) {
		change = 1;
		beta = -beta;
	}

	swtch = 1;
	if (beta > PI / 2.) {
		beta = PI - beta;
		swtch = -1;
	}

	if (beta < 1.e-04 && beta > -1.e-04)
		beta = 0.0;
	if (beta < 0.0)
 STOP("rotmatrix:illegal beta");
/*-------------------------------*/

	if (beta == 0.0) {
		for (N = -Nmax; N <= Nmax; N++) {
			sign = 1.0;
			for (m = -l; m <= l; m++) {
				d[N][m] = 0.0;
				if ((N * swtch) == m)
					d[N][m] = sign;
				sign *= swtch;
			}
		}
		return;
	}
/*-------------------------------*/
	for (N = -Nmax; N <= Nmax; N++)	/*zero l.h.s. of matrix */
		for (m = -l; m <= 0; m++)
			d[N][m] = 0.0;

	/*set up parameters */
	shbt = sin((D) (0.5 * beta));
	chbt = cos((D) (0.5 * beta));
	sbt = 2.0 * shbt * chbt;
	cbt = 2.0 * chbt * chbt - 1.0;
	logf = log10((D) (chbt / shbt));
	logs = log10((D) shbt);

/*iterate from last column using 1. as starting value*/
	for (N = -Nmax; N <= Nmax; N++) {
		fN = (F) N;
		d[N][l] = 1.0;
		for (m = l - 1; m >= abs(N); m--) {
			if (m + 1 == l)
				t1 = 0.0;
			else
				t1 = -d[N][m + 2] * sqrt((D) ((l + 2 + m) * (l - 1 - m)));
			d[N][m] = t1 - (2.0 / sbt) * (cbt * (F) (m + 1) - fN) * d[N][m + 1];
			d[N][m] /= sqrt((D) ((l + m + 1) * (l - m)));

			if (fabs((D) d[N][m]) >= BIG && m != abs(N)) {
				d[N][-m] = LOGBIG;
				d[N][m] /= BIG;
				d[N][m + 1] /= BIG;
			}
		}
	}

/*set up normalization for last column;
 *using the first column as temp working space to store the logrithms.
 */
	t1 = logs * (F) (2 * l);
	for (N = -l + 1; N <= -Nmax; N++)
		t1 += logf + 0.5 * log10((D) (l + 1 - N) / (D) (N + l));
	d[-Nmax][-l] = t1;
	for (N = -Nmax + 1; N <= Nmax; N++)
		d[N][-l] = d[N - 1][-l] + logf + 0.5 * log10((D) (l - N + 1) / (D) (l + N));
/*renormalize rows*/
	sign = 1.0;
	if ((l - Nmax) % 2)
		sign = -1.0;
	for (N = -Nmax; N <= Nmax; sign = -sign, N++) {
		part = d[N][-l];
		mult = 1;
		while (fabs((D) part) >= LOGBIG) {
			mult *= 2;
			part *= 0.5;
		}
		fac = pow(10.0, (D) part);
		for (m = l; m >= abs(N); m--) {
			if (d[N][-m + 1] != 0.0 && -m + 1 < -abs(N)) {
				part = part * (F) (mult) + d[N][-m + 1];
				mult = 1;
				while (fabs((D) part) >= LOGBIG) {
					mult *= 2;
					part *= 0.5;
				}
				fac = pow(10.0, (D) part);
			}
			for (i = 0; i < mult; i++) {
				d[N][m] *= fac;
				if (fabs((D) d[N][m]) <= SMALL)
					break;
			}
			d[N][m] *= sign;
		}
	}
/* correction for beta<0, using d[N][m](-b)=(-1)**(N+m)*d[N][m](b) */

	if (change) {
		for (N = -Nmax; N <= Nmax; N++) {
			sign = 1.0;
			if ((l + N) % 2)
				sign = -1.0;
			for (m = l; m >= abs(N); m--, sign = -sign)
				d[N][m] *= sign;
		}
	}

/*fill rest of matrix*/
	if (swtch < 0) {
		sign = 1.0;
		if ((l - Nmax) % 2)
			sign = -1.0;
		for (N = -Nmax; N <= Nmax; sign = -sign, N++)
			for (m = -l; m <= -abs(N); m++)
				d[N][m] = sign * d[N][-m];

		for (N = -Nmax; N <= Nmax; N++)
			for (sign = 1.0, m = abs(N); m <= l; sign = -sign, m++) {
				d[N][m] = sign * d[-N][-m];
				if (m <= Nmax) {
					d[-m][-N] = d[N][m];
					d[m][N] = d[-N][-m];
				}
			}
	}

	else {
		for (N = -Nmax; N <= Nmax; N++)
			for (sign = 1.0, m = -abs(N); m >= -l; sign = -sign, m--) {
				d[N][m] = sign * d[-N][-m];
				if (abs(m) <= Nmax) {
					d[-m][-N] = d[N][m];
					d[m][N] = d[-N][-m];
				}
			}
	}
}
