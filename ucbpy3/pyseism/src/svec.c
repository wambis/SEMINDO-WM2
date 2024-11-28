/* calculate S_kN in Woodhouse & Girnius, but with a different sign */

#define FPI 12.566371

#include <math.h>

#include "cmplx.h"
#include "structprem.h"
#include "struct.h"

void svec(src, pm, vec)
source_st *src;
sprem_st *pm;
complex *vec;
{
	float r, h, hn, hn2, hn3;
	float a, b, v, vp, u, up, x, L, k0, k1, k2;
	int l;
	r = 1. - src->dep / 6371000.;
	h = r - pm->r1;
	hn = pm->r2 - pm->r1;
	hn2 = 1. / hn / hn;
	hn3 = hn2 / hn;

	a = hn3 * (hn * (pm->v1p + pm->v2p) + 2. * (pm->v1 - pm->v2));
	b = hn2 * (-hn * (pm->v2p + 2. * pm->v1p) + 3. * (pm->v2 - pm->v1));
	v = pm->v1 + h * (pm->v1p + h * (b + h * a));
	vp = pm->v1p + h * (2. * b + 3. * h * a);
	x = vp - v / r;

	if (pm->typ == 'S') {
		a = hn3 * (hn * (pm->u1p + pm->u2p) + 2. * (pm->u1 - pm->u2));
		b = hn2 * (-hn * (pm->u2p + 2. * pm->u1p) + 3. * (pm->u2 - pm->u1));
		u = pm->u1 + h * (pm->u1p + h * (b + h * a));
		up = pm->u1p + h * (2. * b + 3. * h * a);
		x += u / r;
	}

	l = pm->l;
	L = (float)(l * (l + 1));

	k0 = (float)sqrt((double)(2 * l + 1) / FPI);
	k1 = 0.5 * k0 * (float)sqrt((double)L);
	k2 = 0;
	if (l > 0)
		k2 = 0.5 * k1 * (float)sqrt((double)((l - 1) * (l + 2)));

	if (pm->typ == 'S') {
		vec[0].real = k0 * (up * src->moment[0] + (u - .5 * L * v) * (src->moment[1] + src->moment[2]) / r);
		vec[0].imag = 0.;
		vec[1].real = -k1 * x * src->moment[3];
		vec[1].imag = k1 * x * src->moment[4];
		vec[2].real = k2 * v * (src->moment[1] - src->moment[2]) / r;
		vec[2].imag = -2. * k2 * v * src->moment[5] / r;
	} else {
		vec[0].real = vec[0].imag = 0.;
		vec[1].real = k1 * x * src->moment[4];
		vec[1].imag = k1 * x * src->moment[3];
		vec[2].real = -2. * k2 * v * src->moment[5] / r;
		vec[2].imag = k2 * v * (src->moment[2] - src->moment[1]) / r;
	}
}
