/* calculate R_kN in Woodhouse & Girnius. 
 * v_r=v_t=v_p=1 is assumed; S modes: v_t -> real part; v_p -> imag part
 *                           T modes: v_t -> imag part; v_p -> real part
 */

#define FPI 12.566371

#include <math.h>

#include "structprem.h"
#include "cmplx.h"

void rvec(pm, vec)
sprem_st *pm;
complex *vec;
{
	int l;
	float L, k0, k1;
	l = pm->l;
	L = (float)(l * (l + 1));

	k0 = -(float)sqrt((double)(2 * l + 1) / FPI);
	/* va and ha have opposite signs to 
	 * U[sea floor] and V[sea floor] 
	 */
	k1 = 0.5 * k0 * (float)sqrt((double)L);

	vec[0].real = k0 * pm->va;
	vec[0].imag = 0.;

	vec[1].real = vec[1].imag = -k1 * pm->ha;
	if (pm->typ == 'T')
		vec[1].real = -vec[1].real;
}
