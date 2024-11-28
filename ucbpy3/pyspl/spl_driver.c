#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "bspl.h"
#include "sspl.h"

#define min(a,b) ((a)>(b)?(b):(a))

int
_bspl_driver(
	int deriv, 
	double *knots, int nknots, 
	double *r,     int nr, 
	int k, 
	double *res,   int nres)
{
	int ir = 0;
	double (*f[])(int, int, double *, double) = {spl, d_spl, d2_spl};

	/* check deriv for supported value */
	if (deriv < 0 || deriv > 2) {
		printf("Error [%s]: b-spline derivative of order %i is unsupported\n", __func__, deriv);
		return 0;
	}

	/* loop over radii */
	for (ir = 0; ir < min(nr,nres); ir++)
		res[ir] = (*f[deriv])(k, nknots, knots, r[ir]);

	if (ir < nr)
		printf("Warning [%s]: result array too short - only evaluated first %i radii\n", __func__, ir);

	return ir;
}

int
_sspl_driver(
	int deriv, 
	double *knot_lons,  int nknot_lons, 
	double *knot_lats,  int nknot_lats, 
	double *knot_dists, int nknot_dists, 
	double *lons,       int nlons,
	double *lats,       int nlats, 
	double *resval,     int nresval, 
	int *resind,        int nresind,
	int *resloc,        int nresloc)
{
	int k, l, maxres, nknots, ir = 0;
	double *del;
	double (*f[])(double, double) = {spbsp, gspbsp};

	/* get max: (1) allowable storage for result; (2) number of spline knots */
	maxres = min(nresval, min(nresind, nresloc));
	nknots = min(nknot_dists, min(nknot_lons, nknot_lats));

	/* check deriv for sane value */
	if (deriv < 0 || deriv > 1) {
		printf("Error [%s]: spherical spline derivative of order %i is unsupported\n", __func__, deriv);
		return 0;
	}

	/* allocate storage for knot-distance array */
	if ((del = malloc((size_t)(nknots) * sizeof(*del))) == NULL) {
		printf("Error [%s]: cannot allocate memory for del - %s\n", __func__, strerror(errno));
		return 0;
	}

	/* loop over locations */
	for (l = 0; l < nlons; l++) {
		/* get distance to knots */
		getdel_vec(knot_lats, knot_lons, nknots, lats[l], lons[l], del);
		/* loop over knots */
		for (k = 0; k < nknots; k++) {
			if (del[k] < 2 * knot_dists[k]) {
				if (ir == maxres) {
					printf("Warning [%s]: result arrays are too short - only evaluated first %i non-zero knots\n", __func__, ir);
					return ir;
				}
				resval[ir] = (*f[deriv])(del[k], knot_dists[k]);
				resind[ir] = k;
				resloc[ir] = l;
				ir += 1;
			}
		}
	}

	free(del);

	return ir;
}
