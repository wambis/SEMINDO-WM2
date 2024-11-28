/*
 * b-splines on arbitrary spacing 
 * Modified from code by CHM 12/97
 */

#include <math.h>
#include <stdio.h>

#define TOL 1.5

void		fill_hh   ();

#define NKNOTS_MAX 1000

double 
spl(ord, nknots, knot, xi)
	int		ord       , nknots;
	double          *knot, xi;
/*
 * ord: number of rho(x) nknots : # of knkots = Nx+1 (Nx=index of highest
 * spline) xi: point of interest splines defined as : f_i(x) = a_i(x-x_i)^3 +
 * b_i(x-x_i)^2 + c_i(x-x_i) + d_i
 */
{
	int		ii        , Nx;
	/* double          *hh, rho_x; */
	double           rho_x;
	double           hh[NKNOTS_MAX];
	double		coefa   , coefb, coefc, coefd;

	Nx = nknots - 1;
	/* Compute vector hh of spacings */
	/* hh = farray1(0, Nx - 1); */
	fill_hh(hh, knot, Nx);

	/* Consistency checks */
	if ((xi - (double)TOL) > knot[Nx]) {
		/* printf("xi=%g / knot[%d]=%g", xi, Nx, knot[Nx]);
		stop("spl: xi>knot[Nx]"); */
		return 0.0;
	} else if ((xi + (double)TOL) < knot[0]) {
		/* printf("xi=%g / knot[0]=%g", xi, knot[0]);
		stop("spl: xi<knot[0]"); */
		return 0.0;
	} else if (ord > Nx) {
		/* stop("spl: order > Nx"); */
		printf("Warning [%s]: spl index %i exceeds knot count %i\n", __func__, ord, Nx);
		return 0.0;
	}

	if (ord == 0) {		/* LHS */
		double		denom;
		denom = 3. * hh[ord] * hh[ord] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1];
		if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 0.0;
			coefc = -12 / denom;
			coefd = 4 * (2 * hh[ord] + hh[ord + 1]) / denom;

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		} else if (xi > knot[ord + 1] & xi <= knot[ord + 2]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord + 1] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 12 / ((hh[ord] + hh[ord + 1]) * denom);
			coefc = -12. * hh[ord + 1] / ((hh[ord] + hh[ord + 1]) * denom);
			coefd = 4. * hh[ord + 1] * hh[ord + 1] / ((hh[ord] + hh[ord + 1]) * denom);

			rho_x = coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefb * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefc * (xi - knot[ord + 1]);
			rho_x += coefd;
		} else		/* x>x2 */
			rho_x = 0.0;
	} else if (ord == 1) {	/* LHS+1 */
		double		denom   , denomsum, dd;
		denom = (3. * hh[ord - 1] * hh[ord - 1] + 4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] +
		    2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]);
		denomsum = hh[ord - 1] + hh[ord] + hh[ord + 1];
		dd = denomsum * denom;
		if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x0<=x<=x1 */
			coefa = -4. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) /
				(hh[ord - 1] * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 0.;
			coefc = 12. / denom;
			coefd = 0.;

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x1<=x<=x2 */
			coefa = 4. * (2. * hh[ord - 1] * hh[ord - 1] + 6. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + 3. * hh[ord - 1] * hh[ord + 1] +
				      3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) /
				(hh[ord] * (hh[ord - 1] + hh[ord]) * (hh[ord] + hh[ord + 1]) * dd);
			coefb = -12. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) /
				((hh[ord - 1] + hh[ord]) * dd);
			coefc = 12. * (-2. * hh[ord - 1] * hh[ord - 1] + hh[ord] * hh[ord] + hh[ord] * hh[ord + 1]) /
				((hh[ord - 1] + hh[ord]) * dd);
			coefd = 4. * hh[ord - 1] * (4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + 2. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1]) /
				((hh[ord - 1] + hh[ord]) * dd);

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		} else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x2<=x<=x3 */
			dd *= (hh[ord] + hh[ord + 1]);
			coefa = -4. * (2. * hh[ord - 1] + hh[ord]) / (hh[ord + 1] * dd);
			coefb = 12. * (2. * hh[ord - 1] + hh[ord]) / dd;
			coefc = -12. * (2. * hh[ord - 1] + hh[ord]) * hh[ord + 1] / dd;
			coefd = 4. * (2. * hh[ord - 1] + hh[ord]) * hh[ord + 1] * hh[ord + 1] / dd;

			rho_x = coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefb * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefc * (xi - knot[ord + 1]);
			rho_x += coefd;
		} else		/* x>x3 */
			rho_x = 0.0;
	} else if (ord == Nx - 1) {	/* RHS-1 */
		double		denom   , denomsum, dd;
		denom = hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord];
		denomsum = hh[ord - 2] + hh[ord - 1] + hh[ord];
		dd = denomsum * denom;
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. * (hh[ord - 1] + 2. * hh[ord]) / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * dd);
			coefb = coefc = coefd = 0.0;

			rho_x = coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefb * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefc * (xi - knot[ord - 2]);
			rho_x += coefd;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 3. * hh[ord - 2] * hh[ord] +
				       6. * hh[ord - 1] * hh[ord] + 2. * hh[ord] * hh[ord]) /
				(hh[ord - 1] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 12. * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);
			coefc = 12. * hh[ord - 2] * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);
			coefd = 4. * hh[ord - 2] * hh[ord - 2] * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			dd *= (hh[ord - 1] + hh[ord]);
			coefa = 4. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / (hh[ord] * dd);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / dd;
			coefc = 12. * (-hh[ord - 2] * hh[ord - 1] - hh[ord - 1] * hh[ord - 1] + 2. * hh[ord] * hh[ord]) / dd;
			coefd = 4. * hh[ord] * (3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord]) / dd;

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		} else		/* x>x4 */
			rho_x = 0.0;
	} else if (ord == Nx) {	/* RHS */
		double		denom;
		denom = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1]);
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * denom);
			coefb = coefc = coefd = 0.0;

			rho_x = coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefb * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefc * (xi - knot[ord - 2]);
			rho_x += coefd;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord - 1] * denom);
			coefb = 12 / denom;
			coefc = 12 * hh[ord - 2] / denom;
			coefd = 4. * hh[ord - 2] * hh[ord - 2] / denom;

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		} else		/* x>x2 */
			rho_x = 0.0;
	} else {		/* Away from borders */
		double		denom1  , denom2, denom;
		denom1 = hh[ord - 2] + hh[ord - 1] + hh[ord] + hh[ord + 1];
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * denom1);
			coefb = coefc = coefd = 0.;

			rho_x = coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefb * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefc * (xi - knot[ord - 2]);
			rho_x += coefd;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			denom2 = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]);
			denom = denom1 * denom2;

			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] +
				       4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] + hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]) /
				(hh[ord - 1] * (hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]) * denom);
			coefb = 12. / denom;
			coefc = 12. * hh[ord - 2] / denom;
			coefd = 4. * hh[ord - 2] * hh[ord - 2] / denom;

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			denom2 = (hh[ord - 1] + hh[ord]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = 4. * (hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] +
				      hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) /
				(hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) / denom;
			coefc = 12. * (-hh[ord - 2] * hh[ord - 1] - hh[ord - 1] * hh[ord - 1] + hh[ord] * hh[ord] + hh[ord] * hh[ord + 1]) / denom;
			coefd = 4. * (2. * hh[ord - 2] * hh[ord - 1] * hh[ord] + 2. * hh[ord - 1] * hh[ord - 1] * hh[ord] + hh[ord - 2] * hh[ord] * hh[ord] +
				      2. * hh[ord - 1] * hh[ord] * hh[ord] + hh[ord - 2] * hh[ord - 1] * hh[ord + 1] + hh[ord - 1] * hh[ord - 1] * hh[ord + 1] +
				      hh[ord - 2] * hh[ord] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord] * hh[ord + 1]) / denom;

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		} else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x3<=x<=x4 */
			denom2 = (hh[ord] + hh[ord + 1]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = -4. / (hh[ord + 1] * denom);
			coefb = 12 / denom;
			coefc = -12 * hh[ord + 1] / denom;
			coefd = 4. * hh[ord + 1] * hh[ord + 1] / denom;

			rho_x = coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefb * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefc * (xi - knot[ord + 1]);
			rho_x += coefd;
		} else		/* x>x4 */
			rho_x = 0.0;
	}
	/* free_farray1(hh, 0); */
	return (rho_x);
}

void 
fill_hh(hh, knot, Nx)
	double          *hh, *knot;
	int		Nx;
{
	int		ii;

	for (ii = 0; ii < Nx; ii++)
		hh[ii] = knot[ii + 1] - knot[ii];
}

/* 1st derivative of spl */

double 
d_spl(ord, nknots, knot, xi)
	int		ord       , nknots;
	double          *knot, xi;
/*
 * ord: order of rho(x) nknots : # of knkots = Nx+1 (Nx=index of highest
 * spline) xi: point of interest splines defined as : f_i(x) = a_i(x-x_i)^3 +
 * b_i(x-x_i)^2 + c_i(x-x_i) + d_i
 */
{
	int		ii        , Nx;
	/* double          *hh, rho_x; */
	double           rho_x;
	double           hh[NKNOTS_MAX];
	double		coefa   , coefb, coefc;

	Nx = nknots - 1;
	/* Compute vector hh of spacings */
	/* hh = farray1(0, Nx - 1); */
	fill_hh(hh, knot, Nx);

	/* Consistency checks */
	if ((xi - (double)TOL) > knot[Nx]) {
		/* printf("xi=%g / knot[%d]=%g", xi, Nx, knot[Nx]);
		stop("spl: xi>knot[Nx]"); */
		return 0.0;
	} else if ((xi + (double)TOL) < knot[0]) {
		/* printf("xi=%g / knot[0]=%g", xi, knot[0]);
		stop("spl: xi<knot[0]"); */
		return 0.0;
	} else if (ord > Nx) {
		/* stop("spl: order > Nx"); */
		printf("Warning [%s]: spl index %i exceeds knot count %i\n", __func__, ord, Nx);
		return 0.0;
	}


	if (ord == 0) {		/* LHS */
		double		denom;
		denom = 3. * hh[ord] * hh[ord] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1];
		if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 0.0;
			coefc = -12. / denom;

			rho_x = 3. * coefa * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += 2. * coefb * (xi - knot[ord]);
			rho_x += coefc;
		} else if (xi > knot[ord + 1] & xi <= knot[ord + 2]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord + 1] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 12 / ((hh[ord] + hh[ord + 1]) * denom);
			coefc = -12. * hh[ord + 1] / ((hh[ord] + hh[ord + 1]) * denom);

			rho_x = 3. * coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += 2. * coefb * (xi - knot[ord + 1]);
			rho_x += coefc;
		} else		/* x>x2 */
			rho_x = 0.0;
	} else if (ord == 1) {	/* LHS+1 */
		double		denom   , denomsum, dd;
		denom = (3. * hh[ord - 1] * hh[ord - 1] + 4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] +
		    2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]);
		denomsum = hh[ord - 1] + hh[ord] + hh[ord + 1];
		dd = denomsum * denom;
		if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x0<=x<=x1 */
			coefa = -4. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) /
				(hh[ord - 1] * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 0.;
			coefc = 12. / denom;

			rho_x = 3. * coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += 2. * coefb * (xi - knot[ord - 1]);
			rho_x += coefc;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x1<=x<=x2 */
			coefa = 4. * (2. * hh[ord - 1] * hh[ord - 1] + 6. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + 3. * hh[ord - 1] * hh[ord + 1] +
				      3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) /
				(hh[ord] * (hh[ord - 1] + hh[ord]) * (hh[ord] + hh[ord + 1]) * dd);
			coefb = -12. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) /
				((hh[ord - 1] + hh[ord]) * dd);
			coefc = 12. * (-2. * hh[ord - 1] * hh[ord - 1] + hh[ord] * hh[ord] + hh[ord] * hh[ord + 1]) /
				((hh[ord - 1] + hh[ord]) * dd);

			rho_x = 3. * coefa * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += 2. * coefb * (xi - knot[ord]);
			rho_x += coefc;
		} else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x2<=x<=x3 */
			dd *= (hh[ord] + hh[ord + 1]);
			coefa = -4. * (2. * hh[ord - 1] + hh[ord]) / (hh[ord + 1] * dd);
			coefb = 12. * (2. * hh[ord - 1] + hh[ord]) / dd;
			coefc = -12. * (2. * hh[ord - 1] + hh[ord]) * hh[ord + 1] / dd;

			/*
			 * rho_x=3.*coefa*(xi-knot[ord+1])*(xi-knot[ord+1])*(x
			 * i-knot[ord+1]);
			 */
			rho_x = 3. * coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += 2. * coefb * (xi - knot[ord + 1]);
			rho_x += coefc;
		} else		/* x>x3 */
			rho_x = 0.0;
	} else if (ord == Nx - 1) {	/* RHS-1 */
		double		denom   , denomsum, dd;
		denom = hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord];
		denomsum = hh[ord - 2] + hh[ord - 1] + hh[ord];
		dd = denomsum * denom;
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. * (hh[ord - 1] + 2. * hh[ord]) / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * dd);
			coefb = coefc = 0.0;

			rho_x = 3. * coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += 2. * coefb * (xi - knot[ord - 2]);
			rho_x += coefc;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 3. * hh[ord - 2] * hh[ord] +
				       6. * hh[ord - 1] * hh[ord] + 2. * hh[ord] * hh[ord]) /
				(hh[ord - 1] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 12. * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);
			coefc = 12. * hh[ord - 2] * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);

			rho_x = 3. * coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += 2. * coefb * (xi - knot[ord - 1]);
			rho_x += coefc;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			dd *= (hh[ord - 1] + hh[ord]);
			coefa = 4. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / (hh[ord] * dd);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / dd;
			coefc = 12. * (-hh[ord - 2] * hh[ord - 1] - hh[ord - 1] * hh[ord - 1] + 2. * hh[ord] * hh[ord]) / dd;

			rho_x = 3. * coefa * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += 2. * coefb * (xi - knot[ord]);
			rho_x += coefc;
		} else		/* x>x4 */
			rho_x = 0.0;
	} else if (ord == Nx) {	/* RHS */
		double		denom;
		denom = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1]);
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * denom);
			coefb = coefc = 0.0;

			rho_x = 3. * coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += 2. * coefb * (xi - knot[ord - 2]);
			rho_x += coefc;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord - 1] * denom);
			coefb = 12 / denom;
			coefc = 12 * hh[ord - 2] / denom;

			rho_x = 3. * coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += 2. * coefb * (xi - knot[ord - 1]);
			rho_x += coefc;
		} else		/* x>x2 */
			rho_x = 0.0;
	} else {		/* Away from borders */
		double		denom1  , denom2, denom;
		denom1 = hh[ord - 2] + hh[ord - 1] + hh[ord] + hh[ord + 1];
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * denom1);
			coefb = coefc = 0.;

			rho_x = 3. * coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += 2. * coefb * (xi - knot[ord - 2]);
			rho_x += coefc;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			denom2 = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]);
			denom = denom1 * denom2;

			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] +
				       4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] + hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]) /
				(hh[ord - 1] * (hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]) * denom);
			coefb = 12. / denom;
			coefc = 12. * hh[ord - 2] / denom;

			rho_x = 3. * coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += 2. * coefb * (xi - knot[ord - 1]);
			rho_x += coefc;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			denom2 = (hh[ord - 1] + hh[ord]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = 4. * (hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] +
				      hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) /
				(hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) / denom;
			coefc = 12. * (-hh[ord - 2] * hh[ord - 1] - hh[ord - 1] * hh[ord - 1] + hh[ord] * hh[ord] + hh[ord] * hh[ord + 1]) / denom;

			rho_x = 3. * coefa * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += 2. * coefb * (xi - knot[ord]);
			rho_x += coefc;
		} else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x3<=x<=x4 */
			denom2 = (hh[ord] + hh[ord + 1]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = -4. / (hh[ord + 1] * denom);
			coefb = 12. / denom;
			coefc = -12. * hh[ord + 1] / denom;

			rho_x = 3. * coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += 2. * coefb * (xi - knot[ord + 1]);
			rho_x += coefc;
		} else		/* x>x4 */
			rho_x = 0.0;
	}
	/* free_farray1(hh, 0); */
	return (rho_x);
}

/* 2nd derivative of spl */

double 
d2_spl(ord, nknots, knot, xi)
	int		ord       , nknots;
	double          *knot, xi;
/*
 * ord: order of rho(x) nknots : # of knkots = Nx+1 (Nx=index of highest
 * spline) xi: point of interest splines defined as : f_i(x) = a_i(x-x_i)^3 +
 * b_i(x-x_i)^2 + c_i(x-x_i) + d_i
 */
{
	int		ii        , Nx;
	/* double          *hh, rho_x; */
	double           rho_x;
	double           hh[NKNOTS_MAX];
	double		coefa   , coefb;

	Nx = nknots - 1;
	/* Compute vector hh of spacings */
	/* hh = farray1(0, Nx - 1); */
	fill_hh(hh, knot, Nx);

	/* Consistency checks */
	if ((xi - (double)TOL) > knot[Nx]) {
		/* printf("xi=%g / knot[%d]=%g", xi, Nx, knot[Nx]);
		stop("spl: xi>knot[Nx]"); */
		return 0.0;
	} else if ((xi + (double)TOL) < knot[0]) {
		/* printf("xi=%g / knot[0]=%g", xi, knot[0]);
		stop("spl: xi<knot[0]"); */
		return 0.0;
	} else if (ord > Nx) {
		/* stop("spl: order > Nx"); */
		printf("Warning [%s]: spl index %i exceeds knot count %i\n", __func__, ord, Nx);
		return 0.0;
	}


	if (ord == 0) {		/* LHS */
		double		denom;
		denom = 3. * hh[ord] * hh[ord] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1];
		if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 0.0;

			rho_x = 6. * coefa * (xi - knot[ord]);
			rho_x += 2. * coefb;
		} else if (xi > knot[ord + 1] & xi <= knot[ord + 2]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord + 1] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 12. / ((hh[ord] + hh[ord + 1]) * denom);

			rho_x = 6. * coefa * (xi - knot[ord + 1]);
			rho_x += 2. * coefb;
		} else		/* x>x2 */
			rho_x = 0.0;
	} else if (ord == 1) {	/* LHS+1 */
		double		denom   , denomsum, dd;
		denom = (3. * hh[ord - 1] * hh[ord - 1] + 4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] +
		    2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]);
		denomsum = hh[ord - 1] + hh[ord] + hh[ord + 1];
		dd = denomsum * denom;
		if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x0<=x<=x1 */
			coefa = -4. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) /
				(hh[ord - 1] * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 0.;

			rho_x = 6. * coefa * (xi - knot[ord - 1]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x1<=x<=x2 */
			coefa = 4. * (2. * hh[ord - 1] * hh[ord - 1] + 6. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + 3. * hh[ord - 1] * hh[ord + 1] +
				      3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) /
				(hh[ord] * (hh[ord - 1] + hh[ord]) * (hh[ord] + hh[ord + 1]) * dd);
			coefb = -12. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) /
				((hh[ord - 1] + hh[ord]) * dd);

			rho_x = 6. * coefa * (xi - knot[ord]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x2<=x<=x3 */
			dd *= (hh[ord] + hh[ord + 1]);
			coefa = -4. * (2. * hh[ord - 1] + hh[ord]) / (hh[ord + 1] * dd);
			coefb = 12. * (2. * hh[ord - 1] + hh[ord]) / dd;

			rho_x = 6. * coefa * (xi - knot[ord + 1]);
			rho_x += 2. * coefb;
		} else		/* x>x3 */
			rho_x = 0.0;
	} else if (ord == Nx - 1) {	/* RHS-1 */
		double		denom   , denomsum, dd;
		denom = hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord];
		denomsum = hh[ord - 2] + hh[ord - 1] + hh[ord];
		dd = denomsum * denom;
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. * (hh[ord - 1] + 2. * hh[ord]) / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * dd);
			coefb = 0.;

			rho_x = 6. * coefa * (xi - knot[ord - 2]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 3. * hh[ord - 2] * hh[ord] +
				       6. * hh[ord - 1] * hh[ord] + 2. * hh[ord] * hh[ord]) /
				(hh[ord - 1] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 12. * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);

			rho_x = 6. * coefa * (xi - knot[ord - 1]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			dd *= (hh[ord - 1] + hh[ord]);
			coefa = 4. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / (hh[ord] * dd);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / dd;

			rho_x = 6. * coefa * (xi - knot[ord]);
			rho_x += 2. * coefb;
		} else		/* x>x4 */
			rho_x = 0.0;
	} else if (ord == Nx) {	/* RHS */
		double		denom;
		denom = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1]);
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * denom);
			coefb = 0.;

			rho_x = 6. * coefa * (xi - knot[ord - 2]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord - 1] * denom);
			coefb = 12. / denom;

			rho_x = 6. * coefa * (xi - knot[ord - 1]);
			rho_x += 2. * coefb;
		} else		/* x>x2 */
			rho_x = 0.0;
	} else {		/* Away from borders */
		double		denom1  , denom2, denom;
		denom1 = hh[ord - 2] + hh[ord - 1] + hh[ord] + hh[ord + 1];
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * denom1);
			coefb = 0.;

			rho_x = 6. * coefa * (xi - knot[ord - 2]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			denom2 = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]);
			denom = denom1 * denom2;

			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] +
				       4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] + hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]) /
				(hh[ord - 1] * (hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]) * denom);
			coefb = 12. / denom;

			rho_x = 6. * coefa * (xi - knot[ord - 1]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			denom2 = (hh[ord - 1] + hh[ord]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = 4. * (hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] +
				      hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) /
				(hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) / denom;

			rho_x = 6. * coefa * (xi - knot[ord]);
			rho_x += 2. * coefb;
		} else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x3<=x<=x4 */
			denom2 = (hh[ord] + hh[ord + 1]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = -4. / (hh[ord + 1] * denom);
			coefb = 12. / denom;

			rho_x = 6. * coefa * (xi - knot[ord + 1]);
			rho_x += 2. * coefb;
		} else		/* x>x4 */
			rho_x = 0.0;
	}
	/* free_farray1(hh, 0); */
	return (rho_x);
}
