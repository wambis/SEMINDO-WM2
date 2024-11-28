#include <math.h>

/* angular distance between two points on the sphere (degrees) */
double
getdel(double a0, double o0, double a, double o)
{
	double d2r, r2d;
	double del;
	double q0, sq0, cq0, q, sq, cq, ff, sff, cff, arg;

	d2r = 3.1415926535 / 180.0;
	r2d = 180.0 / 3.1415926535;

	q0 = (90.0 - a0) * d2r;
	sq0 = sin(q0);
	cq0 = cos(q0);

	q = (90.0 - a) * d2r;
	sq = sin(q);
	cq = cos(q);

	ff = (o - o0) * d2r;
	sff = sin(ff);
	cff = cos(ff);

	arg = cq * cq0 + sq * sq0 * cff;
	if (arg > 1.0)
		arg = 1.0;
	if (arg < -1.0)
		arg = -1.0;
	del = r2d * acos(arg);

	return del;
}

/* array version */
void
getdel_vec(double *a0_vec, double *o0_vec, int N, double a, double o, double *res)
{
        int i;
	double d2r, r2d;
	double del;
	double q0, sq0, cq0, q, sq, cq, ff, sff, cff, arg;

	d2r = 3.1415926535 / 180.0;
	r2d = 180.0 / 3.1415926535;

	q = (90.0 - a) * d2r;
	sq = sin(q);
	cq = cos(q);

        for (i = 0; i < N; i++) {
                double a0 = a0_vec[i], o0 = o0_vec[i];

                q0 = (90.0 - a0) * d2r;
                sq0 = sin(q0);
                cq0 = cos(q0);

                ff = (o - o0) * d2r;
                sff = sin(ff);
                cff = cos(ff);

                arg = cq * cq0 + sq * sq0 * cff;
                if (arg > 1.0)
                        arg = 1.0;
                if (arg < -1.0)
                        arg = -1.0;
                res[i] = r2d * acos(arg);
        }
}

/* spherical spline */
double 
spbsp(double hdel, double ahdel)
{
	double x;
	if (hdel < ahdel)
		x = 0.75 * (hdel / ahdel) * (hdel / ahdel) * (hdel / ahdel) - 1.5 * (hdel / ahdel) * (hdel / ahdel) + 1.0;
	else if (hdel <= ahdel * 2.0) {
		x = 0.25 * (2.0 - hdel / ahdel) * (2.0 - hdel / ahdel) * (2.0 - hdel / ahdel);
		if (x < 0)
			x = 0.0;
	} else
		x = 0.0;
	return x;
}

/* spherical spline derivative (w.r.t. delta) */
double 
gspbsp(double hdel, double ahdel)
{
	double x;
	if (hdel < ahdel)
		x = (2.25 * (hdel / ahdel) * (hdel / ahdel) - 3.0 * (hdel / ahdel)) / ahdel;
	else if (hdel <= ahdel * 2.0)
		x = -0.75 * (2.0 - hdel / ahdel) * (2.0 - hdel / ahdel) / ahdel;
	else
		x = 0.0;
	return x;
}
