
#ifndef _SPL_DRIVER_H

#define _SPL_DRIVER_H

int
_bspl_driver(
        int deriv, 
        double *knots, int nknots, 
        double *r,     int nr, 
        int k, 
        double *res,   int nres);

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
        int *resloc,        int nresloc);

#endif
