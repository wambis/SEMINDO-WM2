
#ifndef _BSPL_H

#define _BSPL_H

double 
spl(int ord, int nknots, double *knot, double xi);

double 
d_spl(int ord, int nknots, double *knot, double xi);

double 
d2_spl(int ord, int nknots, double *knot, double xi);

#endif
