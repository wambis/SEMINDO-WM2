#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cminos.h"

#define min(a,b) ((a)>(b)?(b):(a))

/* compile-time storage size for n, l, w, and U */
#define MODE_MAX 40000

/* global (unit) scope arrays for n, l, w, c, U, q, and eigenfunctions */
int mode_n[MODE_MAX];
int mode_l[MODE_MAX];
double mode_w[MODE_MAX];
double mode_c[MODE_MAX];
double mode_U[MODE_MAX];
double mode_q[MODE_MAX];
double mode_eig[MODE_MAX][6 * MODEL_NPTS_MAX];
int mode_eig_nvec[MODE_MAX];

/* static variables ... */

/* mode and eigenfunction counters */
static int nmode = 0;
static int neig = 0;

/* sanity checks */
static int minos_is_configured = 0;
static int minos_has_model = 0;
static int minos_stores_eigs = 0;

/* persistent configuration and model storage */
static config_t cminos_config;
static model_t cminos_ref;

/* prototype for the minos calling interface */
extern void run_minos_( );


void
_configure_minos( double eps, double wgrav, int jcom, int lmin, int lmax, 
                  double wmin, double wmax, int nmin, int nmax )
{
#ifdef _DEBUG_CMINOS
  printf( " eps   = %e\n", eps );
  printf( " wgrav = %f\n", wgrav );
  printf( " jcom  = %i\n", jcom );
  printf( " lmin  = %i\n", lmin );
  printf( " lmax  = %i\n", lmax );
  printf( " nmin  = %i\n", nmin );
  printf( " nmax  = %i\n", nmax );
  printf( " wmin  = %f\n", wmin );
  printf( " wmax  = %f\n", wmax );
#endif
  cminos_config.eps   = eps;
  cminos_config.wgrav = wgrav; // mHz
  cminos_config.jcom  = jcom;
  cminos_config.lmin  = lmin;
  cminos_config.lmax  = lmax;
  cminos_config.nmin  = nmin;
  cminos_config.nmax  = nmax;
  cminos_config.wmin  = wmin; // mHz
  cminos_config.wmax  = wmax; // mHz
  minos_is_configured = 1;
#ifdef _DEBUG_CMINOS
  printf( "%s: done\n", __func__ );
#endif
}

int
_run_minos( char *model_name, double tref, int npts, int icb, int cmb, int noc, 
            double *model_vec, int ndata,
            int *n_, int ndim_, 
            int *l_, int ldim_,
            double *w_, int wdim_,
            double *c_, int cdim_,
            double *U_, int Udim_,
            double *q_, int qdim_ )
{
  int min_dim = min( ndim_, min( ldim_, min( wdim_, min( cdim_, min( Udim_, qdim_ ) ) ) ) );
  if ( ndata / 9 > MODEL_NPTS_MAX )
    {
      printf( "Error: [ %s ] MODEL_NPTS_MAX is too small!\n", __func__ );
      return 0;
    }
#ifdef _DEBUG_CMINOS
  printf( "%s: model %s has %i layers\n", __func__, model_name, ndata / 9 );
  printf( "%s: output arrays are of dimension %i\n", __func__, ndim_ );
#endif
  strncpy( cminos_ref.model_name, model_name, MODEL_LINE_LEN );
  cminos_ref.ifanis = 1;
  cminos_ref.tref   = tref;
  cminos_ref.ifdeck = 1;
  cminos_ref.npts   = ndata / 9;
  cminos_ref.icb    = icb;
  cminos_ref.cmb    = cmb;
  cminos_ref.noc    = noc;
  int ipt;
  for ( ipt = 0; ipt < cminos_ref.npts; ipt++ )
    {
      cminos_ref.r[ipt]   = model_vec[9 * ipt    ];
      cminos_ref.rho[ipt] = model_vec[9 * ipt + 1];
      cminos_ref.vpv[ipt] = model_vec[9 * ipt + 2];
      cminos_ref.vsv[ipt] = model_vec[9 * ipt + 3];
      cminos_ref.qk[ipt]  = model_vec[9 * ipt + 4];
      cminos_ref.qmu[ipt] = model_vec[9 * ipt + 5];
      cminos_ref.vph[ipt] = model_vec[9 * ipt + 6];
      cminos_ref.vsh[ipt] = model_vec[9 * ipt + 7];
      cminos_ref.eta[ipt] = model_vec[9 * ipt + 8];
    }
#ifdef _DEBUG_CMINOS
  printf( " %f %f %f %f %f %f %f %f %f\n",
        cminos_ref.r[cminos_ref.npts - 1],
        cminos_ref.rho[cminos_ref.npts - 1],
        cminos_ref.vpv[cminos_ref.npts - 1],
        cminos_ref.vsv[cminos_ref.npts - 1],
        cminos_ref.qk[cminos_ref.npts - 1],
        cminos_ref.qmu[cminos_ref.npts - 1],
        cminos_ref.vph[cminos_ref.npts - 1],
        cminos_ref.vsh[cminos_ref.npts - 1],
        cminos_ref.eta[cminos_ref.npts - 1] );
#endif
  if ( minos_is_configured != 1 )
    {
      printf( "Error: [ %s ] has been called without first calling configure_minos( )\n", __func__ );
      return 0;
    }
  nmode = 0;
  minos_stores_eigs = 0;
  run_minos_( );
  if ( nmode == 0 )
    {
      printf( "[ %s ] Warning: nmode = 0\n", __func__ );
      return 0;
    }
  if ( nmode > min_dim )
    {
      printf( "Error: [ %s ] nmode is larger than smallest dimension of output arrays (nmode = %i)\n", __func__, nmode );
      return 0;
    }
  int imode;
  for ( imode = 0; imode < nmode; imode++ )
    {
      n_[imode] = mode_n[imode];
      l_[imode] = mode_l[imode];
      w_[imode] = mode_w[imode];
      c_[imode] = mode_c[imode];
      U_[imode] = mode_U[imode];
      q_[imode] = mode_q[imode];
    }
  return nmode;
}

int
_run_minos_eig( char *model_name, double tref, int npts, int icb, int cmb, int noc, 
                double *model_vec, int ndata,
                int *n_, int ndim_, 
                int *l_, int ldim_,
                double *w_, int wdim_,
                double *c_, int cdim_,
                double *U_, int Udim_,
                double *q_, int qdim_, 
                double *eig_, int eigdim_ )
{
  int min_dim = min( ndim_, min( ldim_, min( wdim_, min( cdim_, min( Udim_, qdim_ ) ) ) ) );
  if ( ndata / 9 > MODEL_NPTS_MAX )
    {
      printf( "Error: [ %s ] MODEL_NPTS_MAX is too small!\n", __func__ );
      return 0;
    }
#ifdef _DEBUG_CMINOS
  printf( "%s: model %s has %i layers\n", __func__, model_name, ndata / 9 );
  printf( "%s: output arrays are of dimension %i\n", __func__, ndim_ );
#endif
  strncpy( cminos_ref.model_name, model_name, MODEL_LINE_LEN );
  cminos_ref.ifanis = 1;
  cminos_ref.tref   = tref;
  cminos_ref.ifdeck = 1;
  cminos_ref.npts   = ndata / 9;
  cminos_ref.icb    = icb;
  cminos_ref.cmb    = cmb;
  cminos_ref.noc    = noc;
  int ipt;
  for ( ipt = 0; ipt < cminos_ref.npts; ipt++ )
    {
      cminos_ref.r[ipt]   = model_vec[9 * ipt    ];
      cminos_ref.rho[ipt] = model_vec[9 * ipt + 1];
      cminos_ref.vpv[ipt] = model_vec[9 * ipt + 2];
      cminos_ref.vsv[ipt] = model_vec[9 * ipt + 3];
      cminos_ref.qk[ipt]  = model_vec[9 * ipt + 4];
      cminos_ref.qmu[ipt] = model_vec[9 * ipt + 5];
      cminos_ref.vph[ipt] = model_vec[9 * ipt + 6];
      cminos_ref.vsh[ipt] = model_vec[9 * ipt + 7];
      cminos_ref.eta[ipt] = model_vec[9 * ipt + 8];
    }
#ifdef _DEBUG_CMINOS
  printf( " %f %f %f %f %f %f %f %f %f\n",
        cminos_ref.r[cminos_ref.npts - 1],
        cminos_ref.rho[cminos_ref.npts - 1],
        cminos_ref.vpv[cminos_ref.npts - 1],
        cminos_ref.vsv[cminos_ref.npts - 1],
        cminos_ref.qk[cminos_ref.npts - 1],
        cminos_ref.qmu[cminos_ref.npts - 1],
        cminos_ref.vph[cminos_ref.npts - 1],
        cminos_ref.vsh[cminos_ref.npts - 1],
        cminos_ref.eta[cminos_ref.npts - 1] );
#endif
  if ( minos_is_configured != 1 )
    {
      printf( "Error: [ %s ] has been called without first calling configure_minos( )\n", __func__ );
      return 0;
    }
  nmode = neig = 0;
  minos_stores_eigs = 1;
  run_minos_( );
  if ( nmode == 0 )
    {
      printf( "[ %s ] Warning: nmode = 0\n", __func__ );
      return 0;
    }
  if ( minos_stores_eigs == 1 && neig != nmode )
    {
      printf( "Error: [ %s ] nmode != neig (nmode = %i, neig = %i)\n", __func__, nmode, neig );
      return 0;
    }
  if ( nmode > min_dim )
    {
      printf( "Error: [ %s ] nmode is larger than smallest dimension of output arrays (nmode = %i)\n", __func__, nmode );
      return 0;
    }
  if ( ( cminos_config.jcom == 2 && eigdim_ < 2 * cminos_ref.npts ) ||
       ( cminos_config.jcom == 3 && eigdim_ < 6 * cminos_ref.npts ) )
    {
      printf( "Error: [ %s ] eigenfunction output array is too short\n", __func__ );
      return 0;
    }
  int imode, i, ieig = 0;
  for ( imode = 0; imode < nmode; imode++ )
    {
      n_[imode] = mode_n[imode];
      l_[imode] = mode_l[imode];
      w_[imode] = mode_w[imode];
      c_[imode] = mode_c[imode];
      U_[imode] = mode_U[imode];
      q_[imode] = mode_q[imode];
      for ( i = 0; i < mode_eig_nvec[imode]; i++ )
        eig_[ieig++] = mode_eig[imode][i];
      if ( ! ( mode_eig_nvec[imode] == 6 * cminos_ref.npts || mode_eig_nvec[imode] == 2 * cminos_ref.npts ) )
        {
          printf( "Error: [ %s ] mode_eig_nvec does not equal 6 * npts or 2 * npts\n", __func__ );
          return 0;
        }
    }
  return nmode;
}



/* ---- fortran interfaces ---- */

/*
 * Fortran call
      call initialize_model( 
     +     ifanis, tref, ifdeck, n, nic, noc, 
     +     r, rho, vpv, vsv, qkappa, qshear, vph, vsh, eta )
*/
void initialize_model_( int *ifanis, double *tref, int *ifdeck, int *n, int *nic, int *noc,
                        double *r, double *rho, double *vpv, double *vsv, 
                        double *qkappa, double *qshear, double *vph, double *vsh, double *eta )
{
#ifdef _DEBUG_CMINOS
  printf( "%s( )\n", __func__ );
#endif
  *ifanis = cminos_ref.ifanis;
  *tref   = cminos_ref.tref;
  *ifdeck = cminos_ref.ifdeck;
  *n      = cminos_ref.npts;
  *nic    = cminos_ref.icb;
  *noc    = cminos_ref.cmb;
  int ipt;
  for ( ipt = 0; ipt < cminos_ref.npts; ipt++ )
    {
      r[ipt]   = cminos_ref.r[ipt];
      rho[ipt] = cminos_ref.rho[ipt];
      vpv[ipt] = cminos_ref.vpv[ipt];
      vsv[ipt] = cminos_ref.vsv[ipt];
      qkappa[ipt] = cminos_ref.qk[ipt];
      qshear[ipt] = cminos_ref.qmu[ipt];
      vph[ipt] = cminos_ref.vph[ipt];
      vsh[ipt] = cminos_ref.vsh[ipt];
      eta[ipt] = cminos_ref.eta[ipt];
    }
}


/*
 * Fortran call
      call initialize_config( 
     +     eps, wgrav, jcom, lmin, lmax, 
     +     wmin, wmax, normin, normax )
*/
void initialize_config_( double *eps, double *wgrav, int *jcom, int *lmin, int *lmax, 
                         double *wmin, double *wmax, int *nmin, int *nmax )
{
#ifdef _DEBUG_CMINOS
  printf( "%s( )\n", __func__ );
#endif
  *eps   = cminos_config.eps;
  *wgrav = cminos_config.wgrav;
  *jcom  = cminos_config.jcom;
  *lmin  = cminos_config.lmin;
  *lmax  = cminos_config.lmax;
  *wmin  = cminos_config.wmin;
  *wmax  = cminos_config.wmax;
  *nmin  = cminos_config.nmin;
  *nmax  = cminos_config.nmax;
}


/*
 * Fortran call
        call store_mode( nord, l, wmhz, gcom )
*/
void store_mode_( int *n, int *l, double *w, double *c, double *U, double *q, int *store_eigs )
{
#ifdef _DEBUG_CMINOS
  printf( "%s( )\n", __func__ );
#endif
  if ( nmode == MODE_MAX )
    {
      printf( "[ %s ] Error: MODE_MAX is too small\n", __func__ );
      exit( 1 );
    }
  *store_eigs = minos_stores_eigs;
  mode_n[nmode] = *n;
  mode_l[nmode] = *l;
  mode_w[nmode] = *w;
  mode_c[nmode] = *c;
  mode_U[nmode] = *U;
  mode_q[nmode] = *q;
  nmode += 1;
}

/*
 * Fortran call
        call store_mode( nord, l, wmhz, gcom )
*/
void store_eig_( float *eig, int *nvec )
{
  int i;
#ifdef _DEBUG_CMINOS
  printf( "%s( )\n", __func__ );
#endif
  if ( neig == MODE_MAX )
    {
      printf( "[ %s ] Error: MODE_MAX is too small\n", __func__ );
      exit( 1 );
    }
  mode_eig_nvec[neig] = *nvec;
  for ( i = 0; i < *nvec; i++ )
    mode_eig[neig][i] = (double) eig[i];
  neig += 1;
}
