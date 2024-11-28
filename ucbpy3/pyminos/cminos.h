
#if ! defined(_CMINOS_H)

#define _CMINOS_H

/* maximum row length of tabular minos-format model entry */
#define MODEL_LINE_LEN 256

/* maximum number of layer entries in a minos-format model */
#define MODEL_NPTS_MAX 2000

/* minos configuration information 
 * note: wgrav, wmin, and wmax are in units of mHz */
struct minos_config
{
  int jcom, lmin, lmax, nmin, nmax;
  double wgrav, eps, wmin, wmax;
};
typedef struct minos_config config_t;

struct minos_model
{
  int ifanis, ifdeck, npts, icb, cmb, noc;
  double tref;
  double r[MODEL_NPTS_MAX], rho[MODEL_NPTS_MAX], vpv[MODEL_NPTS_MAX], vsv[MODEL_NPTS_MAX], 
    qk[MODEL_NPTS_MAX], qmu[MODEL_NPTS_MAX], vph[MODEL_NPTS_MAX], vsh[MODEL_NPTS_MAX], eta[MODEL_NPTS_MAX];
  char model_name[MODEL_LINE_LEN];
};
typedef struct minos_model model_t;

#endif
