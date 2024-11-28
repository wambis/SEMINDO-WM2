typedef struct { int n,type,l;
                 float w,q,cgp,a_vert,a_hori,phis;
                 float u[222],up[222],v[222],vp[222],ph[222],php[222];
               } eigen222_st;

typedef struct { int n,type,l;
                 float w,q,cgp,a_vert,a_hori,phis;
                 float *u,*up,*v,*vp,*ph,*php;
               } eigen_st;

typedef struct { int n,type,l;
                 float w,q,cgp,a_vert,a_hori,phis;
                 float u[NNDLPNT],up[NNDLPNT],v[NNDLPNT],vp[NNDLPNT],ph[NNDLPNT],php[NNDLPNT];
               } eigen_static_st;

typedef struct { int n,type,l;
                 float w,q,cgp,a_vert,a_hori,phis,ell,rot;
                } eigenhs_st;

typedef struct { int n,type,l;
                 float w,q,cgp,a_vert,a_hori,phis,ell,rot;
                } eigenh_ellrot_st;

typedef struct { int n,type,l;
                 float w,q,cgp,a_vert,a_hori,phis,ell,rot;
                 float *u,*up,*v,*vp,*ph,*php;
               } eigen_ellrot_st;


typedef struct { int n,type,l;
  		 float w,q,cgp,a_vert,a_hori,phis;
		} eigenh_st;
