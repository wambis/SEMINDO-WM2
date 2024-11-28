

typedef struct { short nic,noc,n670,n220,moho,midc,nsl,nsurf;
               }	disc_st;
typedef struct { float rn,wn,vn,gn,rhon;
	       }        norm_st;
typedef struct { short ifanis;
		 short ndisc;
		 disc_st majordisc;
		 short alldisc[MAXDISC];
		 norm_st normfac;
		 float radius[NNDLPNT];
		 float g[NNDLPNT];
		 float ell[NNDLPNT];
		 float eta[NNDLPNT];
		 float qmu[NNDLPNT];
		 float qkappa[NNDLPNT];
		 float rho[NNDLPNT],qrho[NNDLPNT][3];
		 float a[NNDLPNT],qa[NNDLPNT][3];
		 float c[NNDLPNT],qc[NNDLPNT][3];
		 float l[NNDLPNT],ql[NNDLPNT][3];
		 float n[NNDLPNT],qn[NNDLPNT][3];
		 float f[NNDLPNT],qf[NNDLPNT][3];
	       } prem_st;

typedef struct { short ifanis;
		 short ndisc;
		 disc_st majordisc;
		 short alldisc[MAXDISC];
		 norm_st normfac;
		 float *radius;
		 float *g;
		 float *ell;
		 float *eta;
		 float *qmu;
		 float *qkappa;
		 float *rho,**qrho;
		 float *a,**qa;
		 float *c,**qc;
		 float *l,**ql;
		 float *n,**qn;
		 float *f,**qf;
	       } dynprem_st; /* not used */



typedef struct { float w,q,va,ha,grv,ell,a6,a7,a8,a9,a10,rot;
                 float aker[10],u1,u1p,v1,v1p,u2,u2p,v2,v2p;
                 float r1,r2;
                 short n,l;
                 char  typ;   /*S or T*/
               } sprem_st;
