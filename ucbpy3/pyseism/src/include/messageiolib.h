#define STOP(x) {perror(x); exit(-1);}
#define stop(x) {printf("\n\a error : <%s>\n",x); exit(-1);}
#define usage(x) {printf("\n\a usage : <%s>\n",x);exit(-1);}
#define warning(x) {printf("\n\a <warning> : <%s>\n",x);}
