#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double S0 ; 
  double S1 ; 
  double S2 ; 
  double S11 ; 
  double S12 ; 
  double S22 ; 
  double m1 ;
  double m2 ;
  double v11 ;
  double v12 ;
  double v22 ;
  double corr ;
  double Z ;
} CORR; 

int calccorr(CORR *corrpt, int mode, int ztrans)  ;
void printcorr(CORR *corrpt)  ;
void clearcorr(CORR *corrpt)  ;
void addcorr(CORR *corrpt, double x1, double x2) ; 
void addcorr2 (CORR * corrpt, double x0, double x1, double x2, double x12, double x11,  double x22) ;
void addcorrn(CORR *corrpt, double x1, double x2, double yn) ; 
void pluscorr(CORR *out, CORR *c1, CORR *c2) ;
void minuscorr(CORR *out, CORR *c1, CORR *c2) ;

#ifdef __cplusplus
}
#endif
