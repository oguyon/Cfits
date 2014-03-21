
#if !defined(CORONAGRAPHS_H)
#define CORONAGRAPHS_H



#define CORONAGRAPHS_ARRAYSIZE 4096


double coronagraph_make_2Dprolate(double masksize, double centralObs, char *outname);

int coronagraph_make_2Dprolateld(double masksizeld, double centralObs, char *outname);

int coronagraph_update_2Dprolate(double masksizeld, double centralObs, double zfactor);

int coronagraph_make_2Dprolate_CS(double masksize, char *outname);

int coronagraph_APLCapo_compile();

int coronagraph_init_PIAA();

int coronagraphs_make_SUBARU_pupil();

int coronagraphs_PIAA_apodize_beam(char *ampl1, char *opd1, char *ampl2, char *opd2);

int coronagraph_simul_AIC(double xld, double yld, char *psfname);

int coronagraph_simul_4QPM(double xld, double yld, char *psfname);

int coronagraph_simul_BL8(double xld, double yld, char *psfname);

int coronagraph_simul_RRPM(double xld, double yld, char *psfname);

int coronagraph_simul_OVC(double xld, double yld, char *psfname);

int coronagraph_simul_CPA(double xld, double yld, char *psfname);

int coronagraph_simul_PIAA(double xld, double yld, char *psfname);

int coronagraph_simul_PIAAC(double xld, double yld, char *psfname);

int coronagraph_simul_AIC_PIAAC(double xld, double yld, char *psfname);

int coronagraph_simul_MULTISTEP_APLC(double xld, double yld, char *psfname);

int coronagraph_simulPSF(double xld, double yld, char *psfname, long coronagraph_type, char *options);

int coronagraph_transm(char *fname, long coronagraph_type, double logcontrast, char *options);

int coronagraph_userfunc();

int coronagraph_compute_limitcoeff();

int CORONAGRAPHS_scanPIAACMC_centObs_perf( double obs0input );

#endif
