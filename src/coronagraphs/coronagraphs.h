#if !defined(CORONAGRAPHS_H)
#define CORONAGRAPHS_H



#define CORONAGRAPHS_ARRAYSIZE 4096


int_fast8_t init_coronagraphs();



double coronagraph_make_2Dprolate(double masksizepix, double beamradpix, double centralObs, const char *outname, long size, const char *pupmask_name);

int coronagraph_make_2Dprolateld(double masksizeld, double beamradpix, double centralObs, const char *outname, long size, const char *pupmask_name);

int coronagraph_update_2Dprolate(double masksizeld, double beamradpix, double centralObs, double zfactor);

int coronagraph_make_2Dprolate_CS(double masksize, double centralObs, const char *outname);

int_fast8_t coronagraph_APLCapo_compile();

int coronagraph_init_PIAA();

int coronagraphs_make_SUBARU_pupil();

int coronagraphs_PIAA_apodize_beam(const char *ampl1, const char *opd1, const char *ampl2, const char *opd2);

int coronagraph_simul_AIC(double xld, double yld, const char *psfname);

int coronagraph_simul_4QPM(double xld, double yld, const char *psfname);

int coronagraph_simul_BL8(double xld, double yld, const char *psfname);

int coronagraph_simul_RRPM(double xld, double yld, const char *psfname);

int coronagraph_simul_OVC(double xld, double yld, const char *psfname);

int coronagraph_simul_CPA(double xld, double yld, const char *psfname);

int coronagraph_simul_PIAA(double xld, double yld, const char *psfname);

int coronagraph_simul_PIAAC(double xld, double yld, const char *psfname);

int coronagraph_simul_AIC_PIAAC(double xld, double yld, const char *psfname);

int coronagraph_simul_MULTISTEP_APLC(double xld, double yld, const char *psfname);

int coronagraph_simulPSF(double xld, double yld, const char *psfname, long coronagraph_type, const char *options);

int coronagraph_transm(const char *fname, long coronagraph_type, double logcontrast, const char *options);

int coronagraph_userfunc();

int coronagraph_compute_limitcoeff();

int CORONAGRAPHS_scanPIAACMC_centObs_perf( double obs0input );

#endif
