#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _add_stdp_reg(void);
extern void _cacumm_CA3_reg(void);
extern void _cagk_CA3_reg(void);
extern void _cagkCA3_reg(void);
extern void _cal2_CA3_reg(void);
extern void _can2_CA3_reg(void);
extern void _cat_CA3_reg(void);
extern void _distr_CA3_reg(void);
extern void _E2_Prob_reg(void);
extern void _E3_NMDA_reg(void);
extern void _h_CA3_reg(void);
extern void _ichan2_reg(void);
extern void _KahpM95_CA3_reg(void);
extern void _kaprox_CA3_reg(void);
extern void _kd_CA3_reg(void);
extern void _kdrca1_CA3_reg(void);
extern void _km_CA3_reg(void);
extern void _log_stdp2_reg(void);
extern void _log_stdp_reg(void);
extern void _mult_stdp_reg(void);
extern void _na3n_CA3_reg(void);
extern void _naxn_CA3_reg(void);
extern void _simple_hebbian_reg(void);
extern void _STDPE2Syn_reg(void);
extern void _STDPE2Syn_prob_reg(void);
extern void _vecstim_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," add_stdp.mod");
    fprintf(stderr," cacumm_CA3.mod");
    fprintf(stderr," cagk_CA3.mod");
    fprintf(stderr," cagkCA3.mod");
    fprintf(stderr," cal2_CA3.mod");
    fprintf(stderr," can2_CA3.mod");
    fprintf(stderr," cat_CA3.mod");
    fprintf(stderr," distr_CA3.mod");
    fprintf(stderr," E2_Prob.mod");
    fprintf(stderr," E3_NMDA.mod");
    fprintf(stderr," h_CA3.mod");
    fprintf(stderr," ichan2.mod");
    fprintf(stderr," KahpM95_CA3.mod");
    fprintf(stderr," kaprox_CA3.mod");
    fprintf(stderr," kd_CA3.mod");
    fprintf(stderr," kdrca1_CA3.mod");
    fprintf(stderr," km_CA3.mod");
    fprintf(stderr," log_stdp2.mod");
    fprintf(stderr," log_stdp.mod");
    fprintf(stderr," mult_stdp.mod");
    fprintf(stderr," na3n_CA3.mod");
    fprintf(stderr," naxn_CA3.mod");
    fprintf(stderr," simple_hebbian.mod");
    fprintf(stderr," STDPE2Syn.mod");
    fprintf(stderr," STDPE2Syn_prob.mod");
    fprintf(stderr," vecstim.mod");
    fprintf(stderr, "\n");
  }
  _add_stdp_reg();
  _cacumm_CA3_reg();
  _cagk_CA3_reg();
  _cagkCA3_reg();
  _cal2_CA3_reg();
  _can2_CA3_reg();
  _cat_CA3_reg();
  _distr_CA3_reg();
  _E2_Prob_reg();
  _E3_NMDA_reg();
  _h_CA3_reg();
  _ichan2_reg();
  _KahpM95_CA3_reg();
  _kaprox_CA3_reg();
  _kd_CA3_reg();
  _kdrca1_CA3_reg();
  _km_CA3_reg();
  _log_stdp2_reg();
  _log_stdp_reg();
  _mult_stdp_reg();
  _na3n_CA3_reg();
  _naxn_CA3_reg();
  _simple_hebbian_reg();
  _STDPE2Syn_reg();
  _STDPE2Syn_prob_reg();
  _vecstim_reg();
}
