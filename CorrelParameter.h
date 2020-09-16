// parameter holding data needed for correlation routine

#ifndef CORRELPAR_H
#define CORRELPAR_H

#include "Go4StatusBase/TGo4Parameter.h"
#include <iomanip>
#include <fstream>
class CorrelParameter : public TGo4Parameter {

public:
          CorrelParameter();
          CorrelParameter(const Text_t* name);
  virtual ~CorrelParameter();
  Int_t    PrintParameter(Text_t *buf, Int_t);
  Bool_t   UpdateFrom(TGo4Parameter *);
  int      IsData(std::ifstream &f);

  ///Define correlations of interest
 Bool_t GSetup_corr_FRS_Ge, GSetup_corr_FRS_Ge_long, GSetup_corr_FRS_fat, GSetup_corr_FRS_Aida ;
 Bool_t GSetup_corr_Beta_Gamma, GSetup_corr_Beta_Gamma_Gamma, GSetup_corr_Beta_Gamma_Fatima;
   ///Define Long isomer Ge conditions
 Int_t GLongGate;
 Int_t GFRS_Ge_LongIso_THigh, GFRS_Ge_LongIso_TLow, GFRS_Ge_LongIso_HBin,GFRS_Ge_LongIso_HLow, GFRS_Ge_LongIso_HHigh, GFRS_Ge_LongIso_TScale ;
 Bool_t GFRS_Ge_LongIso_incprmt;

 ///Define Beta-delayed gamma conditions
 Float_t GAidaImpDecT_Low, GAidaImpDecT_High;
 Int_t GAidaImpDecT_HBin;
 
 Int_t GAidaFB_dT, GAidaFB_dE;
 
 Int_t GAIDA_DecEFront, GAIDA_DecEBack;
 
 Int_t GAida_Ge_WRdT_Low, GAida_Ge_WRdT_High;
 
 Int_t GAida_Fat_WRdT_Low, GAida_Fat_WRdT_High;

  ///Define other conditions 

 Int_t GFat_TRawgate_low, GFat_TRawgate_high;
 Int_t GFRS_Fat_TLow, GFRS_Fat_THigh;
 Int_t GFRS_Ge_TLow, GFRS_Ge_THigh;

 Int_t GFRS_AIDA_TLow, GFRS_AIDA_THigh;
 Int_t GAIDA_bPlas_TLow, GAIDA_bPlas_THigh;
 Int_t GbPlast_Egate_low, GbPlast_Egate_high;
// Int_t GbPlast_GAL_TLow, GbPlast_GAL_THigh;
 


// Bool_t GFRS_AIDA_DSSD1,  GFRS_AIDA_DSSD2, GFRS_AIDA_DSSD3;
//  Double_t Aff[60], Bff[60];
//  Double_t Ab[32], Bb[32];
//  Double_t Ag[4], Bg[4];

ClassDef(CorrelParameter,6)
};

#endif
