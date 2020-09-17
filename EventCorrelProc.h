// $Id: EventCorrelProc.h 755 2011-05-20 08:04:11Z linev $
//Mistry 10.04.19
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef EVENTCORRELPROCESSOR_H
#define EVENTCORRELPROCESSOR_H

#include "TGo4EventProcessor.h"

#include "CorrelParameter.h"
#include "CalibParameter.h"
#include "EventUnpackStore.h"
#include "EventAnlStore.h"
#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "Go4ConditionsBase/TGo4WinCond.h"
#include "Go4ConditionsBase/TGo4PolyCond.h"
#include "DESPEC_Array_Sizes.h"

#include "TTree.h"
#include "TFile.h"

//Max lengths
#define  maxdef(a,b) ( ((a) > (b)) ? (a) : (b) )
//Absolute value
#define ABS(x)  ((x)>=0 ? (x):-(x))

class EventAnlStore;
class EventCorrelStore;
class TSCNParameter;

class AidaHitPID : public TObject
{
public:
  int Event;

  int DSSD;
  // Coordinate in Strip & Position
  double StripX;
  double StripY;
  double PosX;
  double PosY;
  

  double Energy;
  double EnergyFront;
  double EnergyBack;

  // Cluster Size for correlations
  int StripXMin;
  int StripXMax;
  int StripYMin;
  int StripYMax;
  int ClusterSizeX;
  int ClusterSizeY;

  int64_t Time;
  int64_t TimeFront;
  int64_t TimeBack;
  int64_t FastTime;
  int64_t FastTimeFront;
  int64_t FastTimeBack;
  
  // Add extra PID stuff for implants
  
  Int_t GatePass;
  Float_t FRS_Z;
  Float_t FRS_AoQ;

  // If the implant is believed to have stopped here (AIDA only veto)
  bool Stopped;

  AidaHitPID() { Zero(); }
  void Zero() {}
  virtual ~AidaHitPID() {}

  ClassDef(AidaHitPID, 3)
};


class EventCorrelProc : public TGo4EventProcessor {
   public:
      EventCorrelProc();
      EventCorrelProc(const char * name);
      virtual ~EventCorrelProc();
       CorrelParameter *fCorrel;
       CalibParameter *fCal;
 
      virtual Bool_t BuildEvent(TGo4EventElement* dest);
        std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
        AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);
        std::multimap<int64_t, AidaHitPID> implantMap;
        
    void get_used_systems();
    void Ge_2DPromptFlashCut();
    void Fat_2DPromptFlashCut();
    void Fat_TimeCorrection(EventAnlStore* cInput);
    void FRS_Gates_corrProc();
    int  IsData(std::ifstream &f);
            
    void Make_FRS_AIDA_Histos();
    //void Make_FRS_Prompt_AIDA_FATIMA_Ge_Histos();
    void Make_FRS_Prompt_Ge_Histos();
    void Make_FRS_LongIso_Ge_Histos();
    void Make_FRS_Prompt_Fat_Histos();
    void Make_bPlast_Fat_Histos();
    void Make_bPlast_Ge_Histos();
    
    //void Make_FRS_AIDA_bPlast_Histos();
    void Make_Beta_Gamma_Histos();
    
    void Make_FRS_Delayed_AIDA_Gamma_Histos();
    void Make_FRS_AIDA_bPlast_FATIMA_Histos();
    void Make_Timemachine_Histos();
     
    void Process_FRS_AIDA(EventAnlStore* cInput, EventCorrelStore* cOutput);
    //void Process_FRS_Prompt_AIDA_FATIMA_Ge(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    void Process_FRS_Prompt_Fat(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_FRS_Prompt_Ge(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    void Process_FRS_LongIso_Ge(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    void Process_bPlast_Fat(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_bPlast_Ge(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    
    void Process_Beta_Gamma(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    void Process_FRS_Delayed_AIDA_Gamma(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_FRS_AIDA_FATIMA(EventAnlStore* cInput, EventCorrelStore* cOutput);
    void Process_Timemachine_Histos(EventAnlStore* cInput, EventCorrelStore* cOutput);
    
    TGo4WinCond* MakeWindowCond(const char* foldername,
                  const char* condname,
                  float left = 0.,
                  float right = 4096.,
                  const char* HistoName = 0);
                                 
        TGo4PolyCond* MakePolyCond(const char* foldername,
                 const char* condname,
                 Int_t size,
                 Float_t (*points)[2],
                 const char* HistoName = 0);
     
    //  std::vector<TH2*> hA_FRS_Z1Z2_implants_strip_xy;
      std::vector<TH2*> hA_FRS_Z1Z2_implants_pos_xy;
      std::vector<TH1*> hA_FRS_Z1Z2_implants_e;
      std::vector<TH2*> hA_FRS_Z1Z2_implants_e_xy;
      
      //For the long isomer correlations
        std::vector<float> tag_all; 
        std::vector<Long64_t> ts_all;
        std::vector<double> GeE_all;
        std::vector<Long64_t> GeT_all;
        Long64_t dT_frsge_long;
        Long64_t dT_frsge_prompt;
        Long64_t dT_frsge_mult_long;
        Long64_t dT_frsge_mult_prompt;
        Long64_t dT_GeT_long;
        Long64_t dT_GeT_prompt;
        Long64_t Ge_FirstT_long;
        Long64_t Ge_FirstT_prompt;
        Int_t Ge_mult_long;
        Int_t Ge_mult_prompt;
        
        double GeE_Long[GALILEO_MAX_HITS];
        double GeT_Long[GALILEO_MAX_HITS];
        double GeE_Prm_Long[GALILEO_MAX_HITS];
        double GeT_Prm_Long[GALILEO_MAX_HITS];
        Long64_t ts;
      
//       std::vector<std::vector<TH2*>> hA_FRS_ZAoQ_implants_strip_xy;
//       std::vector<std::vector<TH1*>> hA_FRS_ZAoQ_implants_e_gate;
      TH2 *hA_FRS_Z1Z2_implants_strip_xy[8][3];
      TH2 *hA_FRS_ZAoQ_implants_strip_xy[8][3];
      TH2 *hA_FRS_Z1Z2_implants_strip_xy_stopped[8][3];
      TH2 *hA_FRS_ZAoQ_implants_strip_xy_stopped[8][3];
      TH1 *hA_FRS_ZAoQ_implants_e[8][3];
      TH1 *hA_FRS_ZAoQ_implants_e_stopped[8][3];
      TH2 *hA_FRS_Z1Z2_x2x4AoQ_implants_strip_xy[8][3];
      TH1 *hA_FRS_Z1Z2_x2x4AoQ_implants_e[8][3];
      TH2 *hA_FRS_Z1Z2_x2x4AoQ_implants_strip_xy_stopped[8][3];
      TH1 *hA_FRS_Z1Z2_x2x4AoQ_implants_e_stopped[8][3];
       
      std::vector<TH2*> hA_FRS_ZAoQ_implants_pos_xy;
     // std::vector<TH1*> hA_FRS_ZAoQ_implants_e;
      
     // std::vector<TH2*> hA_FRS_ZAoQ_implants_strip_xy_gate1;
     
     // std::vector<TH2*> hA_FRS_Z1Z2_x2x4AoQ_implants_e_xy;
      std::vector<TH1*> hA_FRS_ZAoQ_implants_time_delta;
      std::vector<TH1*> hA_FRS_ZAoQ_implants_strip_1d;
      std::vector<TH1*> hA_FRS_ZAoQ_implants_per_event;
      std::vector<TH2*> hA_FRS_ZAoQ_implants_strip_xy_dssdg;
      TH1 *hA_FRS_dT;
      TH1 *hA_bPlast_dT;
      TH1 *hADec_bPlast_dT;
      
      TH1 *hFRS_bPlast_dT;
      
      TH1 *hA_impdec_dT;
      TH1 *hA_impdec_dT_FRS_gated;
      
      
      TGo4WinCond  *fconHis1;
      TGo4WinCond  *fWinCon1;
      
      TH1 *hA_implant_FatdT;
      TH1 *hA_implant_FatE;
      TH2 *hA_implant_FatdT_FatE;
      
      TH1 *hA_FRSWR_FatWR;
      TH1 *hA_FRS_FatE;
      TH1 *hA_FRS_PID_FatE[8];
      TH2 *hA_FRS_FatEvsT[8];
      TH2 *hA_FRS_FatE1vsE2[8];
      
      ///Prompt fatima LT analysis 
      TH1 *hA_FRS_Fat_LT1_start_stop;
      TH1 *hA_FRS_Fat_LT1_start_stop_ns;
      TH1 *hA_FRS_Fat_LT2_start_stop;
      TH1 *hA_FRS_Fat_LT2_start_stop_ns;
      TH1 *hA_FRS_Fat_LT3_start_stop;
      TH1 *hA_FRS_Fat_LT3_start_stop_ns;
      
      TH1 *hA_FRS_Fat_LT1_stop_start;
      TH1 *hA_FRS_Fat_LT1_stop_start_ns;
      TH1 *hA_FRS_Fat_LT2_stop_start;
      TH1 *hA_FRS_Fat_LT2_stop_start_ns;
      TH1 *hA_FRS_Fat_LT3_stop_start;
      TH1 *hA_FRS_Fat_LT3_stop_start_ns;
      
      
      TH1 *hA_FRS_Z1Z2_X2AoQX4AoQ_FatE[8];
      
      TH1 *hA_FRSWR_GeWR;
      TH1 *hA_FRS_GeE;
      TH1 *hA_FRS_PID_GeE[8];
      TH1 *hA_FRS_Z1Z2_X2AoQX4AoQ_GeE[8];
      TH2 *hA_FRS_Z1Z2_X2AoQX4AoQ_GeEvsT[8];
      TH2 *hA_FRS_GeEvsT;
      TH2 *hA_FRS_GeEvsT_Gated[8];
      TH2 *hA_FRS_PID_GeE1_GeE2[8];
      
      TH1 *hA_FRS_PID_GeE_LongIso;
      TH2 *hA_FRS_GeEvsT_LongIsoGated;
      TH2 *hA_FRS_GeE1vsGeE2_LongIsoGated;
      
      TH1 *hGe_BetaGamma;     
      TH1 *hAidaImpDecdT;
      TH1 *hA_Dec_Ge_dT;
      TH1 *hA_Dec_bPlas_dT;
      TH1 *hA_Dec_Fatima_dT;
      TH1 *hGe_BetaGamma_E[8];
      TH1 *hGe_BetaGamma_GeE1_GeE2[8];
      TH1 *hGe_BetaGamma_dT[8];
      TH2 *hGe_BetaGamma_EdT[8];
      TH1 *hGe_BetaGamma_GeE1_GatedTrans;
      
      TH1 *hFat_LT1_start_stop;
      TH1 *hFat_LT1_start_stop_ns;
      TH1 *hFat_LT2_start_stop;
      TH1 *hFat_LT2_start_stop_ns;
      
      TH1 *hFat_LT1_stop_start;
      TH1 *hFat_LT1_stop_start_ns;
      TH1 *hFat_LT2_stop_start;
      TH1 *hFat_LT2_stop_start_ns;
      
      
      TH1 *hFat_95Rh_All_LT;
      TH1 *hFat_BetaGamma_E[8];
      TH2 *hFat_BetaGamma_E1_E2[8];
      
      
      TH1 *hA_bPlastWR_GeWR;
      TH1 *hA_bPlastDet1_GeE;
      TH1 *hA_bPlastDet2_GeE;
      
      TH1 *hA_bPlastWR_FatWR;
      TH1 *hA_bPlastDet1_FatE;
      TH1 *hA_bPlastDet2_FatE;

      TH1 *hA_implant_bPlast_GeE;
      TH1 *hA_implant_bPlast_dT;
      
      std::vector<TH2*> hA_bPlast_decays_strip_xy;
      std::vector<TH1*> hA_bPlast_decays_e;
      TH1 *hA_bPlast_Energy_AidaImpGated[2][16];
      TH1 *hA_bPlast_Energy_AidaDecayGated[2][16];
      TH1 *hA_bPlast_E;
      TH1 *hADecay_bPlast_E;
      TH1 *hA_FRSgated_bPlast_E;
      TH1 *hA_gatedE_bPlast_E;
      
      TH1 *hA_Ge_WRdT;
      TH1 *hA_dT_GeE;
      TH2 *hA_dT_imp_decay_vs_GeE;
      TH1 *hA_dT_FRS_Gated_GeE;
      TH2 *hA_dT_FRS_Gated_imp_decay_vs_GeE;
      
      TH1 *hA_Fat_WRdT;
      TH1 *hA_dT_FatE;
      TH2 *hA_dT_imp_decay_vs_FatE;
      TH1 *hA_dT_FRS_Gated_FatE;
      TH2 *hA_dT_FRS_Gated_imp_decay_vs_FatE;
     
      TH1 *hbPlast_FRS[3][16];
      TH1 *dT_bPlast_FRS;
      TH1 *hA_bPlast_E_Ch[3][16];
      TH1 *hA_gatedE_bPlast_E_Ch[3][16];
      
    //  bool AIDA_implantation_gate;
      TGo4PolyCond  *cAIDA_IMPgate_DSSD1;
      TGo4PolyCond  *cAIDA_IMPgate_DSSD2;
      TGo4PolyCond  *cAIDA_IMPgate_DSSD3;
      TGo4PolyCond  *cGe_EdT_cut[8];
      TGo4PolyCond  *cFat_EdT_cut[8];

      int event_number;
      int Used_Systems[7];
      Long64_t AIDA_WR;
      Long64_t FRS_WR;
      Long64_t bPLAS_WR;
      Long64_t FAT_WR;
      Long64_t GAL_WR;
      Long64_t dT_AIDA_FRS;
      Long64_t dT_AIDA_bPlast;
      Long64_t dT_FRS_bPlast;
      Long64_t dT_AIDA_Ge;
      Long64_t dT_AIDA_Fat;
      TSCNParameter *fParam1;
      
      Int_t Ge_dT_cut_num;
      Float_t X_Ge_EdT_cut[8][8],Y_Ge_EdT_cut[8][8];
      
      Int_t Fat_dT_cut_num;
      Float_t X_Fat_EdT_cut[8][8],Y_Fat_EdT_cut[8][8];
      
      Int_t C_ZAoQgnum,C_Z1Z2gnum,C_X2AoQgnum,C_X4AoQgnum,C_dEdeggnum;
      Float_t C_X_ZAoQ[8][8],C_Y_ZAoQ[8][8];
      
      Int_t Fat_det_num;
      Float_t fat_timecorrection[FAT_MAX_VME_DET];
      
       TGo4PolyCond  *cID_Z_AoQ_corrstep_fat[8]; 
       TGo4PolyCond  *cID_Z_AoQ_corrstep_Ge[8]; 
       TGo4PolyCond  *cID_Z_AoQ_corrstep_BDG[8];     
      
      //Long64_t*** aidatime;

       Long64_t aida_imptime[3][128][128];
       Long64_t aida_imptime_FRS_gated[3][128][128];
      
      ClassDef(EventCorrelProc, 1)
	};
#endif 
