// $Id: EventAnlProc.h 755 2011-05-20 08:04:11Z linev $
// Adapted for DESPEC by A.K.Mistry 2020
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef EVENTANLPROCESSOR_H
#define EVENTANLPROCESSOR_H

#include "TGo4EventProcessor.h"

#include "EventUnpackStore.h"
#include "EventAnlStore.h"
#include "CalibParameter.h"
#include "CorrelParameter.h"
#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "TCutG.h"
#include "TGraph.h"
#include "Go4ConditionsBase/TGo4WinCond.h"
#include "Go4ConditionsBase/TGo4PolyCond.h"
#include "DESPEC_Array_Sizes.h"

#define CYCLE_TIME  (Double_t) 5000
#define COARSE_CT_RANGE  0x800  // 11 bits

///Temp definitions for combined crate
#define bPLAS1_TAMEX_NUM 3
#define bPLAS2_TAMEX_NUM 3
#define bPLAS_TAMEX_ID 2
#define FAT_TAMEX_ID 0
  

class EventAnlStore;


class EventAnlProc : public TGo4EventProcessor {
   public:
      EventAnlProc();
      EventAnlProc(const char * name);
      virtual ~EventAnlProc();
      CalibParameter *fCal;
      CorrelParameter *fCorrel;
    

      virtual Bool_t BuildEvent(TGo4EventElement* dest);
      virtual void UserPostLoop();
      TGo4PolyCond *fCond_FingToTvsStrip;
        void get_used_systems();
        
        double lead_lead_bplas[48][100], trail_trail_bplas[48][100];
        double lead_lead_fat[48][100];
        
        double ToT_bplas[3][16][10] ;
        double lead_bplas[3][16][10] ;
        double trail_bplas[3][16][10] ;
        double lead_lead_bplas_gated[48][100]; 
        double lead_lead_bplas_Ref1[16][10];
        double lead_lead_bplas_Ref2[16][10];
        
        
        double lead_fat[48][10], trail_fat[48][10];
        double lead_lead_fat_Ref1[48][10];
        double ToT_fat[48][10];
        double FatTam_RefCh0[10];
        double SC41L_ANA_lead_fat[10], SC41R_ANA_lead_fat[10];
        double SC41L_DIG_lead_fat[10], SC41R_DIG_lead_fat[10];
        double bPlasDet1_coin_lead_Fat[10], bPlasDet2_coin_lead_Fat[10] ;
        int hits_fat_lead;
        
        double SC41L_ANA_lead_bPlas[48][100], SC41R_ANA_lead_bPlas[48][100];
        double SC41L_DIG_lead_bPlas[48][100], SC41R_DIG_lead_bPlas[48][100];
       
        
        
        double lead_lead_fat_Ref0[48][100];
        int hits_bplas_lead = 0, hits_bplas_trail=0;
        

        double bPlas_TAM_SC41L_ANA;
        double bPlas_TAM_SC41R_ANA;
        double bPlas_TAM_SC41L_DIG;
        double bPlas_TAM_SC41R_DIG;
       
        double bPlas_RefCh0_Det1[10];
        double bPlas_RefCh0_Det2[10];
        
        double Fat_TAM_SC41L_ANA;
        double Fat_TAM_SC41R_ANA;
        double Fat_TAM_SC41L_DIG;
        double Fat_TAM_SC41R_DIG;
        double Fat_RefCh[10];
        double bPlas_AND_Coinc[100];
      
      //FRS Histograms and conditions setup
        TH1I* MakeH1I(const char* foldername,
        const char* histoname,
        Int_t nbins,
        Float_t xmin, Float_t xmax,
        const char* xtitle = "channels",
        Color_t linecolor = 2,
        Color_t fillcolor = 6,
        const char* ytitle = 0);

      TH2I* MakeH2I(const char* foldername,
        const char* histoname,
        Int_t nbinsx, Float_t xmin, Float_t xmax,
        Int_t nbinsy, Float_t ymin, Float_t ymax,
        const char* xtitle, const char* ytitle,
        Color_t marker);

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

         Bool_t Check_PolyCond_Multi_X_Y(Float_t X, Float_t Y, Float_t*** V, int n, int cond_num);
         Bool_t Check_PolyCond_X_Y(Float_t X, Float_t Y, Float_t** V, int n );
         void FRS_Gates();

      void ProcessAida(EventUnpackStore* pInput, EventAnlStore* pOutput);


      int PrcID_Conv[7];
      int Used_Systems[7];
      int event_number;
      bool VMEorTAMEX_bPlas;
      bool VMEorTAMEX_fatima;
      bool VMEandTAMEX_fatima;

      void checkTAMEXorVME();

      Float_t  FRS_Music_dE[2],  FRS_Music_dE_corr[2];
      Float_t FRS_Music_E1[8], FRS_Music_E2[8],FRS_Music_T1[8],FRS_Music_T2[8];
      Float_t  FRS_sci_l[12],  FRS_sci_r[12],  FRS_sci_e[12],  FRS_sci_tx[12],  FRS_sci_x[12];
      Float_t FRS_tpc_x[7],FRS_tpc_y[7];
      Float_t  FRS_sci_tofll2,  FRS_sci_tofll3, FRS_sci_tof2, FRS_sci_tofrr2, FRS_sci_tofrr3, FRS_sci_tof3;
      Float_t  FRS_ID_x2, FRS_ID_y2, FRS_ID_a2, FRS_ID_b2;
      Float_t  FRS_ID_x4, FRS_ID_y4, FRS_ID_a4, FRS_ID_b4;
      Int_t    FRS_sci_dt_21l_21r, FRS_sci_dt_41l_41r, FRS_sci_dt_42l_42r, FRS_sci_dt_43l_43r;
      Int_t    FRS_sci_dt_21l_41l, FRS_sci_dt_21r_41r, FRS_sci_dt_21l_42l, FRS_sci_dt_21r_42r;
      Float_t  FRS_ID_brho[2], FRS_ID_rho[2];
      Float_t  FRS_beta, FRS_beta3, FRS_gamma;
      Float_t  FRS_AoQ, FRS_AoQ_corr;
      Float_t  FRS_z, FRS_z2, FRS_z3;
      Float_t  FRS_znocorr, FRS_z2nocorr;
      Float_t  FRS_dEdeg, FRS_dEdegoQ;
      Float_t  FRS_timestamp, FRS_ts, FRS_ts2;
      bool FRS_spill;
      
      Float_t  X_ZAoQ_gate[8][24];   
      Int_t ZAoQgnum,Z1Z2gnum,X2AoQgnum,X4AoQgnum,dEdeggnum;
      Float_t X_ZAoQ[8][8],Y_ZAoQ[8][8];
      Float_t X_ZZ2[8][6],Y_ZZ2[8][6];
      Float_t XX4_AoQ[8][6], YX4_AoQ[8][6];
      Float_t XX2_AoQ[8][6], YX2_AoQ[8][6];
      Float_t X_dEdeg[8][6], Y_dEdeg[8][6];
       // Float_t** cID_x4AoQ_Z;    
//         Float_t*** cID_x2AoQ;
//         Float_t*** cID_x4AoQ;
//         Float_t** cID_Z_AoQ_test;
//         Float_t*** cID_Z_Z2;
        
//         Bool_t       id_g_x2AoQ[5];
//         Bool_t       id_g_x4AoQ[5];
//         Bool_t       id_g_x4AoQ_Z[5]; 
//         Bool_t      id_g_z_AoQ_test;
//         Bool_t     id_g_z_z2[5];
      int       Aida_Fired;
     // Long64_t WR_Aida_Det_diff[10000];
  
      int   Fatmult;
      int    SC40mult;
      int    SC41mult;
      int   Fat_bPlasmult;
      int    stdcmult;
      int    sqdcmult;
      int    Fat_QDC_ID[50];
      double Fat_QDC_E[50];
      double Fat_QDC_E_Raw[50];

      Long64_t Fat_QDC_T_coarse[50];
      double Fat_QDC_T_fine[50];
      
      int    Fat_TDC_ID[50];
      Long64_t Fat_TDC_T[50];
      Long64_t Fat_TDC_T_Raw[50];
      //Long64_t FatTDC_TS[50][50];
      int    FatTDC_Multipl[50];
      double Fat_Cha_Ref_TDC;
      Long64_t Fat_WR;
      double SC40[50];
      double SC41[50];
      double Fat_bPlas[100];

      int Fat_TDC_Singles_ID[50];
      int Fat_QDC_Singles_ID[50];
      double Fat_QDC_Singles_E[50];
      double Fat_QDC_Singles_E_Raw[50];
      Long64_t Fat_TDC_Singles_t[50];
      Long64_t Fat_TDC_Singles_t_Raw[50];
      Long64_t Fat_QDC_Singles_t_coarse[50];
      double_t Fat_QDC_Singles_t_fine[50];

      
      
       // From unpacker
          int    GalFired;
          int    GalDet[GALILEO_MAX_HITS];
          int    GalCrys[GALILEO_MAX_HITS];
          double GalE[GALILEO_MAX_HITS];
          double GalE_Cal[GALILEO_MAX_HITS];
          double GalEventT[GALILEO_MAX_HITS];
          double GalT[GALILEO_MAX_HITS];
          bool GalPileUp[GALILEO_MAX_HITS];
          bool GalOverFlow[GALILEO_MAX_HITS];

     
      Long64_t Gal_WR;
      Long64_t Gal_WR_Store[100];

      int    Fing_firedTamex;
      int    Fing_leadHits[4];
      int    Fing_trailHits[4];
      int    Fing_iterator[4];
      double Fing_trig[4];
      int    Fing_tamex_ch[4][32];
      int    Fing_leadChan[4][32];
      int    Fing_trailChan[4][32];
      double Fing_lead_coarse[4][32];
      double Fing_lead_fine[4][32];
      double Fing_trail_coarse[4][32];
      double Fing_trail_fine[4][32];

      double Fing_leadT[4][32];
      double Fing_trailT[4][32];
      double Fing_TOT[4][32];
      double Fing_TOT_added[4][32];
      int    Fing_chID[4][32];
      double dataSetPerEvent[50];
      double pmtSetPerEvent[50];
      double totaltimeEvent[50];
      double maxToT;
      int    maxToTChan;
      double maxToT_added;
      int    maxToT_added_Chan;

     void Make_FRS_Histos();
     void Make_Aida_Histos();
     void Make_Plastic_VME_Histos();
     void Make_Plastic_Tamex_Histos();
     void Make_Fatima_Tamex_Histos();
     void Make_Fatima_Histos();
     void Make_Fatima_VME_Tamex_Histos();
     void Make_Galileo_Histos();
     void Make_Finger_Histos();
     void Make_WR_Histos();

     void Make_Fat_Plas_Histos();
     void Make_Fing_Plas_Histos();

     void Do_FRS_Histos(EventAnlStore* pOutput);
     void Do_Plastic_VME_Histos(EventAnlStore* pOutput);
     void Do_Plastic_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
     void Do_Fatima_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
     void Do_Fatima_Histos(EventAnlStore* pOutput);
     
     void Do_Fatima_VME_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
     void Do_Galileo_Histos(EventAnlStore* pOutput);
     void Do_Finger_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput);
     void Do_Fing_Plas_Histos(EventAnlStore* pOutput);
     void Do_Fat_Plas_Histos(EventAnlStore* pOutput);
     void Do_WR_Histos(EventUnpackStore* pInput);
      // TH1 *GermaniumCal;

     void read_setup_parameters();
     void FAT_det_pos_setup();
     double distance_between_detectors(double, double, double, double, double, double);
     double angle_between_detectors(double, double, double);

    long lastTime;
    int ID;
    AidaEvent old;
    AidaEvent evt;
    long startTime;
    long stopTime;

         //   GalEvent gal;

      /* Multiplexer correction */


            int totalEvents;
            int implantEvents;
            int goodImplantEvents;
            int stoppedEvents;
            int decayEvents;
            int pulserEvents;
            int nonsenseEvents;

            std::vector<TH2*> implants_strip_xy;
            std::vector<TH2*> implants_pos_xy;
            std::vector<TH1*> implants_e;
            std::vector<TH2*> implants_e_xy;
            std::vector<TH1*> implants_time_delta;
            std::vector<TH1*> implants_strip_1d;
            std::vector<TH1*> implants_per_event;
            std::vector<TH1*> implants_channels;
            std::vector<TH2*> implants_x_ex;
            std::vector<TH2*> implants_y_ey;
            #ifdef AIDA_PULSER_ALIGN
            TH2* aida_pulser_time;
            #endif

            std::vector<TH2*> decays_strip_xy;
            std::vector<TH2*> decays_pos_xy;
            std::vector<TH1*> decays_e;
            std::vector<TH2*> decays_e_xy;
            std::vector<TH1*> decays_time_delta;
            std::vector<TH1*> decays_strip_1d;
            std::vector<TH1*> decays_per_event;
            std::vector<TH1*> decays_channels;

             std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
            AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);
            
            int      IsData(std::ifstream &f);

          

            //WR histograms
            TH1 *hAida_Fat_WRdT;
            TH1 *hAida_Gal_WRdT;
            TH1 *hAida_bPlas_WRdT;
            TH1 *hbPlas_Fat_WRdT;
            TH1 *hbPlas_Gal_WRdT;
            TH1 *hFat_Gal_WRdT;
            TH1 *hFRS_Gal_WRdT;
            TH1 *hFRS_bPlas_WRdT;

            //FRS Histograms
             TH2I *hID_x2AoQ;
             TH2I *hID_x4AoQ;
             TH2I *hID_x4AoQ_zsame;
             TH2I *hID_Z_AoQ;
             TH2I *hID_Z_AoQ_zsame;
             TH2I *hID_Z_AoQ_corr;
             TH2I *hID_Z_AoQ_overlap;
             TH2I *hID_Z_Z2;
             TH2I *hID_Z_Z2_nocorr;
             TH2I *hdEdegoQ_Z;
             TH2I *hdEdegoQtimesBeta_Z;
             TH2I *hdEdeg_Z;
             TH2I *hdEdegoQ_Beta;
             TH2I *hdEdegoQ_AoQ;
             TH2I *hdEdegoQtimesBeta_AoQ;
              
             TH2I *hID_x2AoQ_x2AoQgate[8];
             TH2I *hID_x2AoQ_x4AoQgate[8];
             
             TH2I *hID_x4AoQ_x2AoQgate[8];
             TH2I *hID_x4AoQ_x4AoQgate[8];
             
             TH2I *hID_ZAoQ_x2AoQgate[8];
             TH2I *hID_ZAoQ_x4AoQgate[8];
             
             TH2I *hID_x2AoQ_Z1Z2gate[8];
             TH2I *hID_x4AoQ_Z1Z2gate[8];
             TH2I *hID_ZAoQ_Z1Z2gate[8];
             TH2I *hID_SC43Z1_Z1Z2gate[8];
             
             TH2I *hID_Z1_AoQ_dEdegZgate[8];
             TH2I *hID_Z1_AoQ_zsame_dEdegZgate[8];
             TH2I *hID_Z1_AoQcorr_dEdegZgate[8];
             TH2I *hID_Z1_AoQcorr_zsame_dEdegZgate[8];
             
             TGo4PolyCond  *cID_Z_AoQ[8];       
             TGo4PolyCond  *cID_Z_Z2gate[8];
             TGo4PolyCond  *cID_x2AoQ[8];
             TGo4PolyCond  *cID_x4AoQ[8];
             TGo4PolyCond  *cID_dEdeg_Z1[8];
        
             //TH1I *hID_Z3_gate[5];
             TH2I *hID_Z_AoQgate[8];  
             TH2I *hID_x2AoQ_ZAoQgate[8];  
             TH2I *hID_x4AoQ_ZAoQgate[8];  
             TH2I *hID_Z_Z2_ZAoQgate[8];  
             
             
            //bPlast Histograms           
            TH1 *hbPlas_ToT_det[3][32];       
            TH1 *hbPlas_Lead_T[3][16];
            TH1 *hbPlas_Trail_T[3][16];
            TH1 *hbPlas_Multiplicity_Det1;
            TH1 *hbPlas_Multiplicity_Det2; 
           
            //TH1 *hbPlas_lead_lead[3][32];
            TH1 *hbPlas_lead_lead_ref_det1[3][16];
            TH1 *hbPlas_lead_lead_ref_det2[3][16];
            TH1 *hbPlas_lead_lead_gated[3][32];       
            TH1 *hbPlas_ToT[3][32];
//             TH1 *hbPlas_ToT_coinc[3][32];
//             TH1 *hbPlas_ToT_coinc_dets[3][32];
//             TH1 *hbPlas_ToT_coinc_dets_perchan[3][32];
            TH1 *hbPlas_trail_trail[3][32];
            TH1 *hbPlas_Energy_Calib[3][32];
             TH1 *hbPlas_SC41L_Anal_lead[3][32];
             TH1 *hbPlas_SC41R_Anal_lead[3][32];
             TH1 *hbPlas_SC41L_Digi_lead[3][32];
             TH1 *hbPlas_SC41R_Digi_lead[3][32];
             TH1 *hbPlas_ToT_Sum[3];
 
            TH1 *hbPlas_hit_pattern_det[3];
            TH1 *hbPlas_num_fired_chans;
            TH1 *hSC41_Analogue_Tamex;
            TH1 *hSC41_Digital_Tamex;

             TH1 *hFATTEMP_lead_lead_ref[32];
             TH1 *hFATTEMP_lead_lead[3][32];
             TH1 *hFATTEMP_lead_lead_gated[3][32];
             TH1 *hFATTEMP_trail_trail[3][32];
             TH1 *hFATTEMP_ToT[3][32];
             TH1 *hFATTEMP_Energy_Calib[3][32];
             TH1 *hFATTEMP_SC41L_Anal_lead[3][32];
             TH1 *hFATTEMP_SC41R_Anal_lead[3][32];
             TH1 *hFATTEMP_SC41L_Digi_lead[3][32];
             TH1 *hFATTEMP_SC41R_Digi_lead[3][32];

             TH1 *hFat_Lead_T[48];
             TH1 *hFat_Trail_T[48];
             TH1 *hFat_lead_lead_ref[48];
             TH1 *hFat_ToT_det[48];
             TH1 *hFat_ToT_Sum;
             TH1 *hFat_tamex_hit_pattern;
             TH1 *hFat_tamex_multiplicity;


            TH1 *hScalar_hit_pattern;
            //Fatima Histograms
            TH1 *hFAT_Energy[50];
            TH1 *hFAT_QDCdt[50];
            TH1 *hFAT_EnergySum;
            TH1 *hFAT_hits_QDC;
            TH1 *hFAT_TDCdt_refCha[50];
            
            TH1 *hFAT_Multipl;
            TH1 *hFAT_hits_TDC;

            TH1 *hFAT_SC41_check;

            TH1 *hFAT_lead_lead[50];
            TH1 *hFAT_lead_lead_ref[50];
            TH1 *hFAT_lead_lead_gated[50];
            TH1 *hFAT_lead_lead_gated1[50];
            TH1 *hFAT_lead_lead_gated2[50];
            TH2 *hFAT_lead_lead_energy[50];
            TH1 *hFAT_ToT[50];
            TH1 *hFAT_trail_trail[50];
            TH2 *hFAT_gamma_gamma;

            TH1* hFAT_lead_lead_QDC_Gate[50];

            ///Galileo Histograms
            TH1 *hGAL_Chan_E[GALILEO_MAX_DETS][GALILEO_CRYSTALS];
            //TH1 *hGAL_Chan_E2;
            TH1 *hGAL_Chan_Egate;
            TH1 *hGAL_AddbackSum;
            TH1 *hGAL_ESum_SC41;
            TH1 *hGAL_Chan_Time_Diff[32];
            TH1 *hGAL_Time_Diff_vs_Energy[32];
            TH1 *hGAL_ESum;
            TH1 *hGAL_ESum_largerange_OF;
            TH1 *hGAL_ESum_largerange_PU;
            TH1 *hGAL_Hit_Pat;
            TH2 *hGAL_Chan_E_Mat;
            TH2 *hGAL_Chan_E_vsCrys;
            TH2 *hGAL_Time_Diff_vs_Energy_sum;
            TH2 *hGAL_Chan_E_vsDet;
             TH1 *hGAL_Chan_E_gate[32];
             TH1 *hGAL_ge_time_difference_gg;
        

            //Finger Histograms
            TH1 *hFING_Hit_Pat;
//             TH1 *hFING_lead_lead[52];
//             TH1 *hFING_trail_trail[52];
//             TH1 *hFING_ToT[52];
//             TH1 *hFING_MaxToT[52];
//             TH1 *hFING_ToTMaxAdded[52];
//             TH1 *hFING_trig_lead[52];
//             TH1 *hFING_SC41_ToT;
//             TH2 *hFING_ToT_StripID;
//             TH2 *hFING_ToT_StripID_UpDown;
//             TH2 *hFING_ToT_StripID_UpDown_max;
//             TH2 *hFING_Pos_ToT;
//             TH2 *hFING_Pos_ToT_Max;
//             TH1 *hFING_Pos_ToT_gated;
//             TH2 *hFING_MaxToT_StripID;
//             TH2 *hFING_MaxToTExp_StripID;
//             TH2 *hFING_ToT_lead_lead[52];
//             TH2 *hFING_Pos;                 //E.Sahin 24.05
//             TH2 *hFING_ToT_StripID_Exp;
//             TH1 *hFING_fcoarse[52];
//             TH1 *hFING_ffine[52];
//             TH2* hFING_LeadLead_StripID;
            TH1 *hFING_Multiplicity;
            TH1 *hFING_ToT_Strip[52];
            TH1 *hFING_MaxToT_Strip[52];
            TH2 *hFING_ToT_v_Strip;
            TH2 *hFING_MaxToT_v_Strip;

            TH1 *hFING_ToT_PMT[52];
            TH2 *hFING_ToT_v_PMT;
            TH1 *hFING_lead_lead[52];
            TH1 *hFING_trail_trail[52];
            TH1 *hFING_Sc41lead_leadmaxtot[52];
            TH2* hFing_ToTRatio_v_Strip;
            TH2* hFing_MaxToTRatio_v_Strip;

            TH1 *hFING_SC41_lead_lead;
            TH1 *hFING_SC41_trail_trail;
            TH1 *hFING_SC41_tot;
            TH1 *hFat_minus_plasticTDC[50];
            TH1 *hFat_minus_plasticTDCFatGate[50];
            TH1 *hFat_bplas_Corr_RawTDCFAT_RawTDCbP_Dets[50];

            TH1 *hFat_bplas_Corr_SC41_Dets[50];
            TH1 *hFat_bplas_Corr_SC41_Dets_gated[50];
            TH1 *hFat_bplas_Corr_SiPM_Dets[50];
            TH1 *hFat_bplas_Corr_SiPM_Dets_gated[50];
            TH2 *hFat_bplas_Corr_EFatEbPlas;
            TH1 *hFat_bplas_Corr_SiPM_Gated_Dets[50];
            TH1 *hSC41FatPlas;

            TH2 *h_FingStrip_PlasID;
            TH2 *h_FingToT_PlasE;
            TH2 *h_FingToT_PlasT[32];
            TH1 *h_FingToT_PlasQDC;
            TH2 *hFing_MaxToT_v_Strip_bPgated;

            TGo4Picture *picPMT;
private :
            bool FAT_dist_corr_used; // Read from General Setup File
            int  FAT_exclusion_dist; // Read from General Setup File
            int  FAT_num_TDC_modules; // Read from General Setup File
            bool FAT_nearest_neighbour_exclusion; // Read from General Setup File
            bool same_ring_exclusion; // Read from General Setup File
            bool output_position_matrix; // Read from General Setup File

            bool VME_TAMEX_bPlas,VME_TAMEX_Fatima, VME_AND_TAMEX_Fatima;
            void clear_gate_arrays();

      ClassDef(EventAnlProc, 1)
	};
#endif //EVENTANLPROCESSOR_H