// Adapted for DESPEC by A.K.Mistry 2020
//
//---------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//       Experiment Data Processing at EE department, GSI
//---------------------------------------------------------------
//
//Copyright (C) 2000- Gesellschaft f. Schwerionenforschung, GSI
//                    Planckstr. 1, 64291 Darmstadt, Germany
//Contact:            http://go4.gsi.de
//----------------------------------------------------------------
//This software can be used under the license agreements as stated
//in Go4License.txt file which is part of the distribution.
//----------------------------------------------------------------

#ifndef UNPACKEVENTSTORE_H
#define UNPACKEVENTSTORE_H

#include "AIDA_Event.h"
#include "TGo4EventElement.h"

struct AidaUnpackData {
     uint64_t AIDATime;
     int AIDAHits;
       
      std::vector<AidaEvent> ImplantEvents;
      std::vector<AidaEvent> DecayEvents;

      std::vector<AidaHit> Implants;
      std::vector<AidaHit> Decays;
};

class EventUnpackStore : public TGo4EventElement {
public:
        EventUnpackStore();
        EventUnpackStore(const char* name);
        virtual ~EventUnpackStore();

        /** Method called by the framework to clear the event element. */
        void Clear(Option_t *t="");
        //Clear the members
        //void ClearEvent();       

      
        //FRS Output
   Float_t fFRS_Music_dE[3], fFRS_Music_dE_corr[3];
   Float_t fFRS_Music_E1[8], fFRS_Music_E2[8], fFRS_Music_T1[8], fFRS_Music_T2[8];
   Float_t fFRS_sci_l[12], fFRS_sci_r[12], fFRS_sci_e[12], fFRS_sci_tx[12], fFRS_sci_x[12];
   
   Double_t fTRaw_vftx_21l,fTRaw_vftx_21r,fTRaw_vftx_22l,fTRaw_vftx_22r,fTRaw_vftx_41l,fTRaw_vftx_41r,fTRaw_vftx_42l,fTRaw_vftx_42r;
  //  Float_t fFRS_sci_e[12];
  //  Float_t fFRS_sci_tofll2, fFRS_sci_tofll3, fFRS_sci_tof2, fFRS_sci_tofrr2, fFRS_sci_tofrr3, fFRS_sci_tof3;
    Float_t fFRS_ID_x2, fFRS_ID_y2, fFRS_ID_a2, fFRS_ID_b2;
    Float_t fFRS_ID_x4, fFRS_ID_y4, fFRS_ID_a4, fFRS_ID_b4;
   
  //  Int_t fFRS_sci_dt_21l_21r, fFRS_sci_dt_41l_41r, fFRS_sci_dt_42l_42r, fFRS_sci_dt_43l_43r;
   // Int_t fFRS_sci_dt_21l_41l, fFRS_sci_dt_21r_41r, fFRS_sci_dt_21l_42l, fFRS_sci_dt_21r_42r;
    Float_t fFRS_ID_brho[2];
    //, fFRS_ID_rho[2];
    Float_t fFRS_beta, fFRS_gamma;
    Float_t fFRS_AoQ, fFRS_AoQ_corr;
    Float_t fFRS_z, fFRS_z2, fFRS_znocorr,fFRS_z2nocorr;
    Float_t fFRS_dEdeg, fFRS_dEdegoQ;
  //  Float_t fFRS_timestamp, fFRS_ts, fFRS_ts2;  
    Long64_t fFRS_WR;
    Float_t fFRS_TPC_x[7];
    Float_t fFRS_TPC_y[7];
    UInt_t fFRS_scaler[64];
    UInt_t fFRS_scaler_delta[64];
        //AIDA output 
        std::vector<AidaUnpackData> Aida;
        
         int fAIDAHits;
         int AIDAHits;
         Long64_t AIDATime; 
         Long64_t fAIDA_WR; 
       
        //General Output
        Int_t  fevent_number;
        Int_t  fTrigger;
        Int_t  fProcID[7];
        Int_t    fScalar_fired;
        Int_t    fScalar_ID;
        
         UInt_t fFat_mult;
         
        UInt_t fFat_SC40mult;
        UInt_t fFat_SC41mult; 
        UInt_t fFat_bPlasmult;
        UInt_t fFat_QDC_Multiplicity;
        Int_t  fFat_QDC_ID[50];
        Double_t fFat_QDC_E[50];
        Double_t fFat_QDC_E_Raw[50];
        Long64_t fFat_QDC_T_coarse[50];
        Double_t fFat_QDC_T_fine[50];
    
        Int_t fFat_TDC_ID[50];
        UInt_t fFat_TDC_Multiplicity[50];
        Long64_t fFat_WR;      
        double_t fSC40[10];
        double_t fSC41[10];
        double_t fFat_bPlas[100];
        Long64_t fFat_TDC_Time[50];
        Long64_t fFat_TDC_Time_Raw [50];
        
        Int_t fFat_TDC_Singles_ID[50];
        Long64_t fFat_TDC_Singles_t[50];
        Long64_t fFat_TDC_Singles_t_Raw[50]; 
        Int_t fFat_QDC_Singles_ID[50];
        Double_t fFat_QDC_Singles_E[50];
        Double_t fFat_QDC_Singles_E_Raw[50];
        Long64_t fFat_QDC_Singles_t_coarse[50];
        Double_t fFat_QDC_Singles_t_fine[50];
        UInt_t fFat_tdcsinglescount;
        UInt_t fFat_qdcsinglescount; 
        
        ///Germanium
         Int_t fGal_fired;
         Int_t fGal_Detector[GALILEO_MAX_HITS];
         Int_t fGal_Crystal[GALILEO_MAX_HITS];
         Double_t fGal_E[GALILEO_MAX_HITS];
         ULong64_t fGal_T[GALILEO_MAX_HITS];
         ULong64_t fGal_Event_T[GALILEO_MAX_HITS];
         bool fGal_Pileup[GALILEO_MAX_HITS];
         bool fGal_Overflow[GALILEO_MAX_HITS];
         ULong64_t fGal_WR;
         
         ///Finger
        Long64_t fFinger_WR;
       Int_t fFing_Strip_N[52];
       Double_t fFing_Lead_Up[52][100];
       Double_t fFing_Trail_Up[52][100];
       Double_t fFing_Lead_Down[52][100];
       Double_t fFing_Trail_Down[52][100];
       Int_t fFing_Strip_N_LU[52];
       Int_t fFing_Strip_N_TU[52];
       Int_t fFing_Strip_N_LD[52];
       Int_t fFing_Strip_N_TD[52];
       Double_t fFing_SC41_lead[2][50];
       Double_t fFing_SC41_trail[2][50];
       Int_t fFing_PMT_Lead_N[52];
       Int_t fFing_PMT_Trail_N[52];
       Double_t fFing_Lead_PMT[52][10];
       Double_t fFing_Trail_PMT[52][10];
       
       //FATIMA TAMEX 
      // Int_t fFat_Strip_N[50];     
//        Double_t fFat_SC41_lead[2][50];
//        Double_t fFat_SC41_trail[2][50];
       Int_t fFat_chan;
       Int_t fFat_PMT_Lead_N[50];
       Int_t fFat_PMT_Trail_N[50];
       Double_t fFat_Lead_PMT[50][10];
       Double_t fFat_Trail_PMT[50][10];
       Long64_t fFat_Tamex_WR;
       
       
        //bPlas TAMEX 
       Long64_t fbPlas_WR;
       Int_t fbPlaschan;
       Int_t fbPlas_PMT_Lead_N[3][16];
       Int_t fbPlas_PMT_Trail_N[3][16];
       Double_t fbPlas_Lead_PMT[3][16][10];
       Double_t fbPlas_Trail_PMT[3][16][10];

   ClassDef(EventUnpackStore,1)
};
#endif //TEVENT_H




//----------------------------END OF GO4 SOURCE FILE ---------------------
