// $Id: EventAnlProc.cxx 754 2011-05-18 11:04:52Z adamczew $
// Adapted for DESPEC by A.K.Mistry 2020
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum fuer Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
    //-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

// Uncomment this to align the AIDA ASICs with a pulser
//  Only needed if the ASICs didn't align properly
//#define AIDA_PULSER_ALIGN

#include "EventAnlProc.h"

#include <cstdlib>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"

#include "TGo4WinCond.h"
#include "TGo4Analysis.h"

#include "TGo4Picture.h"

#include "EventAnlStore.h"
#include "EventUnpackStore.h"

#include "TAidaConfiguration.h"
#define ABS(x)  ((x)>=0 ? (x):-(x))  // absolute_value(x)
//-----------------------------------------------------------
EventAnlProc::EventAnlProc() :
TGo4EventProcessor()
{
}
//-----------------------------------------------------------
EventAnlProc::EventAnlProc(const char* name) :
TGo4EventProcessor(name)
{
  //Clear up for AIDA
  implantEvents = 0;
  goodImplantEvents = 0;
  stoppedEvents = 0;
  decayEvents = 0;
  pulserEvents = 0;
  nonsenseEvents = 0;

  cout << "**** EventAnlProc: Create" << endl;

  checkTAMEXorVME();

  fCal = (CalibParameter*) GetParameter("CalibPar"); 

  ///FINGER 2D Polygon ToT vs strip gate window test
    Double_t ToTvalues[15]={17,24,30,36,42,48,52,52,48,42,36,30,24,17,17};
    Double_t Stripvalues[15]={0,0,0,0,0,0,0,1,1,1,1,1,1,1,0};
    TCutG* ToTvsStripcut = new TCutG("initialcut",15,ToTvalues,Stripvalues);
    fCond_FingToTvsStrip = new TGo4PolyCond("FING_TOTvsStrip");
    fCond_FingToTvsStrip -> SetValues(ToTvsStripcut);
    AddAnalysisCondition(fCond_FingToTvsStrip);
    fCond_FingToTvsStrip -> Enable();
    delete ToTvsStripcut;

   
     
     
    clear_gate_arrays();
    read_setup_parameters();
    get_used_systems();
    FRS_Gates();
    
   

}
//-----------------------------------------------------------
EventAnlProc::~EventAnlProc()
{
  cout << "**** EventAnlProc: Delete" << endl;

}

void EventAnlProc::UserPostLoop()
{
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  if (!conf->ShowStats()) return;
  std::cout << "AIDA Analysis Statistics" << std::endl;
  std::cout << "Incoming Implant Events: " << implantEvents << std::endl;
  std::cout << "Good Implant Events    : " << goodImplantEvents << " (" << (100. * goodImplantEvents / implantEvents) << "%)" << std::endl;
  std::cout << "Stopped Implant Events : " << stoppedEvents << " (" << (100. * stoppedEvents / implantEvents) << "%)" << std::endl;
}

//-----------------------------------------------------------


Bool_t EventAnlProc::BuildEvent(TGo4EventElement* dest)
{
     for(int i=0; i<7; i++){
        PrcID_Conv[i]=-1;
        }

  Bool_t isValid=kFALSE; // validity of output event

  EventUnpackStore* pInput  = (EventUnpackStore*) GetInputEvent();
  EventAnlStore* pOutput = (EventAnlStore*) dest;

//   pAida.Implants.clear();
//   pAida.Decays.clear();


  if((pInput==0) || !pInput->IsValid()){ // input invalid
    pOutput->SetValid(isValid); // invalid
    return isValid; // must be same is for SetValid
  }
  isValid=kTRUE;

   ///general inputs from the unpacker
    event_number = pInput->fevent_number;
    pOutput->pOnSpill = FRS_spill;
    pOutput->pTrigger = pInput->fTrigger;
    pOutput->pEvent_Number = event_number;
    VMEorTAMEX_bPlas = VME_TAMEX_bPlas;
    VMEorTAMEX_fatima = VME_TAMEX_Fatima;
    VMEandTAMEX_fatima = VME_AND_TAMEX_Fatima;


  for (int i = 0; i<7; i++){
      if(pInput->fProcID[i]>-1){
      PrcID_Conv[i] = pInput->fProcID[i];
      pOutput->pPrcID_Conv[i] = pInput->fProcID[i];

    }
  }
 

   static bool create =false;
  //Create histograms
  if (!create)
  {
    Make_WR_Histos();
    if(Used_Systems[0])  Make_FRS_Histos();
    if(Used_Systems[1]) Make_Aida_Histos();
   // if (Used_Systems[2] && VMEorTAMEX_bPlas==true) Make_Plastic_VME_Histos();
    if(Used_Systems[2] && VMEorTAMEX_bPlas==false){ 
        Make_Plastic_Tamex_Histos();
        Make_Fatima_Tamex_Histos();
        
    }
    if (Used_Systems[3] && VMEorTAMEX_fatima==true) Make_Fatima_Histos();
   // if (Used_Systems[3] && VMEorTAMEX_fatima==false && VMEandTAMEX_fatima==false) Make_Fatima_Tamex_Histos();

   // if ((Used_Systems[3] || Used_Systems[4]) && VMEorTAMEX_fatima==false && VMEandTAMEX_fatima==true) Make_Fatima_VME_Tamex_Histos();


    if (Used_Systems[5]) Make_Galileo_Histos();
    //if (Used_Systems[6]) Make_Finger_Histos();
//     if (Used_Systems[2] && Used_Systems[3]) Make_Fat_Plas_Histos();
//     if (Used_Systems[6] && Used_Systems[2]) Make_Fing_Plas_Histos();
       
        }
        
        create = true;
        Do_WR_Histos(pInput);
                 /** Now extract the data from the stored Unpacker array (root tree)**/
    ///--------------------------------------/**FRS Input**/------------------------------------------///

      if (PrcID_Conv[0]==0&& Used_Systems[0]==1){
         FRS_AoQ=-999;
         FRS_z=-999;
         FRS_z2=-999;
          pOutput->pFRS_WR = pInput->fFRS_WR;
        ///MUSIC
       for(int i =0; i<2; ++i){
            FRS_Music_dE[i] = pInput->fFRS_Music_dE[i];
            FRS_Music_dE_corr[i] = pInput->fFRS_Music_dE_corr[i];
           }
            for(int i =0; i<8; ++i){
           FRS_Music_E1[i] = pInput->fFRS_Music_E1[i];
           FRS_Music_E2[i] = pInput->fFRS_Music_E2[i];
           FRS_Music_T1[i] = pInput->fFRS_Music_T1[i];
           FRS_Music_T2[i] = pInput->fFRS_Music_T2[i];
            }
           
        ///SCI
        
        for(int l=0;l<12;++l){
          
           FRS_sci_l[l] = pInput->fFRS_sci_l[l];
          
           FRS_sci_r[l] = pInput->fFRS_sci_r[l];
           FRS_sci_e[l] = pInput->fFRS_sci_e[l];
           FRS_sci_tx[l] = pInput->fFRS_sci_tx[l];
           FRS_sci_x[l] = pInput->fFRS_sci_x[l];
         
        }
        
        ///VFTX
        pOutput->pTRaw_vftx_21l = pInput->fTRaw_vftx_21l;
        pOutput->pTRaw_vftx_21r = pInput->fTRaw_vftx_21r;
        pOutput->pTRaw_vftx_22l = pInput->fTRaw_vftx_22l;
        pOutput->pTRaw_vftx_22r = pInput->fTRaw_vftx_22r;
        pOutput->pTRaw_vftx_41l = pInput->fTRaw_vftx_41l;
        pOutput->pTRaw_vftx_41r = pInput->fTRaw_vftx_41r;
        pOutput->pTRaw_vftx_42l = pInput->fTRaw_vftx_42l;
        pOutput->pTRaw_vftx_42r = pInput->fTRaw_vftx_42r;
        
          for(int brho = 0; brho < 2; brho++) 
              FRS_ID_brho[brho] = pInput->fFRS_ID_brho[brho];
        
        
        for(int i = 0; i < 7; i++)
        {
         FRS_tpc_x[i] = pInput->fFRS_TPC_x[i];
         FRS_tpc_y[i] = pInput->fFRS_TPC_y[i];
        
        }
        for (int i = 0; i < 64; i++)
        {
          pOutput->pFRS_scaler[i] = pInput->fFRS_scaler[i];
          pOutput->pFRS_scaler_delta[i] = pInput->fFRS_scaler_delta[i];
        }
        if (pOutput->pFRS_scaler_delta[8] > 0)
        {
          FRS_spill = true;
          pOutput->pOnSpill = true;
        }
        if (pOutput->pFRS_scaler_delta[9] > 0)
        {
          FRS_spill = false;
          pOutput->pOnSpill = false;
        }
        ///SCI TOF
       // FRS_sci_tofll2 = pInput->fFRS_sci_tofll2;
       // FRS_sci_tofll3 = pInput->fFRS_sci_tofll3;
        //FRS_sci_tof2 = pInput->fFRS_sci_tof2;
       // FRS_sci_tofrr2 = pInput->fFRS_sci_tofrr2;
        //FRS_sci_tofrr3 = pInput->fFRS_sci_tofrr3;
      //  FRS_sci_tof3 = pInput->fFRS_sci_tof3;
        ///ID 2 4
        FRS_ID_x2 = pInput->fFRS_ID_x2;
        FRS_ID_y2 = pInput->fFRS_ID_y2;
        FRS_ID_a2 = pInput->fFRS_ID_a2;
        FRS_ID_b2 = pInput->fFRS_ID_b2;

        FRS_ID_x4 = pInput->fFRS_ID_x4;
        FRS_ID_y4 = pInput->fFRS_ID_y4;
        FRS_ID_a4 = pInput->fFRS_ID_a4;
        FRS_ID_b4 = pInput->fFRS_ID_b4;
            ///SCI dT
//         FRS_sci_dt_21l_21r = pInput->fFRS_sci_dt_21l_21r;
//         FRS_sci_dt_41l_41r = pInput->fFRS_sci_dt_41l_41r;
//         FRS_sci_dt_42l_42r = pInput->fFRS_sci_dt_42l_42r;
//         FRS_sci_dt_43l_43r = pInput->fFRS_sci_dt_43l_43r;
//
//         FRS_sci_dt_21l_41l = pInput->fFRS_sci_dt_21l_41l;
//         FRS_sci_dt_21r_41r = pInput->fFRS_sci_dt_21r_41r;
//
//         FRS_sci_dt_21l_42l = pInput->fFRS_sci_dt_21l_42l;
//         FRS_sci_dt_21r_42r = pInput->fFRS_sci_dt_21r_42r;
            ///ID B Rho
//         for(int i =0; i<2; ++i){
//        // FRS_ID_brho[i] = pInput->fFRS_ID_brho[i];
//       //  FRS_ID_rho[i] = pInput->fFRS_ID_rho[i];
//         }
        FRS_beta = pInput->fFRS_beta;
//         FRS_beta3 = pInput->fFRS_beta3;
        FRS_gamma  = pInput->fFRS_gamma;
            ///ID Z AoQ
        FRS_AoQ = pInput->fFRS_AoQ;
        FRS_AoQ_corr = pInput->fFRS_AoQ_corr;
//         if(FRS_AoQ_corr>0){

//            // cout<<" " << endl;
//         }
        FRS_z = pInput->fFRS_z;
      // if(FRS_AoQ>2 && FRS_AoQ < 2.12) cout<<"ANL STAGE FRS_AoQ " <<pInput->fevent_number <<" AoQ  " << FRS_AoQ <<" FRS_z " <<FRS_z <<  endl;
        FRS_z2 = pInput->fFRS_z2;
        FRS_znocorr = pInput->fFRS_znocorr;
        FRS_z2nocorr = pInput->fFRS_z2nocorr;
        FRS_dEdeg = pInput->fFRS_dEdeg;
        FRS_dEdegoQ = pInput->fFRS_dEdegoQ;
       
       // FRS_z3 = pInput->fFRS_z3;
            ///ID Timestamp
       // FRS_timestamp = pInput->fFRS_timestamp;
//         FRS_ts = pInput->fFRS_ts;
//         FRS_ts2 = pInput->fFRS_ts2;
 // if(pOutput->pEvent_Number==42715 ){       }
 Do_FRS_Histos(pOutput);
            }


   ///-------------------------------- /**AIDA Input**/ --------------------------------///
      //  if (Used_Systems[1]&&  PrcID_Conv[1]==1) {ProcessAida(pInput);}
        ProcessAida(pInput, pOutput);
        Aida_Fired = 0;
        Aida_Fired = pInput->fAIDAHits;

   ///-------------------------------- /**bPlastic VME Input**/ --------------------------------///
        ///Disabled A.M. 11.12.19
       
///--------------------------------------/**bPlastic TAMEX Input**/------------------------------------------///
        
       
   if (PrcID_Conv[2] ==2){
      
        pOutput->pbPLAS_WR = pInput->fbPlas_WR;

       
       // for (int i = 0; i < fbPlaschan; i++){
            //for(int j=0; j<fOutput->fbPlas_PMT_Lead_N[2][chan];j++){ ///Hits iterator
                bPlas_TAM_SC41L_ANA = pInput->fbPlas_Lead_PMT[2][0][0];
                bPlas_TAM_SC41R_ANA = pInput->fbPlas_Lead_PMT[2][1][0];
                bPlas_TAM_SC41L_DIG = pInput->fbPlas_Lead_PMT[2][2][0];
                bPlas_TAM_SC41R_DIG = pInput->fbPlas_Lead_PMT[2][3][0];

       // bPlas_AND_Coinc[j] = pInput->fFat_Lead_PMT[9][j];
       
           // }
       // }
  // if(pOutput->pEvent_Number==100598) 
    Do_Plastic_Tamex_Histos(pInput,pOutput);
    Do_Fatima_Tamex_Histos(pInput,pOutput); 
   }
  
  ///--------------------------------------/**Fatima TAMEX Input**/------------------------------------------///
//    if (PrcID_Conv[4] ==4){
//        Do_Fatima_Tamex_Histos(pInput,pOutput);  
//         pOutput->pFAT_Tamex_WR = pInput->fFat_Tamex_WR;
//         
//         Fat_TAM_SC41L_ANA = pInput->fFat_Lead_PMT[48][0];
//         Fat_TAM_SC41R_ANA = pInput->fFat_Lead_PMT[49][0];
//         Fat_TAM_SC41L_DIG = pInput->fFat_Lead_PMT[50][0];
//         Fat_TAM_SC41R_DIG = pInput->fFat_Lead_PMT[51][0];
//    }

   ///--------------------------------------/**Fatima VME Input**/------------------------------------------///
  
        SC40mult = 0;
        SC41mult = 0;
        Fatmult = 0;
        Fat_bPlasmult = 0;
        Fat_Cha_Ref_TDC = 0;
        Fat_WR = 0;
        stdcmult = 0;
        sqdcmult = 0;
        
       
      for (int i=0; i<100; i++)  Fat_bPlas[i]=0;
      for (int i=0; i<50; i++){
         //Fat_QDC_E[i] = 0;
         //Fat_QDC_T[i] =0;
         SC40[i] = 0;
         SC41[i] = 0;
         
         Fat_TDC_ID[i] = -2;
         Fat_QDC_ID[i] = -1;
         Fat_TDC_T[i] = 0;
         
         Fat_TDC_Singles_ID[i] = -2;
         Fat_QDC_Singles_ID[i] = -1;
         Fat_QDC_Singles_E[i] = 0;
         Fat_QDC_Singles_t_coarse[i] = 0;
         Fat_QDC_Singles_t_fine[i] = 0;
         Fat_TDC_Singles_t[i] = 0;      
         
}///End of initialise arrays

  if ( PrcID_Conv[3]==3 && VMEorTAMEX_fatima==true){
      
    Fat_WR = pInput->fFat_WR;

    Fatmult =  pInput->fFat_mult;
    
    //cout << " anl proc fat mult " << Fatmult << endl; 
  
    
if(Fatmult > 0){

    //QDC
    for (int i = 0; i<Fatmult; i++){
      Fat_QDC_ID[i] = pInput->fFat_QDC_ID[i];
      Fat_QDC_E[i] = pInput->fFat_QDC_E[i];
      Fat_QDC_E_Raw[i] = pInput->fFat_QDC_E_Raw[i];
      Fat_QDC_T_coarse[i] = pInput->fFat_QDC_T_coarse[i];
      Fat_QDC_T_fine[i] = pInput->fFat_QDC_T_fine[i];
     }//End QDC for loop

    //TDC
      
    for (int i = 0; i<Fatmult; i++){
      Fat_TDC_ID[i] = pInput->fFat_TDC_ID[i];
      Fat_TDC_T[i] = pInput->fFat_TDC_Time[i];
      Fat_TDC_T_Raw[i] = pInput->fFat_TDC_Time_Raw[i];
       if(Fat_TDC_ID[i]==1)Fat_Cha_Ref_TDC = Fat_TDC_T[i]; 
       if(Fat_TDC_ID[i]==42) cout<<"bPlast is here" << endl;
    }//End TDC for loop
}
//                  INPUTS

      SC40mult =  pInput->fFat_SC40mult;
      SC41mult =  pInput->fFat_SC41mult; 
 
        for(int i = 0; i < SC40mult; i++){
            
        SC40[i] = pInput->fSC40[i]; 
            
        }//End sc41 for fill
        
        for(int i = 0; i < SC41mult; i++){
            
            SC41[i] = pInput->fSC41[i]; 
            
        }//End sc41 for fill
       
       Fat_bPlasmult = pInput->fFat_bPlasmult;
        for(int i = 0; i < Fat_bPlasmult; i++){
            
            Fat_bPlas[i] = pInput->fFat_bPlas[i]; 
            
        }//End bPlas Fatima coinc for fill
//     for (int i = 0; i<Fatmult; i++){
//     
//         //cout << " fat mult " << Fatmult << " i: " << i << " QDC ID: " << Fat_QDC_ID[i] << " TDC ID: " << Fat_TDC_ID[i] << " QDC Energy: " << FatQDC[i]  << " TDC Time: " << Fat_TDC_T[i] << endl;
//     }

		//Inputting singles TDC and QDC hits
        stdcmult = pInput->fFat_tdcsinglescount;
        sqdcmult = pInput->fFat_qdcsinglescount;

       for (int i = 0; i<stdcmult; i++){
        Fat_TDC_Singles_t[i] = pInput->fFat_TDC_Singles_t[i];
        Fat_TDC_Singles_t_Raw[i] = pInput->fFat_TDC_Singles_t_Raw[i];
        Fat_TDC_Singles_ID[i] = pInput->fFat_TDC_Singles_ID[i];
       }

       for (int i = 0; i<sqdcmult; i++){
        Fat_QDC_Singles_E[i] = pInput->fFat_QDC_Singles_E[i];
        Fat_QDC_Singles_E_Raw[i] = pInput->fFat_QDC_Singles_E_Raw[i];
        Fat_QDC_Singles_ID[i] = pInput->fFat_QDC_Singles_ID[i];
        Fat_QDC_Singles_t_coarse[i] = pInput->fFat_QDC_Singles_t_coarse[i];
        Fat_QDC_Singles_t_fine[i] = pInput->fFat_QDC_Singles_t_fine[i];
       }




        Do_Fatima_Histos(pOutput);
  }///End of proc ID


  ///--------------------------------------/**Galileo Input**/------------------------------------------///
       GalFired = -1;
       
       //Gal_WR = 0;
       for(int g = 0; g<GALILEO_MAX_HITS; g++){
          GalDet[g] = -1;
          GalCrys[g] = -1;
          GalE[g] = -1;
          GalE_Cal[g] = -1;
          GalT[g] =-1;
          GalPileUp[g] = false;
          GalOverFlow[g] = false;
       }
    
       if ( PrcID_Conv[5]==5){

        GalFired =  pInput->fGal_fired;
    //    GalPileup = pInput->fGal_Pileup;
        Gal_WR = pInput->fGal_WR;
    
        for (int i = 0; i<GalFired; i++)
        {
          GalDet[i] = pInput->fGal_Detector[i];
          GalCrys[i] = pInput->fGal_Crystal[i];
          GalPileUp[i] = pInput->fGal_Pileup[i];
          GalOverFlow[i] = pInput->fGal_Overflow[i];
          GalE[i] = pInput->fGal_E[i];
          GalEventT[i]=pInput->fGal_Event_T[i];
          GalT[i] = pInput->fGal_T[i];
          int id = GalDet[i] * 3 + GalCrys[i];
          GalE_Cal[i] = (fCal->AGal[id]* pow( GalE[i],2) + fCal->BGal[id]*  GalE[i] + fCal->CGal[id]);
        }

            Do_Galileo_Histos(pOutput);
      }

 ///--------------------------------------/**Finger Input**/------------------------------------------///
      


                        /** End of Unpack Tree input**/

  
  pOutput->SetValid(isValid);
  return isValid;

 }  //End of BuildEvent

///End of Input from Unpacker ///
 
///-----------------------------------------------------------------------------------------------------------------///
  void EventAnlProc::read_setup_parameters(){

    // unused // const char* format = "%s %d";

    ifstream file("Configuration_Files/DESPEC_General_Setup/Detector_System_Setup_File.txt");

    if(file.fail()){
        cerr << "Could not find File for setup parameters!" << endl;
        exit(0);
    }

    string line;
    string var_name;
    // unused //int dummy_var;
    //file.ignore(256,'GENERAL_CONFIGURATION');

    file.ignore(256,':');
    file >> FAT_exclusion_dist;//dummy_var;

    file.ignore(256,':');
    file >> FAT_nearest_neighbour_exclusion;//dummy_var;

    file.ignore(256,':');
    file >> same_ring_exclusion;//dummy_var;

    file.ignore(256,':');
    file >> output_position_matrix;//dummy_var;

    cout<<endl;
    cout<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Setup Parameters List Analysis Proc: "<<endl;
    if(FAT_exclusion_dist > 0) cout<<"FATIMA Detectors Excluded if Linear Difference Exceeds "<<FAT_exclusion_dist<<" mm"<<endl;
    else if(FAT_exclusion_dist == 0) cout<<"'Nearest Neighbour Exclusion': Disabled (Distance set to 0)"<<endl;
    cout<<"output_position_matrix " <<output_position_matrix<< endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<endl;
    cout<<endl;
    /*while(file.good()){
        getline(file,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&var_name,&dummy_var);

        cout<<"Hello Again?"<<endl;

        if (var_name == "White_Rabbit_Enabled:" && dummy_var == 1)  WHITE_RABBIT_USED = true;
        else if (var_name == "White_Rabbit_Enabled:" && dummy_var == 0)  WHITE_RABBIT_USED = false;

        if (var_name == "FATIMA_Gain_Match_Enabled:" && dummy_var == 1)  FAT_gain_match_used = true;
        else if (var_name == "FATIMA_Gain_Match_Enabled:" && dummy_var == 0) FAT_gain_match_used  = false;

    }*/

}

/**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------    White Rabbit   ---------------------------------------------**/
 /**----------------------------------------------------------------------------------------------**/
  void EventAnlProc::Make_WR_Histos(){
   hAida_Fat_WRdT = MakeTH1('I',"WR/Aida_Fatima_dT","White Rabbit Aida-Fatima",1000,-1000,1000);
   hAida_Gal_WRdT = MakeTH1('I',"WR/Aida_Galileo_dT","White Rabbit Aida-Galileo",1000,-1000,1000);
   hAida_bPlas_WRdT = MakeTH1('I',"WR/Aida_bPlas_dT","White Rabbit Aida-bPlas",1000,-1000,1000);
   hbPlas_Fat_WRdT = MakeTH1('I',"WR/bPlas_Fatima_dT","White Rabbit bPlas-Fatima",1000,-1000,1000);
   hbPlas_Gal_WRdT = MakeTH1('I',"WR/bPlas_Galileo_dT","White Rabbit bPlas-Galileo",1000,-1000,1000);
   hFat_Gal_WRdT = MakeTH1('I',"WR/Fatima_Galileo_dT","White Rabbit Fatima_Galileo",1000,-1000,1000);
   hFRS_Gal_WRdT = MakeTH1('I',"WR/FRS_Galileo_dT","White Rabbit FRS_Galileo",1000,-1000,1000);
   hFRS_bPlas_WRdT = MakeTH1('I',"WR/FRS_bPlas_dT","White Rabbit FRS_bPlas",1000,-1000,1000);



}
 void EventAnlProc::Do_WR_Histos(EventUnpackStore* pInput){
      
  if(pInput->fAIDA_WR>0 && pInput->fFat_WR>0)hAida_Fat_WRdT->Fill(pInput->fAIDA_WR -pInput->fFat_WR );
  if(pInput->fAIDA_WR>0 && pInput->fGal_WR>0)hAida_Gal_WRdT->Fill(pInput->fAIDA_WR -pInput->fGal_WR );
  if(pInput->fAIDA_WR>0 && pInput->fbPlas_WR>0)hAida_bPlas_WRdT->Fill(pInput->fAIDA_WR -pInput->fbPlas_WR);
  if(pInput->fFat_WR>0 && pInput->fbPlas_WR>0)hbPlas_Fat_WRdT->Fill(pInput->fbPlas_WR -pInput->fFat_WR);
  if(pInput->fGal_WR>0 && pInput->fbPlas_WR>0)hbPlas_Gal_WRdT->Fill(pInput->fbPlas_WR - pInput->fGal_WR);
  if(pInput->fGal_WR>0 && pInput->fFat_WR>0)hFat_Gal_WRdT->Fill(pInput->fFat_WR - pInput->fGal_WR);
    if(pInput->fGal_WR>0 && pInput->fFRS_WR>0)hFRS_Gal_WRdT->Fill(pInput->fFRS_WR - pInput->fGal_WR);
    if(pInput->fbPlas_WR>0 && pInput->fFRS_WR>0)hFRS_bPlas_WRdT->Fill(pInput->fFRS_WR - pInput->fbPlas_WR);

 }


/**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------    FRS   ---------------------------------------------**/
 /**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_FRS_Histos(){
    // hFRS_z1_z2 = MakeTH2('D',"FRS/z1vsz2","FRS Z1 vs Z2", 500, 0, 5000,  500, 0, 5000);
    hID_x2AoQ = MakeH2I("FRS/ID","ID_x2AoQ", 1500,1.5,2.5,500,-100,100.,"A/Q s2-s4", "X at S2 [mm]", 2);
    hID_x4AoQ = MakeH2I("FRS/ID","ID_x4AoQ", 1500,1.5,2.5,500,-100,100.,"A/Q s2-s4", "X at S4 [mm]", 2);
    
    hID_x4AoQ_zsame = MakeH2I("FRS/ID","ID_x4AoQ_zsame", 1500,1.5,2.5,500,-100,100,"A/Q s2-s4", "X at S4 [mm]", 2);

    hID_Z_AoQ = MakeH2I("FRS/ID","ID_Z_AoQ", 1500,1.9,2.7, 1500,25.,60.,"A/Q s2-s4", "Z s2-s4", 2);
    hID_Z_AoQ_zsame = MakeH2I("FRS/ID","ID_Z_AoQ_zsame", 1500,1.9,2.7, 1500,25.,60,"Z1==Z2 A/Q s2-s4", "Z s2-s4", 2);
    
    hID_Z_AoQ_corr = MakeH2I("FRS/ID","ID_Z_AoQ_S2-S4corr", 1500,1.3,2.8, 1500,30.,60., "A/Q s2-s4", "Z s2-s4", 2);
    
    
    hID_Z_AoQ_overlap = MakeH2I("FRS/ID","ID_Z_AoQ_overlap", 1500,1.9,2.7, 1500,25.,60.,"A/Q s2-s4", "Z s2-s4", 2);
    
    //High-Z charge state
    hdEdegoQ_Z = MakeH2I("FRS/ID/ID_S4","ID_dEdegoQ_Z", 1500,30,60, 1000, 0.1,0.8, "Z from MUSIC41", "dE(S2deg)/Q [a.u.]", 2);
    
    hdEdegoQtimesBeta_Z = MakeH2I("FRS/ID/ID_S4","ID_dEdegoQtimesBeta_Z", 1500,30,60, 1000, 0.1,0.8, "Z from MUSIC41", "(dE(S2deg)/Q)*Beta [a.u.]", 2);
    
    hdEdeg_Z   = MakeH2I("FRS/ID/ID_S4","ID_dEdeg_Z"  , 1500,30,60, 1000, 10.,50., "Z from MUSIC41", "dE(S2deg) [a.u.]", 2);
   
    hdEdegoQ_Beta = MakeH2I("FRS/ID/ID_S4","ID_dEdegoQ_Beta", 2000,0,1, 1000, 0.1,0.8, "Beta", "dE(S2deg)/Q [a.u.]", 2);
    
    hdEdegoQ_AoQ = MakeH2I("FRS/ID/ID_S4","ID_dEdegoQ_AoQ", 1500,1.3,2.8, 1000, 0.1,0.8, "AoQ", "dE(S2deg)/Q [a.u.]", 2);
    
    hdEdegoQtimesBeta_AoQ = MakeH2I("FRS/ID/ID_S4","ID_dEdegoQtimesBeta_AoQ", 1500,1.3,2.8, 1000, 0.1,0.8, "AoQ", "(dE(S2deg)/Q)*Beta [a.u.]", 2);
    
    hID_Z_Z2 = MakeH2I("FRS/ID","ID_Z_Z2", 2000,35,70, 2000,35.,70.,"Z", "Z2", 2);
    
    hID_Z_Z2_nocorr = MakeH2I("FRS/ID","ID_Z_Z2_nocorr", 2000,0,20000, 2000,0.,20000.,"Z", "Z2", 2);

    int num_ID_x2AoQ = {6};
    int num_ID_x4AoQ = {6};
    int num_ID_Z_Z2 = {6};
    int num_ID_Z_AoQ = {6};
    int num_ID_dEdeg_Z1{6};
     Float_t init_ID_x2AoQ[8][6][2];
     Float_t init_ID_x4AoQ[8][6][2];
     Float_t init_ID_Z_Z2[8][6][2];
     Float_t init_ID_Z_AoQ[8][8][2];
     Float_t init_dEdeg_Z1[8][6][2];
     ///FRS gates initialisation
   for(int i=0; i<8; i++){
       for(int j=0; j<6; j++){
        init_ID_x2AoQ[i][j][0] = XX2_AoQ[i][j];
        init_ID_x2AoQ[i][j][1] = YX2_AoQ[i][j];
        init_ID_x4AoQ[i][j][0] = XX4_AoQ[i][j];
        init_ID_x4AoQ[i][j][1] = YX4_AoQ[i][j];
        init_ID_Z_Z2[i][j][0] = X_ZZ2[i][j];
        init_ID_Z_Z2[i][j][1] = Y_ZZ2[i][j];
        
        init_dEdeg_Z1[i][j][0] =X_dEdeg[i][j];
        init_dEdeg_Z1[i][j][1] =Y_dEdeg[i][j];
        
            }
        for(int j=0; j<8; j++){
        init_ID_Z_AoQ[i][j][0] =X_ZAoQ[i][j];
        init_ID_Z_AoQ[i][j][1] =Y_ZAoQ[i][j];
       
    }
   }
   
  
    char name[50], title[100];
  ///Z vs AoQ
       
    for(int i=0; i<8; i++){
      sprintf(name,"cID_Z_AoQ%d",i);
      cID_Z_AoQ[i] = MakePolyCond("FRS_ID_Gated", name, num_ID_Z_AoQ, init_ID_Z_AoQ[i], hID_Z_AoQ->GetName());
     
      sprintf(name,"ID_Z_AoQgate%d",i);
      hID_Z_AoQgate[i] = MakeH2I("FRS/ID_Gated/ZAoQ_ZAoQgate",name,  1500,1.9,2.7, 1500,25.,60.,"A/Q s2-s4", "Z s2-s4", 2);
          
     sprintf(name,"ID_x2AoQ_ZAoQgate%d",i);
      hID_x2AoQ_ZAoQgate[i] = MakeH2I("FRS/ID_Gated/x2AoQ_ZAoQgate", name, 500,1.5,2.5,500,-100,100,"A/Q s2-s4", "gate on Z AoQ X at S2 [mm]", 2);

      sprintf(name,"ID_x4AoQ_ZAoQgate%d",i);
      hID_x4AoQ_ZAoQgate[i] = MakeH2I("FRS/ID_Gated/x4AoQ_ZAoQgate", name, 500,1.5,2.5,500,-100,100,"A/Q s2-s4", "gate on Z    X at S4 [mm]", 2);

      sprintf(name,"ID_ZAoQ_ZAoQgate%d",i);
      hID_Z_Z2_ZAoQgate[i] = MakeH2I("FRS/ID_Gated/Z1Z2_ZAoQgate", name, 2000,35,70, 2000,35.,70.,"Z music1", " Z music2", 2);
      
      
       ///Z vs Z2
        
      sprintf(name,"cID_Z_Z2%d",i);
      cID_Z_Z2gate[i] = MakePolyCond("FRS_Z1_Z2_Gated",name,num_ID_Z_Z2,init_ID_Z_Z2[i], hID_Z_Z2 ->GetName());
       
      sprintf(name,"ID_x2AoQ_Z1Z2gate%d",i);
      hID_x2AoQ_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated/x2AoQ_Z1Z2gate", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on Z    X at S2 [mm]", 2);

      sprintf(name,"ID_x4AoQ_Z1Z2gate%d",i);
      hID_x4AoQ_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated/x4AoQ_Z1Z2gate", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on Z    X at S4 [mm]", 2);

      sprintf(name,"ID_ZAoQ_Z1Z2gate%d",i);
      hID_ZAoQ_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated/ZAoQ_Z1Z2gate", name, 300,1.4,2.5, 400,1.,20.,"A/Q s2-s4", " Z music2", 2);

      sprintf(name,"ID_SC43Z1_Z1Z2gate%d",i);
      hID_SC43Z1_Z1Z2gate[i] = MakeH2I("FRS/ID_Gated/SC43Z1_Z1Z2gate", name, 300,1.,2.4, 400,30.,90.,"SC41 dE", " Z music1", 2);
        

              
///X2 vs AoQ
      sprintf(name,"cID_x2AoQ%d",i);
      cID_x2AoQ[i] = MakePolyCond("FRS_ID_Gated",name,num_ID_x2AoQ,init_ID_x2AoQ[i], hID_x2AoQ->GetName());
      
      sprintf(name,"hID_x2AoQ_x2AoQgate%d",i);
      hID_x2AoQ_x2AoQgate[i] = MakeH2I("FRS/ID_Gated/x2AoQ_x2AoQgate", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X2: X at S2 [mm]", 2);

      sprintf(name,"hID_x4AoQ_x2AoQgate%d",i);
      hID_x4AoQ_x2AoQgate[i] = MakeH2I("FRS/ID_Gated/x4AoQ_x2AoQgate", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X2: X at S4 [mm]", 2);

      ///X4 vs AoQ
      sprintf(name,"cID_x4AoQ%d",i);
      cID_x4AoQ[i] = MakePolyCond("FRS_ID_Gated",name,num_ID_x4AoQ,init_ID_x4AoQ[i], hID_x4AoQ->GetName());
      
      sprintf(name,"hID_x2AoQ_x4AoQgate%d",i);
      hID_x2AoQ_x4AoQgate[i] = MakeH2I("FRS/ID_Gated/x2AoQ_x4AoQgate", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X4: X at S2 [mm]", 2);

      sprintf(name,"hID_x4AoQ_x4AoQgate%d",i);
      hID_x4AoQ_x4AoQgate[i] = MakeH2I("FRS/ID_Gated/x4AoQ_x4AoQgate", name, 300,1.,2.4, 200,-100.,100.,"A/Q s2-s4", "gate on FRS AoQ, ID X4: X at S4 [mm]", 2);

      sprintf(name,"ID_ZAoQ_x2AoQgate%d",i);
      hID_ZAoQ_x2AoQgate[i] = MakeH2I("FRS/ID_Gated/ZAoQ_x2AoQgate", name, 300,1.,2.4, 400,30.,90.,"A/Q s2-s4", " Z music2", 2);

      sprintf(name,"ID_ZAoQ_x4AoQgate%d",i);
      hID_ZAoQ_x4AoQgate[i] = MakeH2I("FRS/ID_Gated/ZAoQ_x4AoQgate", name, 1500,1.9,2.7, 1500,25.,60.,"A/Q s2-s4", " Z music2", 2);
        
      ///dE_deg (high charge states)
      
      cID_dEdeg_Z1[i] = MakePolyCond("FRS_ID_Gated",name, num_ID_dEdeg_Z1,init_dEdeg_Z1[i], hdEdeg_Z->GetName());
       
      sprintf(name,"ID_Z1_AoQ_cdEdegZ %d",i);
      hID_Z1_AoQ_dEdegZgate[i]   = MakeH2I("FRS/ID_Gated/ZAoQ_IDdEdegZ1gate",name, 600,1.3,2.8, 600,30.,60.,"A/Q s2-s4", "Z from MUSIC41", 2);
      
      sprintf(name,"ID_Z1_AoQ_zsame_cdEdegZ %d",i);
      hID_Z1_AoQ_zsame_dEdegZgate[i]   = MakeH2I("FRS/ID_Gated/ZAoQSame_IDdEdegZ1gate",name, 600,1.3,2.8, 600,30.,60.,"A/Q s2-s4", "Z from MUSIC41", 2);
     
      sprintf(name,"ID_Z1_AoQcorr_cdEdegZ %d",i);
      hID_Z1_AoQcorr_dEdegZgate[i]   = MakeH2I("FRS/ID_Gated/ZAoQCorr_IDdEdegZ1gate",name, 600,1.3,2.8, 600,30.,60.,"A/Q s2-s4 corr.", "Z from MUSIC41", 2);
    
      sprintf(name,"ID_Z1_AoQcorr_zsame_cdEdegZ %d",i);
      hID_Z1_AoQcorr_zsame_dEdegZgate[i]   = MakeH2I("FRS/ID_Gated/ZAoQCorrSame_IDdEdegZ1gate",name, 600,1.3,2.8, 600,30.,60.,"A/Q s2-s4 corr.", "Z from MUSIC41", 2);

        }
      ///////////////////////////////////////////////////////
 
      //------------------------------------------------------------------------------------------------//

   
      
}

void EventAnlProc::Do_FRS_Histos(EventAnlStore* pOutput){
 
    
     pOutput->pFRS_AoQ = FRS_AoQ;
     pOutput->pFRS_ID_x2 = FRS_ID_x2;
     pOutput->pFRS_ID_y2 = FRS_ID_y2;
     pOutput->pFRS_ID_a2 = FRS_ID_a2;
     pOutput->pFRS_ID_b2 = FRS_ID_b2;
     pOutput->pFRS_ID_x4 = FRS_ID_x4;
     pOutput->pFRS_ID_y4 = FRS_ID_y4;
     pOutput->pFRS_ID_a4 = FRS_ID_a4;
     pOutput->pFRS_ID_b4 = FRS_ID_b4;
     
     pOutput->pFRS_z = FRS_z;
     pOutput->pFRS_z2 = FRS_z2;
     pOutput->pFRS_dEdeg = FRS_dEdeg;
     pOutput->pFRS_dEdegoQ = FRS_dEdegoQ;
     pOutput->pFRS_beta = FRS_beta;
     
     for(int i =0; i<2; ++i){
           pOutput->pFRS_Music_dE[i] = FRS_Music_dE[i]; 
           pOutput->pFRS_Music_dE_corr[i] = FRS_Music_dE_corr[i] ;
          //  cout<<"pOutput->pFRS_Music_dE[i] " <<pOutput->pFRS_Music_dE[i] << " i " << i << endl;
           }
            for(int i =0; i<8; ++i){
          pOutput->pFRS_Music_E1[i] = FRS_Music_E1[i]; 
          pOutput->pFRS_Music_E2[i] = FRS_Music_E2[i];
          pOutput->pFRS_Music_T1[i] = FRS_Music_T1[i]; 
          pOutput->pFRS_Music_T2[i] = FRS_Music_T2[i];
         
            }
            
     for(int l=0; l<12;l++){
     
     pOutput->pFRS_sci_l[l] = FRS_sci_l[l];
     pOutput->pFRS_sci_r[l] = FRS_sci_r[l];
     pOutput->pFRS_sci_e[l] = FRS_sci_e[l];
     pOutput->pFRS_sci_tx[l] = FRS_sci_tx[l];
     pOutput->pFRS_sci_x[l] = FRS_sci_x[l];
      
     }
 for(int i = 0; i < 7; i++)
        {
          pOutput->pFRS_tpc_x[i] = FRS_tpc_x[i];
          pOutput->pFRS_tpc_y[i] =  FRS_tpc_y[i];
  
        }
     
    ///AoQ vs X
    
   if(FRS_AoQ>0 && FRS_ID_x2>-100 && FRS_ID_x2<100)  hID_x2AoQ->Fill(FRS_AoQ, FRS_ID_x2);
   if(FRS_AoQ>0 && FRS_ID_x4>-100 && FRS_ID_x4<100)  hID_x4AoQ->Fill(FRS_AoQ, FRS_ID_x4);

     ///AoQ vs Z
   if(FRS_AoQ>0 && FRS_z>0)  {hID_Z_AoQ->Fill(FRS_AoQ, FRS_z);
   // cout<<  "  FRS_AoQ " << FRS_AoQ <<endl;
     
   }
   if(FRS_AoQ_corr>0 && FRS_z>0)   hID_Z_AoQ_corr->Fill(FRS_AoQ_corr, FRS_z);  //S2-S4 correlated
    
     ///Z1 Z2
   if(FRS_z>0 && FRS_z2>0) hID_Z_Z2->Fill(FRS_z,FRS_z2);
    ///Z1 Z2 no velocity correction
   if(FRS_znocorr>0 && FRS_z2nocorr>0) hID_Z_Z2_nocorr->Fill(FRS_znocorr,FRS_z2nocorr);
     
     ///High charge states
   if(FRS_z>0 && FRS_dEdegoQ!=0)   hdEdegoQ_Z->Fill(FRS_z, FRS_dEdegoQ);
   if(FRS_z>0 && FRS_dEdeg!=0)   hdEdeg_Z->Fill(FRS_z, FRS_dEdeg);
   if(FRS_z>0 && FRS_dEdeg!=0)   hdEdegoQtimesBeta_Z->Fill(FRS_z, FRS_dEdegoQ*FRS_beta);
   if(FRS_beta>0 && FRS_dEdeg!=0)   hdEdegoQ_Beta->Fill(FRS_beta,FRS_dEdegoQ);
     
   if(FRS_AoQ>0 && FRS_dEdeg!=0)  hdEdegoQ_AoQ->Fill(FRS_AoQ,FRS_dEdegoQ);
   if(FRS_AoQ>0 && FRS_dEdeg!=0)  hdEdegoQtimesBeta_AoQ->Fill(FRS_AoQ,FRS_dEdegoQ*FRS_beta);
    
    //  if(TMath::Abs(FRS_z-FRS_z2-0.3)<0.6)
      if(TMath::Abs(FRS_z-FRS_z2)<0.4)
            {
              hID_Z_AoQ_zsame->Fill(FRS_AoQ, FRS_z);
              hID_x4AoQ_zsame->Fill(FRS_AoQ, FRS_ID_x4);
            }

    
        ///GATE: ID vs x2AoQ
     for(int i=0; i<8; i++){
      if(cID_x2AoQ[i]->Test(FRS_AoQ, FRS_ID_x2)==true)
        {
          pOutput->pFRS_x2AoQ_pass[i] = true;
           if(FRS_AoQ>2 && FRS_AoQ<2.4 &&   FRS_ID_x2 > -100 && FRS_ID_x2<100){

          hID_x2AoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x2);


        }
          if(FRS_AoQ>2 && FRS_AoQ<2.4 &&   FRS_ID_x4 > -100 && FRS_ID_x4<100){
          hID_x4AoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x4);}
          if (FRS_z) hID_ZAoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_z);
            }
     }
        ///GATE: ID vs x4AoQ
        for(int i=0; i<8; i++){
      if(cID_x4AoQ[i]->Test(FRS_AoQ, FRS_ID_x4)==true)
        {
                pOutput->pFRS_x4AoQ_pass[i] = true;
                hID_x2AoQ_x4AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x2);
                hID_x4AoQ_x4AoQgate[i]->Fill(FRS_AoQ, FRS_ID_x4);
          if (FRS_z) hID_ZAoQ_x4AoQgate[i]->Fill(FRS_AoQ, FRS_z);
            }
        }
        for(int i=0; i<8; i++){
          ///GATE: Z1 vs Z2
            if(cID_Z_Z2gate[i]->Test(FRS_z, FRS_z2)==true)
            {
                 pOutput->pFRS_Z_Z2_pass[i] = true;


                 if(FRS_AoQ>2 && FRS_AoQ<2.4 &&   FRS_ID_x2 > -100 && FRS_ID_x2<100){
                 
                hID_x2AoQ_Z1Z2gate[i]->Fill(FRS_AoQ, FRS_ID_x2);
                }
                   if(FRS_AoQ>2 && FRS_AoQ<2.4   && FRS_ID_x4 > -100 && FRS_ID_x4<100){
                   
                hID_x4AoQ_Z1Z2gate[i]->Fill(FRS_AoQ, FRS_ID_x4);
                   }
                   if(FRS_AoQ>2 && FRS_AoQ<2.4){
                   hID_ZAoQ_Z1Z2gate[i] ->Fill(FRS_AoQ, FRS_z);
                   }
                 
                hID_SC43Z1_Z1Z2gate[i]->Fill(FRS_sci_e[7], FRS_z);
          if (FRS_z2) hID_ZAoQ_x2AoQgate[i]->Fill(FRS_AoQ, FRS_z2);
            }
        }
     
        ///GATE: AoQ vs Z
        if(cID_Z_AoQ[2]->Test(FRS_AoQ, FRS_z)==true && cID_Z_AoQ[3]->Test(FRS_AoQ, FRS_z)==true) hID_Z_AoQ_overlap->Fill(FRS_AoQ, FRS_z);
        
        for(int g=0; g<8; g++){
       if( cID_Z_AoQ[g]->Test(FRS_AoQ, FRS_z)==true){
          pOutput->pFRS_ZAoQ_pass[g] =true;
	
         hID_Z_AoQgate[g]->Fill(FRS_AoQ, FRS_z);
         hID_x2AoQ_ZAoQgate[g]->Fill(FRS_AoQ, FRS_ID_x2);
         hID_x4AoQ_ZAoQgate[g]->Fill(FRS_AoQ, FRS_ID_x4);
         hID_Z_Z2_ZAoQgate[g]->Fill(FRS_z, FRS_z2);
            }
        }
        
     ///High charge states
      for(int g=0; g<8; g++){
      if(cID_dEdeg_Z1[g]->Test(FRS_z,FRS_dEdeg)==true){ // dEdeg_vs_Z check
        hID_Z1_AoQ_dEdegZgate[g]->Fill(FRS_AoQ, FRS_z);;
        hID_Z1_AoQcorr_dEdegZgate[g]->Fill(FRS_AoQ_corr, FRS_z);
   
    if(fabs(FRS_z2-FRS_z)<0.4){ // in addition Z same
      hID_Z1_AoQ_zsame_dEdegZgate[g]->Fill(FRS_AoQ, FRS_z);
      hID_Z1_AoQcorr_zsame_dEdegZgate[g]->Fill(FRS_AoQ_corr, FRS_z);
            }
        }
      }
    }


/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------    AIDA   ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

void EventAnlProc::Make_Aida_Histos(){
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  implants_strip_xy.resize(conf->DSSDs());
  implants_pos_xy.resize(conf->DSSDs());
  implants_e.resize(conf->DSSDs());
  implants_e_xy.resize(conf->DSSDs());
  implants_time_delta.resize(conf->DSSDs());
  implants_strip_1d.resize(conf->DSSDs());
  implants_per_event.resize(conf->DSSDs());
  decays_strip_xy.resize(conf->DSSDs());
  decays_pos_xy.resize(conf->DSSDs());
  decays_e.resize(conf->DSSDs());
  decays_e_xy.resize(conf->DSSDs());
  decays_time_delta.resize(conf->DSSDs());
  decays_strip_1d.resize(conf->DSSDs());
  decays_per_event.resize(conf->DSSDs());
  implants_channels.resize(conf->DSSDs());
  decays_channels.resize(conf->DSSDs());
  implants_x_ex.resize(conf->DSSDs());
  implants_y_ey.resize(conf->DSSDs());
#ifdef AIDA_PULSER_ALIGN
  aida_pulser_time = MakeTH2('I', "AIDA/Pulser_Time", "AIDA Pulser Time Comparison", 768, 0, 768, 2000, -4000, 4000);
#endif
  for (int i = 0; i < conf->DSSDs(); ++i)
  {
    implants_strip_xy[i] = MakeTH2('I', Form("AIDA/Implants/DSSD%d_implants_strip_XY", i+1), Form("DSSD %d implant hit pattern", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
    //implants_pos_xy[i] = MakeTH2('D', Form("AIDA/Implants/DSSD%d_implants_pos_XY", i+1), Form("DSSD %d implant position", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    implants_e[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_energy", i+1), Form("DSSD %d implant energy", i+1), 1000, 0, 10000, "Implant Energy/MeV");
    //implants_e_xy[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_energy_XY", i+1), Form("DSSD %d implant front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 10000, "X Energy", "Y Energy");
    implants_time_delta[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_time_delta", i+1), Form("DSSD %d implant front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    implants_strip_1d[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_strip_1d", i+1), Form("DSSD %d implant 1D hit pattern", i+1), 256, 0, 256, "Strip number");
    implants_per_event[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_per_event", i+1), Form("DSSD %d implants per event", i+1), 100, 0, 100, "Number of implants");
    implants_x_ex[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_x_ex", i+1), Form("DSSD %d Ex vs X position", i+1), 128, 0, 128, 1000, 0, 10000, "X Strip", "X Energy");
    implants_y_ey[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_y_ey", i+1), Form("DSSD %d Ey vs Y position", i+1), 128, 0, 128, 1000, 0, 10000, "Y Strip", "Y Energy");

    decays_strip_xy[i] = MakeTH2('I', Form("AIDA/Decays/DSSD%d_decays_strip_XY", i+1), Form("DSSD %d decay hit pattern", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
   // decays_pos_xy[i] = MakeTH2('D', Form("AIDA/Decays/DSSD%d_decays_pos_XY", i+1), Form("DSSD %d decay position", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    decays_e[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_energy", i+1), Form("DSSD %d decay energy", i+1), 1000, 0, 20000, "Decay Energy/keV");
    //decays_e_xy[i] = MakeTH2('F', Form("AIDA/Decays/DSSD%d_decays_energy_XY", i+1), Form("DSSD %d decay front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 20000, "X Energy", "Y Energy");
    decays_time_delta[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_time_delta", i+1), Form("DSSD %d decay front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    decays_strip_1d[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_strip_1d", i+1), Form("DSSD %d decay 1D hit pattern", i+1), 256, 0, 256, "Strip number");
    decays_per_event[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_per_event", i+1), Form("DSSD %d decays per event", i+1), 100, 0, 100, "Number of decays");

    implants_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_implants_channels", i+1), Form("DSSD %d number of implant channels", i+1), 769, 0, 769);
    decays_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_decays_channels", i+1), Form("DSSD %d number of decay channels", i+1), 769, 0, 769);
  }
}

///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///

void EventAnlProc::ProcessAida(EventUnpackStore* pInputMain, EventAnlStore* pOutput) {
 // int Aida_hits =0;
//       double bPlasQDCGainMatch_AIDA[32] ={0};
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();


  int ev = 0;
  for(AidaUnpackData& pInputD : pInputMain->Aida)
  {
    ev++;
    AidaUnpackData* pInput = &pInputD;
    pOutput->pAIDA_WR = pInputMain->fAIDA_WR;
    AidaAnlData aida;
    if (pInput->ImplantEvents.size() > 1)
    {

      //         cout << " pInput->ImplantEvents.size() " << pInput->ImplantEvents.size() <<  endl;

      implantEvents++;


      // Cluster events on adjecent strips into one
      std::vector<AidaCluster> clusters = EventsToClusters(pInput->ImplantEvents);
      //
      //     // Match front-back clusters which define physical hits on the detector
      std::vector<std::pair<AidaCluster, AidaCluster>> hits;
      //
      //
      // Track the number of implants in every DSSD
      std::vector<int> counts(conf->DSSDs(), 0);
      int max_dssd = 0;

      for (auto& i : clusters)
      {
        if (i.DSSD == -1) continue;

        counts[i.DSSD - 1]++;
        if (i.DSSD > max_dssd) max_dssd = i.DSSD;

        if(i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        for (auto& j : clusters)
        {
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          /// Gates (set in TAidaConfiguration)
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyH() && i.IsGoodTime(j, conf->FrontBackWindow()))
          {
            hits.push_back({i, j});
          }
        }
      }

      int channels[768] = {0};
      for (auto& i : pInput->ImplantEvents)
      {
        channels[i.Module * 64 + i.Channel]++;
      }
      int channelM = 0;
      for (int i = 0; i < 768; ++i){
        if (channels[i]) ++channelM;
      }
      implants_channels[0]->Fill(channelM);

      //     // Generate stored data for hits and plot the histograms

      for (auto& i : hits)
      {
        AidaHit hit = ClusterPairToHit(i);
        hit.Event = ev;

        hit.Stopped = (hit.DSSD == max_dssd);

        // Check that every DSSD before has at least one implant event
        for(int j = hit.DSSD - 1; j > 0; j--)
        {
            if (counts[j - 1] == 0) hit.Stopped = false;
        }
        //
        //
        if (hit.Stopped) stoppedEvents++;
//         pInput->Implants.push_back(hit);
// 
//         pOutput->pAida.Implants.push_back(hit);
        
        aida.Implants.push_back(hit);

        implants_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
       // implants_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        implants_e[hit.DSSD - 1]->Fill(hit.Energy);
    //   cout<<"ANALYSIS AIDA " <<pOutput->pEvent_Number<< " Energy " <<  hit.Energy << endl;

      //  implants_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //implants_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        implants_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        implants_x_ex[hit.DSSD - 1]->Fill(hit.StripX, hit.EnergyFront);
        implants_y_ey[hit.DSSD - 1]->Fill(hit.StripY, hit.EnergyBack);

        int channel = i.first.Strip;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + 128;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
      }
      //
      if (hits.size() > 0)
      {
        implants_per_event[0]->Fill(hits.size());
        goodImplantEvents++;
      }
    }
    else if (pInput->DecayEvents.size() > 1)
    {
      decayEvents++;

      int channels[768] = {0};
#ifdef AIDA_PULSER_ALIGN
      int64_t wr_base = 0;
#endif
      for (auto& i : pInput->DecayEvents)
      {
        channels[i.Module * 64 + i.Channel]++;
#ifdef AIDA_PULSER_ALIGN
        if(i.Module == 0 && i.Channel == 0) wr_base = i.Time;
#endif
      }
      int channelM = 0;
      for (int i = 0; i < 768; ++i)
      if (channels[i]) ++channelM;
      decays_channels[0]->Fill(channelM);

      if (channelM > 400)
      {
        decayEvents--;
        pulserEvents++;
#ifdef AIDA_PULSER_ALIGN
        if(channelM > 750)
        {
          std::cout << "Identified a pulser event!" << std::endl;
          std::vector<int> offsets;
          offsets.resize(conf->FEEs());
          for (int i = 0; i < offsets.size(); i++) offsets[i] = 0;

          for (auto& i : pInput->DecayEvents)
          {
            if (i.Energy < 1000) continue;
            std::cout << i.Module << ", " << i.Channel << ", " << i.Energy << ", " << std::hex << i.Time << "," << i.FastTime << std::dec << std::endl;
            int offset = (i.Time - wr_base) % 2000;
            if (offsets[i.Module] == 0)
            {
              offsets[i.Module] = offset;
            }
            else if (offset > offsets[i.Module])
            {
              // confirm the offset is 2us out
              if (abs(offset - offsets[i.Module]) % 2000 != 0)
                std::cout << "LOGICAL MISTAKE IN ALIGNMENT" << std::endl;
            }
            else if (offset < offsets[i.Module])
            {
              if (abs(offset - offsets[i.Module]) % 2000 != 0)
                std::cout << "LOGICAL MISTAKE IN ALIGNMENT" << std::endl;
              offsets[i.Module] = offset;
            }

            aida_pulser_time->Fill(i.Module * 64 + i.Channel, offsets[i.Module]);
          }

          std::cout << std::endl;

          std::cout << "Put this into AIDA_Times.txt" << std::endl;
          for (int i  = 0; i < offsets.size(); i++)
          {
            std::cout << i << " " << offsets[i] << std::endl;
          }

          throw TGo4UserException(3,"");
        }
#endif
        return;
      }

      // Clean up huge event buffers - for now we just destroy them
      if (pInput->DecayEvents.size() > 400)
      {
        decayEvents--;
        nonsenseEvents++;
//        pInput->SetValid(kFALSE); ///NEEDED!?
        return;
      }

      std::vector<AidaCluster> clusters = EventsToClusters(pInput->DecayEvents);
      std::vector<std::pair<AidaCluster, AidaCluster>> hits;

      for (auto& i : clusters)
      {

        if(i.DSSD == -1 || i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        //if(i.Energy < 100) continue;
        for (auto& j : clusters)
        {
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          //if(j.Energy < 100) continue;
          // Gates
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyL() && i.IsGoodTime(j, conf->FrontBackWindow()))
          hits.push_back({i, j});
        }
      }

      for (auto& i : hits)
      {
        AidaHit hit = ClusterPairToHit(i);
        hit.Event = ev;

      

         if (hit.DSSD != 3 || (hit.EnergyBack > 150 && hit.EnergyFront > 150))
          aida.Decays.push_back(hit);
        decays_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
        //decays_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        decays_e[hit.DSSD - 1]->Fill(hit.Energy);
       // decays_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //decays_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        decays_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        int channel = i.first.Strip;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + 128;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
      }

      if (clusters.size() > 0)
      {
        decays_per_event[0]->Fill(clusters.size());
      }
      //   cout <<" decayEvents " << decayEvents << endl;
    }
    else
    {
      nonsenseEvents++;
      //pInput->SetValid(kFALSE);
    }
        pOutput->pAida.push_back(aida);
  }
}



std::vector<AidaCluster> EventAnlProc::EventsToClusters(std::vector<AidaEvent> const& events)
{
  std::vector<AidaCluster> clusters;
  for (auto& i : events)
  {
    // Don't cluster invalid events
    if (i.DSSD == -1) continue;

    bool added = false;
    AidaCluster* cluster;

    // Try to add the event to an existing cluster
    // Weird loop because we can remove clusters too
    auto it = std::begin(clusters);
    while (it != std::end(clusters))
    {
      auto& j = *it;
      if(j.IsAdjacent(i) && j.IsGoodTime(i))
      {
        // Add the event to the cluster
        if (!added)
        {
          j.AddEvent(i);
          cluster = &j;
          added = true;
        }
        // If we match again the two clusters are merged into one
        // The old cluster is then removed
        else
        {
          cluster->AddCluster(j);
          it = clusters.erase(it);
          continue;
        }
      }
      ++it;
    }

    // Otherwise make a new cluster for the event
    if (!added)
    {
      AidaCluster c_test;
      c_test.AddEvent(i);
      clusters.push_back(c_test);
    }
  }
  return clusters;
}

AidaHit EventAnlProc::ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const& i)
{
  AidaHit hit;
  hit.DSSD = i.first.DSSD;

  hit.StripX = i.first.Strip;
  hit.StripY = i.second.Strip;
  hit.PosX = 75.6 * i.first.Strip / 128. - 37.75;
  hit.PosY = 75.6 * i.second.Strip / 128. - 37.75;

  hit.StripXMin = i.first.StripMin;
  hit.StripXMax = i.first.StripMax;
  hit.StripYMin = i.second.StripMin;
  hit.StripYMax = i.second.StripMax;
  hit.ClusterSizeX  = i.first.N;
  hit.ClusterSizeY = i.second.N;

  hit.Energy = (i.first.Energy + i.second.Energy) / 2;
  hit.EnergyFront = i.first.Energy;
  hit.EnergyBack = i.second.Energy;

  hit.Time = std::min(i.first.Time, i.second.Time);
  hit.TimeFront = i.first.Time;
  hit.TimeBack = i.second.Time;
  hit.FastTime = std::min(i.first.FastTime, i.second.FastTime);
  hit.FastTimeFront = i.first.FastTime;
  hit.FastTimeBack = i.second.FastTime;

  return hit;
}
                                                    ///End of Aida ///

 /**----------------------------------------------------------------------------------------------**/   
     /**--------------------------------------  bPlastic  TAMEX ----------------------------------------------   **/
     //**----------------------------------------------------------------------------------------------**/  
        
    void EventAnlProc::Make_Plastic_Tamex_Histos(){ 
        
         for (int i =1; i<3; i++)   
            { 
                
                 hbPlas_ToT_Sum[i] = MakeTH1('D', Form("bPlastic/ToT_Sum_Det.%2d",i), Form("bPlastic Sum Gainmatched ToT Det. %2d",i-FATIMA_TAMEX_MODULES), 2500, 0, 1500000);
  
                  hbPlas_hit_pattern_det[i]= MakeTH1('D', Form("bPlastic/Stats/HitPattern_Det.%2d",i), Form("bPlastic Hit pattern Det. %2d",i), 17, 0, 17);
             
          for(int j=0; j<16; j++){  
              
             hbPlas_Lead_T[i][j] = MakeTH1('D', Form("bPlastic/LeadTime/Lead T Plas Det. %2d Ch.%2d",  i,j), Form("Lead - Time Det %2d Ch. %2d", i,j),2500, 0, 2000);
             hbPlas_Trail_T[i][j] = MakeTH1('D', Form("bPlastic/TrailTime/Trail T Plas Det. %2d Ch.%2d",  i,j), Form("Trail - Time Det %2d Ch. %2d", i,j),2500, 0, 2000);
               
               
            hbPlas_ToT_det[i][j] = MakeTH1('D', Form("bPlastic/ToT/ToT Plas Det. %2d Ch. %2d",  i,j), Form("ToT Det. %2d Ch. %2d", i,j),10000, 0., 100000.);   
            
           hbPlas_Energy_Calib[i][j] = MakeTH1('D', Form("bPlastic/Energy_Calib/Energy Calib Plas Det. %2d Ch. %2d",  i,j), Form("Energy Calib Det. %2d Ch. %2d", i,j),10000, 0., 100000.);    
            
          //if(i==bPLAS_TAMEX_ID){    
            
//         hbPlas_lead_lead[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead/Lead-Lead Plas Det. %2d Ch.%2d",  i,j), Form("Lead - Lead Det %2d Ch. %2d", i,j),2500, -50000., 50000);  
        
        hbPlas_lead_lead_ref_det1[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Ref/Lead-Lead Plas Det. %2d RefCh. %2d", i,j), Form("Lead Ref Ch.0 - Lead Det.%2d Ch. %2d", i,j),2500, -50000., 50000.);
        hbPlas_lead_lead_ref_det2[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Ref/Lead-Lead Plas Det. %2d RefCh. %2d", i,j), Form("Lead Ref Ch.0 - Lead Det.%2d Ch. %2d", i,j),2500, -50000., 50000.);
            
        hbPlas_lead_lead_gated[i][j] = MakeTH1('D', Form("bPlastic/Lead-Lead_Egated/Lead-Lead Egated Plas Det. %2d Ch. %2d",  i,j), Form("Lead - Lead Energy gated Det. %2d Ch.  %2d", i,j),2500, -50000., 50000.); 
            
       // hbPlas_SC41L_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41-Lead_Plas/SC41_Lead Plas Det. %2d Ch.%02d", i,j), Form("SC41 Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
         hbPlas_SC41L_Anal_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41L_Anal-Lead_bPlas/SC41L_Anal_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41L Analogue Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
            
         hbPlas_SC41R_Anal_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41R_Anal-Lead_bPlas/SC41R_Anal_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41R Analogue Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);    
            
         hbPlas_SC41L_Digi_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41L_Digi-Lead_bPlas/SC41L_Digi_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41L Digital Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.); 
            
         hbPlas_SC41R_Digi_lead[i][j] = MakeTH1('D', Form("bPlastic/SC41R_Digi-Lead_bPlas/SC41R_Digi_Lead bPlas Det. %2d Ch.%02d", i,j), Form("SC41R Digital Lead - bPlas Lead Det. %2d Lead Ch. %2d ", i,j), 4002, -100000., 100000.);             
        
          }
            }
        
                                 
        hSC41_Analogue_Tamex = MakeTH1('D',"bPlastic/SC41/Analogue L-R","SC41 Analogue L - R",4002, -100000., 100000.); 
        hSC41_Digital_Tamex = MakeTH1('D',"bPlastic/SC41/Digital L-R","SC41 Analogue L - R",4002, -100000., 100000.);   
           
            
        hbPlas_Multiplicity_Det1 = MakeTH1('D',"bPlastic/Stats/Multiplicity_Det1","bPlastic Multiplicity Det 1",32,0,32);   
        hbPlas_Multiplicity_Det2 = MakeTH1('D',"bPlastic/Stats/Multiplicity_Det2","bPlastic Multiplicity Det 2",32,0,32);   
       
      
        
    }   
    /////////////////////////////////////////////////// 
    void EventAnlProc::Do_Plastic_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){   
           
        bool fired_det1=false, fired_det2=false;    
        int bPlas_tot_hits = 0;
       
         for(int a=1; a<3; a++){
                 for (int b = 0; b < 16; b++){  
                     for(int k=0; k<10; k++){
                        lead_bplas[a][b][k]=0; 
                        ToT_bplas[a][b][k] = 0;   
                     }
                 }
         }
         for(int i=0; i<10; i++){
         bPlas_RefCh0_Det1[i] =0;
         bPlas_RefCh0_Det2[i] =0;
         }
        for (int i = 0; i < 48; i++)    
            {   
                for(int j=0; j<10;j++){    

                 if(i < 16) lead_lead_bplas_Ref1[i][j]=0;  
                 lead_lead_fat_Ref0[i][j]=0;    
  
                 SC41L_ANA_lead_bPlas[i][j] = 0;    
                 SC41R_ANA_lead_bPlas[i][j] = 0;    
                 SC41L_DIG_lead_bPlas[i][j] = 0;    
                 SC41R_DIG_lead_bPlas[i][j] = 0;    
                
                }   
            }
          
            
        
     ///**---------------------------------------------LEAD -------------------------------------------------**/        
           ///Loop on channels First    
                           
              for(int i=1; i<3; i++){ ///Detector number
                 for (int j = 0; j < 16; j++){  ///Channel number 
                     
                for(int k=0; k< pInput->fbPlas_PMT_Lead_N[i][j]; k++){ 
                    //Fat_RefCh[j] = pInput->fFat_Lead_PMT[1][j]; 
                    bPlas_RefCh0_Det1[k] = pInput->fbPlas_Lead_PMT[1][0][k];
                    bPlas_RefCh0_Det2[k] = pInput->fbPlas_Lead_PMT[2][0][k];
            
                        }      
                    }
              }
            ////////////////////////////    
              ///Loop over channels 
              for(int a=1; a<3; a++){ ///Detector number
                 for (int b = 0; b < 16; b++){  ///Channel number 
                     
 ///**---------------------------------------------Plastic Lead Time ----------------------------------**/    
                                                            
                pOutput->pbPlas_PMT_Lead_N[a][b] = pInput->fbPlas_PMT_Lead_N[a][b]; 
              //cout<<"pInput->fbPlas_PMT_Lead_N[a][b] " << pInput->fbPlas_PMT_Lead_N[a][b] <<endl;
               for(int j=0; j< pInput->fbPlas_PMT_Lead_N[a][b]; j++){ ///Hits 
                
                    
                    lead_bplas[a][b][j] = pInput->fbPlas_Lead_PMT[a][b][j];  
                    hbPlas_Lead_T[a][b]->Fill(lead_bplas[a][b][j]);
                    hits_bplas_lead++;  
                    pOutput->pbPlas_LeadT[a][b][j] = lead_bplas[a][b][j];    
                    pOutput->pbPlas_LeadHits = hits_bplas_lead; 
                    pOutput->pbPlas_LeadT_Avg = lead_bplas[a][b][j]/hits_bplas_lead;
                // cout<<"event nunber "<<pOutput->pEvent_Number<<" fbPlas_Lead_PMT[a][b][j] "<<pInput->fbPlas_Lead_PMT[a][b][j] << " a " << a << " b " << b << " j " << j <<endl;
        
        
    ///**---------------------------------------------Plastic Lead Ref dT ----------------------------------**/          
                          
//                     if(i>15 && pInput->fbPlas_Lead_PMT[16][j]>0 && pInput->fbPlas_Lead_PMT[a][b][j]>0) {   
                
        if(bPlas_RefCh0_Det1[j]>0 && lead_bplas[1][b][j]>0){
            lead_lead_bplas_Ref1[b][j] = (bPlas_RefCh0_Det1[j] -  lead_bplas[1][b][j])*CYCLE_TIME; 
        }
        if(bPlas_RefCh0_Det2[j]>0 && lead_bplas[2][b][j]>0){
                      lead_lead_bplas_Ref2[b][j] = (bPlas_RefCh0_Det2[j] -  lead_bplas[2][b][j])*CYCLE_TIME;
        }
                      
          //    if(lead_lead_bplas_Ref1[b][j]!=0) hbPlas_lead_lead_ref_det1[1][j] ->Fill(lead_lead_bplas_Ref1[b][j]);
//              if(lead_lead_bplas_Ref2[b][j]!=0) hbPlas_lead_lead_ref_det2[2][j] ->Fill(lead_lead_bplas_Ref2[b][j]);
           
                }
              
                                
      ///**---------------------------------------------Plastic Trail ----------------------------------**/  
                                                            
                pOutput->pbPlas_PMT_Trail_N[a][b] = pInput->fbPlas_PMT_Trail_N[a][b]; 
                
               for(int j=0; j< pInput->fbPlas_PMT_Trail_N[a][b]; j++){ ///Hits 

                    
                    trail_bplas[a][b][j] = pInput->fbPlas_Trail_PMT[a][b][j];  
                    hbPlas_Trail_T[a][b]->Fill(trail_bplas[a][b][j]);
                    hits_bplas_trail++;  
                    pOutput->pbPlas_TrailT[a][b][j] = trail_bplas[a][b][j];  
                     }
   ///**---------------------------------------------Plastic ToT ----------------------------------**/  
              for(int j=0; j< pInput->fbPlas_PMT_Lead_N[a][b]; j++){ 
                  
                  if(pInput->fbPlas_Trail_PMT[a][b][j] >0 && pInput->fbPlas_Lead_PMT[a][b][j]>0){ 
                
        ToT_bplas[a][b][j] = (pInput->fbPlas_Trail_PMT[a][b][j] - pInput->fbPlas_Lead_PMT[a][b][j]);   
                    
                ///Correction for overflows 
                if(ABS(ToT_bplas[a][b][j]) >(double)(COARSE_CT_RANGE>>1)) {   
                    
                       ToT_bplas[a][b][j] = CYCLE_TIME*(ToT_bplas[a][b][j] + COARSE_CT_RANGE);    
                      } 
                 else{  
                           ToT_bplas[a][b][j]= CYCLE_TIME*ToT_bplas[a][b][j];                         
                       }    
                       ///Gain matching  
               // pOutput-> pbPlas_ToTCalib[a][b][j] = fCal->Abplas_TAMEX_ZAoQ[i]* ToT_bplas[a][b][j] + fCal->Bbplas_TAMEX_ZAoQ[i];
               pOutput-> pbPlas_ToTCalib[a][b][j] =ToT_bplas[a][b][j];
              
                       if(ToT_bplas[a][b][j]>0) {
                        hbPlas_ToT_det[a][b] ->Fill(ToT_bplas[a][b][j]);   
                        hbPlas_ToT_Sum[a]->Fill(ToT_bplas[a][b][j]);   
                        hbPlas_hit_pattern_det[a]->Fill(b);
                     //   cout<<"Event " << pOutput->pEvent_Number  <<" a " << a << " b " << b << endl;
                    
                          bPlas_tot_hits++; 
          
                         if(a==1) hbPlas_Multiplicity_Det1->Fill(bPlas_tot_hits);
                         if(a==2) hbPlas_Multiplicity_Det2->Fill(bPlas_tot_hits);
                            }
                        }         
                    }
                 }
              }
            }
    
    /**-----------------------------------------------------------------------------------------------**/
    /**--------------------------------------  FATIMA TAMEX ------------------------------------------**/
    /**-----------------------------------------------------------------------------------------------**/
 void EventAnlProc::Make_Fatima_Tamex_Histos(){
     
     for(int i=0; i<48; i++){ 
         hFat_Lead_T[i] =  MakeTH1('D', Form("FATIMA_TAMEX/LeadT/Lead Time Ch.%2d",i), Form("Lead time. %2d",i), 5000,0,5000);
         hFat_Trail_T[i] =  MakeTH1('D', Form("FATIMA_TAMEX/TrailT/Trail Time Ch.%2d",i), Form("Trail time. %2d",i), 5000,0,5000);
         hFat_lead_lead_ref[i] =   MakeTH1('D', Form("FATIMA_TAMEX/LeadRef-Lead/Lead-Lead Time Ref Ch1- Ch.%2d",i), Form("RefLead-Lead time. %2d",i), 1000,-10000,10000);
         hFat_ToT_det[i]=   MakeTH1('D', Form("FATIMA_TAMEX/ToT/ToT Ch.%2d",i), Form("ToT %2d",i), 15000,0,150000);        
     }
          hFat_ToT_Sum= MakeTH1('D', "FATIMA_TAMEX/ToTSum", "ToT LaBr (all detectors)",15000,0,150000);
          hFat_tamex_hit_pattern =  MakeTH1('D', "FATIMA_TAMEX/Fatima_Hitpattern", "Fatima Hit pattern",48,0,48);  
          hFat_tamex_multiplicity =  MakeTH1('D', "FATIMA_TAMEX/Fatima_Multiplicity", "Fatima Multiplicity",48,0,48);  
 }
 
 
//-----------------------------------------------------------------------------------------------//
void EventAnlProc::Do_Fatima_Tamex_Histos(EventUnpackStore* pInput, EventAnlStore* pOutput){   
          
        //bool fired_det1=false, fired_det2=false;    
        int Fat_tot_hits = 0;
       
        // for(int a=1; a<3; a++){
                 for (int i = 0; i < 48; i++){  
                     for(int j=0; j<10; j++){
                        lead_fat[i][j]=0; 
                        ToT_fat[i][j] = 0; 
                        lead_lead_fat_Ref1[i][j]=0;
//                         SC41L_ANA_lead_fat[i][j] = 0;  
//                         SC41R_ANA_lead_fat[i][j] = 0;  
//                         SC41L_DIG_lead_fat[i][j] = 0;  
//                         SC41R_DIG_lead_fat[i][j] = 0;  
                     }
                 }
         
         for(int i=0; i<10; i++){
         FatTam_RefCh0[i] =0;
         SC41L_ANA_lead_fat[i] =0;
         SC41R_ANA_lead_fat[i] =0;
         bPlasDet1_coin_lead_Fat[i] =0;
         bPlasDet2_coin_lead_Fat[i] =0;
         
         }
                
     ///**---------------------------------------------LEAD -------------------------------------------------**/        
           ///Loop on channels First    
                           
              for(int i=0; i<48; i++){ ///Channel number
                 for (int j = 0; j < pInput->fFat_PMT_Lead_N[i]; j++){  ///Hit 
                     
                    
                    FatTam_RefCh0[j] = pInput->fFat_Lead_PMT[1][j];
                    SC41L_ANA_lead_fat[j]=pInput->fFat_Lead_PMT[4][j];
                    SC41R_ANA_lead_fat[j]=pInput->fFat_Lead_PMT[5][j];
                    bPlasDet1_coin_lead_Fat[j]=pInput->fFat_Lead_PMT[6][j];
                    bPlasDet2_coin_lead_Fat[j]=pInput->fFat_Lead_PMT[7][j];
                 
                   
                        }      
              }
            ////////////////////////////    
              ///Loop over channels 
              for(int i=0; i<FATIMA_TAMEX_CHANNELS; i++){ ///Channel number
                
                     
 ///**---------------------------------------------Fatima Lead Time ----------------------------------**/    
   
                for (int j = 0; j < pInput->fFat_PMT_Lead_N[i]; j++){  ///Hit 

                    lead_fat[i][j] = pInput->fFat_Lead_PMT[i][j];  
                    hFat_Lead_T[i]->Fill(lead_fat[i][j]);
                    hits_fat_lead++;  
                    pOutput->pFat_LeadT[i][j] = lead_fat[i][j];
                    //cout<<"pOutput->pFat_LeadT[i][j] " <<pOutput->pFat_LeadT[i][j] << " i " << i << " j " << j << endl;
                    pOutput->pFat_LeadHits = hits_fat_lead; 
                    pOutput->pFat_Tamex_chan = i;
                  // cout<<"lead_fat[i][j] " << lead_fat[i][j]<< " i " << i << " j " << j <<endl;
        
        
    ///**---------------------------------------------Fatima Lead Ref dT ----------------------------------**/          
                          
//                     if(i>15 && pInput->fFat_Lead_PMT[16][j]>0 && pInput->fFat_Lead_PMT[i][j]>0) {   
                
        if(FatTam_RefCh0[j]>0 && lead_fat[i][j]>0){
            lead_lead_fat_Ref1[i][j] = (FatTam_RefCh0[j] -  lead_fat[i][j])*CYCLE_TIME; 
//             cout<<"pOutput->Event_Number " << pOutput->pEvent_Number << " FatTam_RefCh0[j] " << FatTam_RefCh0[j] << " lead_fat[i][j] " << lead_fat[i][j] <<" lead_lead_fat_Ref1[i][j] " <<lead_lead_fat_Ref1[i][j] <<  " i " << i << " j  " << j << endl;
        }
                              
              if(lead_lead_fat_Ref1[i][j]!=0) hFat_lead_lead_ref[i] ->Fill(lead_lead_fat_Ref1[i][j]);
          
           
                }
                         
      ///**---------------------------------------------Fatima Trail ----------------------------------**/  
                                                            
               // pOutput->pFat_PMT_Trail_N[a][b] = pInput->fFat_PMT_Trail_N[a][b]; 
                
               for(int j=0; j< pInput->fFat_PMT_Trail_N[i]; j++){ ///Hits 

                    
                    trail_fat[i][j] = pInput->fFat_Trail_PMT[i][j];  
                    hFat_Trail_T[i]->Fill(trail_fat[i][j]);
//                    hits_fat_trail++;  
                    pOutput->pFat_TrailT[i][j] = trail_fat[i][j];  
                     }
   ///**---------------------------------------------Fatima  ToT ----------------------------------**/  
              for(int j=0; j< pInput->fFat_PMT_Lead_N[i]; j++){ 
                 // cout<< " pInput->fFat_Trail_PMT[i][j] " << pInput->fFat_Trail_PMT[i][j]<< " pInput->fFat_Lead_PMT[i][j] " <<pInput->fFat_Lead_PMT[i][j] <<" i " <<i << " j " << j << endl;
                  if(pInput->fFat_Trail_PMT[i][j] >0 && pInput->fFat_Lead_PMT[i][j]>0){ 
                
        ToT_fat[i][j] = (pInput->fFat_Trail_PMT[i][j] - pInput->fFat_Lead_PMT[i][j]);   
       // cout<<"pInput->fFat_Trail_PMT[i][j] " <<pInput->fFat_Trail_PMT[i][j]<< " pInput->fFat_Lead_PMT[i][j] " <<pInput->fFat_Lead_PMT[i][j]<< endl;
              
                ///Correction for overflows 
                if(ABS(ToT_fat[i][j]) >(double)(COARSE_CT_RANGE>>1)) {   
                    
                       ToT_fat[i][j] = CYCLE_TIME*(ToT_fat[i][j] + COARSE_CT_RANGE);    
                      } 
                 else{  
                           ToT_fat[i][j]= CYCLE_TIME*ToT_fat[i][j];                         
                       }    
                       ///Gain matching  
               // pOutput-> pFat_ToTCalib[i][j] = fCal->Afat_TAMEX_ZAoQ[i]* ToT_fat[i][j] + fCal->Bfat_TAMEX_ZAoQ[i];
               pOutput-> pFat_ToTCalib[i][j] =ToT_fat[i][j];
           //    cout<<"pOutput-> pFat_ToTCalib[i][j] " <<pOutput-> pFat_ToTCalib[i][j]<< endl;
                       if(ToT_fat[i][j]>0) {
                        hFat_ToT_det[i] ->Fill(ToT_fat[i][j]);   
                        hFat_ToT_Sum->Fill(ToT_fat[i][j]);   
                        hFat_tamex_hit_pattern->Fill(i);  
                    
                          Fat_tot_hits++; 
                          hFat_tamex_multiplicity->Fill(Fat_tot_hits);
                            
                            }
                        }         
                    }
                 }
              }
            
  /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------  FATIMA VME ----------------------------------------------**/

void EventAnlProc::Make_Fatima_Histos(){


  for (int i=0; i<50; i++){
    hFAT_Energy[i] =  MakeTH1('D', Form("FATIMA_VME/Energy/EnergyCalib/LaBr_ECalib_Ch.%2d",i), Form("QDC Calib Ch. %2d",i), 4000,0,4000);
   // hFAT_QDCdt[i]   = MakeTH1('D', Form("FATIMA_VME/Timing/QDCdt/QDCdt%2d",i), Form("QDCdT Ch.%2d",i), 3201,-40,40);
    //hFAT_TDCCalib1[i] =  MakeTH1('D', Form("FATIMA/Timing/TDCCalib/LaBr_Tcalib%2d",i), Form("TDC channel Calib %2d",i), 1E5,0,2E5);
   
    hFAT_TDCdt_refCha[i] = MakeTH1('D', Form("FATIMA_VME/Timing/TDCdT_TDC0-TDC/TDCdT_Cha_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel All Multip LaBr%02d - LaBr%02d",0 , i),250,-2E4,2E4);
   
  }

    hFAT_EnergySum = MakeTH1('D', "FATIMA_VME/Energy/Fat_VME_EnergySum", "LaBr Energy (all detectors)",40000,0,40000);
    hFAT_hits_QDC       = MakeTH1('D', "FATIMA_VME/Stats/QDC_FAThits", "FATIMA QDC hit pattern",50,0,50);
    hFAT_hits_TDC       = MakeTH1('D', "FATIMA_VME/Stats/TDC_FAThits", "FATIMA TDC hit pattern",50,0,50);
    hFAT_Multipl       = MakeTH1('D', "FATIMA_VME/Stats/Fatima_VME_Multiplicity", "FATIMA TDC Multiplicity",50,0,50);

    hFAT_SC41_check      = MakeTH1('D', "FATIMA_VME/Timing/dt_SC41_L_R", "SC41 diff L vs R", 4E4,-2E6,2E6);

}
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///
void EventAnlProc::Do_Fatima_Histos(EventAnlStore* pOutput){
    Long64_t Fat_Ch_dT[50];
   
    pOutput->pFAT_WR = Fat_WR;
    pOutput->pFatmult =  Fatmult;
    
    //cout << " anl proc fat mult " << Fatmult << endl; 
  
    
if(Fatmult > 0){
            hFAT_Multipl->Fill(Fatmult);
    //QDC and TDC outputs
    for (int i = 0; i<Fatmult; i++){
      pOutput->pFat_QDC_ID[i] = Fat_QDC_ID[i];
      pOutput->pFat_QDC_E[i] = Fat_QDC_E[i];
      pOutput->pFat_QDC_T_coarse[i] = Fat_QDC_T_coarse[i];
      pOutput->pFat_QDC_T_fine[i] = Fat_QDC_T_fine[i];
      pOutput->pFat_QDC_E_Raw[i] = Fat_QDC_E_Raw[i];

      pOutput->pFat_TDC_ID[i] = Fat_TDC_ID[i];
      pOutput->pFat_TDC_T[i] = Fat_TDC_T[i];
      pOutput->pFat_TDC_T_Raw[i] = Fat_TDC_T_Raw[i];
    
      ///Fatima Energy 
      hFAT_Energy[Fat_QDC_ID[i]]->Fill(Fat_QDC_E[i]);    
      hFAT_EnergySum->Fill(Fat_QDC_E[i]);
      
          
         ///Hit patterns
          hFAT_hits_QDC->Fill(Fat_QDC_ID[i]);
          hFAT_hits_TDC->Fill(Fat_TDC_ID[i]);
             /// Fatima Time PM 1 - PM Ch.x Energy gated
              Fat_Ch_dT[Fat_TDC_ID[i]] =  (Fat_Cha_Ref_TDC - Fat_TDC_T[i]);

                    if( Fat_TDC_ID[i]!=0){
                        hFAT_TDCdt_refCha[i]->Fill(Fat_Ch_dT[i]);
                    }   
                }
      pOutput->pSC40mult =  SC40mult;
      pOutput->pSC41mult =  SC41mult; 
 
        for(int i = 0; i < SC40mult; i++){
            
            pOutput->pSC40[i] = SC40[i];    
            
        }//End sc40 for fill
        
        for(int i = 0; i < SC41mult; i++){
            
            pOutput->pSC41[i] = SC41[i];    
            
        }//End sc41 for fill        
    }     
    
        pOutput->pFat_bPlasmult= Fat_bPlasmult;
        for(int i = 0; i < Fat_bPlasmult; i++){
            
            pOutput->pFat_bPlas[i] = Fat_bPlas[i];
            
        }//End bPlas Fatima coinc for fill
        
       pOutput->pstdcmult =  stdcmult;
       pOutput->pstdcmult =  stdcmult; 


		//Outputting singles TDC and QDC data
      for(int i = 0; i < stdcmult; i++){
            
         pOutput->pFat_TDC_Singles_t[i] = Fat_TDC_Singles_t[i];    
         pOutput->pFat_TDC_Singles_t[i] = Fat_TDC_Singles_t_Raw[i];    
         pOutput->pFat_TDC_Singles_ID[i] = Fat_TDC_Singles_ID[i];

         //cout << " i: " << i << " stdcmult " << stdcmult <<  " singles time " << Fat_TDC_Singles_t[i] << " singles tid " << Fat_TDC_Singles_ID[i] << endl;

      }

      for(int i = 0; i < sqdcmult; i++){
            
         pOutput->pFat_QDC_Singles_E[i] = Fat_QDC_Singles_E[i];    
         pOutput->pFat_QDC_Singles_E_Raw[i] = Fat_QDC_Singles_E_Raw[i];    
         pOutput->pFat_QDC_Singles_ID[i] = Fat_QDC_Singles_ID[i];
         pOutput->pFat_QDC_Singles_t_coarse[i] = Fat_QDC_Singles_t_coarse[i];    
         pOutput->pFat_QDC_Singles_t_fine[i] = Fat_QDC_Singles_t_fine[i];

         //cout << " i: " << i << " sqdcmult " << sqdcmult <<  " singles energy " << Fat_QDC_Singles_E[i] << " singles qid " << Fat_QDC_Singles_ID[i] << endl;

      }
    
    
    
    
    
    
    
    
}//end of fatima

    /**----------------------------------------------------------------------------------------------**/
    /**--------------------------------------  GALILEO  ---------------------------------------------**/
    /**----------------------------------------------------------------------------------------------**/
    void EventAnlProc::Make_Galileo_Histos()
    {
      hGAL_ESum = MakeTH1('I',"Germanium/Energy/Ge_ESum","Germanium Energy Sum",5000,0,5000);
      hGAL_ESum_SC41 = MakeTH1('I',"Germanium/Sum/Ge_ESum_SC41","Germanium Gated SC41 AddbackSum",5000,0,5000);
      hGAL_ESum_largerange_OF = MakeTH1('I',"Germanium/Energy/Ge_ESum_largerange_OF","Germanium Energy Sum (Overflow)",5000,0,5000);
      hGAL_ESum_largerange_PU = MakeTH1('I',"Germanium/Energy/Ge_ESum_largerange_PU","Germanium Energy Sum (Pileup)",5000,0,5000);
      hGAL_Hit_Pat = MakeTH1('I',"Germanium/Stats/Ge_Hit_Pat","Germanium Hit Pattern",36,0,36);
      
      hGAL_Chan_E_Mat = MakeTH2('D',"Germanium/Energy/Germanium_E_Mat","Germanium Energy-Energy Matrix",2500,0,10000,2500,0,10000);
        
      hGAL_Chan_E_vsDet = MakeTH2('D',"Germanium/Energy/Germanium_E_vsDetID","Germanium Energy vs Det ID Matrix",48,0,48,2500,0,10000);
     
      hGAL_AddbackSum = MakeTH1('I',"Germanium/Energy/Germanium_Addback","Germanium Addback Energy Sum",5000,0,5000);
    
      for (int i=0; i<GALILEO_MAX_DETS; i++)
      {
        for (int j = 0; j < GALILEO_CRYSTALS; j++)
        {
          hGAL_Chan_E[i][j] = MakeTH1('D',Form("GALILEO/Energy_Ch/GALILEO_E_Det_%2d_%1d",i, j), Form("GALILEO Channel Energy Detector %2d Crystal %1d",i, j),5000,0,5000);
                }
   
        }
      //for (int k=0; k<32; k++){
    
    
    
        //hGAL_Chan_Time_Diff[k] = MakeTH1('D',Form("GALILEO/Time_diff/GALILEO_Chan_Time_Diff%2d",k), Form("GALILEO Channel Time Difference for %2d",k),2000,-1000,1000);
        //hGAL_Chan_Timedifference_new[k] = MakeTH1('D',Form("GALILEO/Timediff_new/GALILEO_Chan_T_Diff%2d",k), Form("GALILEO Channel T Difference for %2d",k),2000,-1000,1000);
       // hGAL_Time_Diff_vs_Energy[k] = MakeTH2('D',Form("GALILEO/GALILEO_dT_vs_Energy_Spectra/GALILEO_dT_vs_E%2d",k), Form("GALILEO Time Difference Vs Channel Energy Channel %2d",k),5000,0,5000,100,-1000,1000);
      //}
        }
    ///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------///
    void EventAnlProc::Do_Galileo_Histos(EventAnlStore* pOutput)
    {
      
        pOutput->pGAL_WR=Gal_WR;
       
        // Process hits once
       
        for (int i = 0; i < GalFired; i++)
        {
              int det = GalDet[i];
              int crys = GalCrys[i];
              pOutput->pGalPileUp[i] = GalPileUp[i];
              pOutput->pGalOverFlow[i] = GalOverFlow[i];
           // Skip pileup/overflow events
           if (GalPileUp[i] && (det * GALILEO_CRYSTALS + crys)<24)
           {
              hGAL_ESum_largerange_PU->Fill(GalE_Cal[i]);
              continue;
   
            }
  
           if (GalOverFlow[i]&& (det * GALILEO_CRYSTALS + crys)<24)
           {
             hGAL_ESum_largerange_OF->Fill(GalE_Cal[i]);
             continue;
            }
  
          
	   if(det>-1){
           pOutput->pGal_Event_T[det][crys] = GalEventT[i];
          
           pOutput->pGal_T[det][crys] = GalT[i];
           pOutput->pGal_E_Raw[det][crys] = GalE[i];
           pOutput->pGal_E[det][crys] = GalE_Cal[i];
          
 
           hGAL_Hit_Pat->Fill(det * GALILEO_CRYSTALS + crys);
           hGAL_ESum->Fill(GalE_Cal[i]);
           
           hGAL_Chan_E[det][crys]->Fill(GalE_Cal[i]);
           hGAL_Chan_E_vsDet->Fill(det,GalE_Cal[i]);
    
           // 2D Matrix generation
           for (int j = 0; j < GalFired; j++)
           {
              if (i == j) continue;
              hGAL_Chan_E_Mat->Fill(GalE_Cal[i], GalE_Cal[j]);
                }
            }  
        }
        static const long dT_addback = 50;
        static const double energy_thres_addback = 250;
               
        // Detector addback
        for (int i = 0; i < GALILEO_MAX_DETS; i++)
        {
            double E[GALILEO_CRYSTALS] = { 0 };
            long T[GALILEO_CRYSTALS] = { 0 };
            int n[GALILEO_CRYSTALS] = { 0 };
            int v = 0;
            for (int j = 0; j < GALILEO_CRYSTALS; j++)
                {
              if (pOutput->pGal_E[i][j] == 0) continue;
              bool added = false; 
        
              // Try to addback to an existing hit
              for (int k = 0; k < v; k++)
              {
        if (T[k] - pOutput->pGal_T[i][j] < dT_addback){
            if(pOutput->pGal_E[i][j]>energy_thres_addback){
                
                  E[k] += pOutput->pGal_E[i][j];
                  T[k] = (T[k] + pOutput->pGal_T[i][j]) / (n[k] + 1);
                  n[k] += 1;
                  added = true;
                    }
                }
              }
               
              // Add to a new hit
              if (!added)
              {
                T[v] = pOutput->pGal_T[i][j];
                E[v] = pOutput->pGal_E[i][j];
                n[v] = 1;
                v++;
              }
            }         
            // Fill and write to Tree the addback energies
            for (int j = 0; j < v; j++)
            {
                pOutput->pGal_EAddback[i][j] = E[j];
               
                ///SC41 gated gal
             
           if(pOutput->pGal_EAddback[8][0]>0) hGAL_ESum_SC41->Fill(E[j]);
            if(i!=8)   hGAL_AddbackSum->Fill(E[j]);
                
            }
        }
    }//end of Do_Galileo_Histos()
//--------------------------------------------------------------------------------------------------------------------//
TH1I* EventAnlProc::MakeH1I(const char* fname,
                            const char* hname,
                            Int_t nbinsx,
                            Float_t xmin, Float_t xmax,
                            const char* xtitle,
                            Color_t linecolor,
                            Color_t fillcolor,
                            const char* ytitle) {
//    TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetHistogram, fname, hname);
//    if (res!=0) return dynamic_cast<TH1I*>(res);

   TH1I* histo = new TH1I(hname, hname, nbinsx, xmin, xmax);
   histo->SetXTitle(xtitle);
   if (ytitle) histo->SetYTitle(ytitle);
   histo->SetLineColor(linecolor);
   histo->SetFillColor(fillcolor);
   AddHistogram(histo, fname);
   return histo;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TH2I* EventAnlProc::MakeH2I(const char* fname,
                             const char* hname,
                             Int_t nbinsx, Float_t xmin, Float_t xmax,
                             Int_t nbinsy, Float_t ymin, Float_t ymax,
                             const char* xtitle, const char* ytitle,
                             Color_t markercolor) {
//    TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetHistogram, fname, hname);
//    if (res!=0) return dynamic_cast<TH2I*>(res);

   TH2I* histo = new TH2I(hname, hname, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
   histo->SetMarkerColor(markercolor);
   histo->SetXTitle(xtitle);
   histo->SetYTitle(ytitle);
   AddHistogram(histo, fname);
   return histo;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TGo4WinCond* EventAnlProc::MakeWindowCond(const char* fname,
                                           const char* cname,
                                           float left,
                                           float right,
                                           const char* HistoName) {
  // TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4WinCond*>(res);

   TGo4WinCond* cond = new TGo4WinCond((Text_t*)cname);
   cond->SetValues(left, right);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}
//-----------------------------------------------------------------------------------------------------------------------------//

TGo4PolyCond* EventAnlProc::MakePolyCond(const char* fname,
                                          const char* cname,
                                          Int_t size,
                                          Float_t (*points)[2],
                                          const char* HistoName) {
   //TNamed* res = TestObject((getfunc)&TGo4EventProcessor::GetAnalysisCondition, fname, cname);
   //if (res!=0) return dynamic_cast<TGo4PolyCond*>(res);

   Float_t fullx[size+1], fully[size+1];
   int numpoints = size;

   for (int i=0;i<numpoints;i++) {
     fullx[i] = points[i][0];
     fully[i] = points[i][1];
   }

   // connect first and last points
   if ((fullx[0]!=fullx[numpoints-1]) || (fully[0]!=fully[numpoints-1])) {
      fullx[numpoints] = fullx[0];
      fully[numpoints] = fully[0];
      numpoints++;
   }

   TCutG mycat("initialcut", numpoints, fullx, fully);
   TGo4PolyCond* cond = new TGo4PolyCond((Text_t*)cname);
   cond->SetValues(&mycat);
   cond->Enable();
   if (HistoName!=0)
     cond->SetHistogram(HistoName);
   AddAnalysisCondition(cond, fname);
   return cond;
}

void EventAnlProc::clear_gate_arrays(){
    ZAoQgnum=-1;
     for(int i=0;i<8; i++){
        for(int j=0;j<6; j++){
                XX2_AoQ[i][j]=0;
                YX2_AoQ[i][j]=0;
                XX4_AoQ[i][j]=0;
                YX4_AoQ[i][j]=0;
                X_ZZ2[i][j]=0;
                Y_ZZ2[i][j]=0;
                X_dEdeg[i][j]=0;
                Y_dEdeg[i][j]=0;
               
    
        }
        for(int j=0;j<8; j++){
         X_ZAoQ[i][j]=0;
         Y_ZAoQ[i][j]=0;
     }
    }
    
}
///-------------------------------------------------------------------------------------------------------
void EventAnlProc::FRS_Gates(){
  Int_t i;
  ifstream    file;
   file.open("Configuration_Files/2D_Gates/ID_x2AoQ.txt");
    if(!file.good()) cout<<"Configuration_Files/2D_Gates/ID_x2AoQ.txt Not found!"<<endl;
    for (i = 0; i < 8; i++){
         for (int j = 0; j < 6; j++){
       if(IsData(file)) file >>X2AoQgnum>> XX2_AoQ[i][j]>> YX2_AoQ[i][j] ;
         }
    }
  file.close();
  
  
 ///--------------------------------------------------------------------------------
  file.open("Configuration_Files/2D_Gates/ID_x4AoQ.txt");
    if(!file.good()) cout<<"Configuration_Files/2D_Gates/ID_x4AoQ.txt Not found!"<<endl;
    for (i = 0; i < 8; i++){
        for (int j = 0; j < 6; j++){
       if(IsData(file)) file >> X4AoQgnum>>XX4_AoQ[i][j]>> YX4_AoQ[i][j] ;
        }
    }
  file.close();
  
  
 ///--------------------------------------------------------------------------------
  
  file.open("Configuration_Files/2D_Gates/ID_Z_Z2.txt");
  if(!file.good()) cout<<"Configuration_Files/2D_Gates/ID_Z_Z2.txt Not found!"<<endl;
 for (i = 0; i < 8; i++){
    for (int j = 0; j < 6; j++){
       if(IsData(file)) file >>Z1Z2gnum>> X_ZZ2[i][j]>> Y_ZZ2[i][j] ;
    }
    }
  file.close();
  
  
 ///--------------------------------------------------------------------------------
      file.open("Configuration_Files/2D_Gates/ID_ZvsAoQ.txt");
    if(!file.good()) cout<<"Configuration_Files/2D_Gates/ID_ZvsAoQ.txt Not found!"<<endl;
   
    for (i = 0; i < 8; i++){
        for(int j=0; j<8; j++){
       if(IsData(file)) file >>ZAoQgnum >> X_ZAoQ[i][j]>> Y_ZAoQ[i][j] ;
 //   cout<<"ZAoQgnum " << ZAoQgnum << " X_ZAoQ[i][j] " <<X_ZAoQ[i][j]<< " Y_ZAoQ[i][j] " << Y_ZAoQ[i][j]<<" i " << i << " j " << j <<  endl;
        }
    }
  file.close();
  ///--------------------------------------------------------------------------------
      file.open("Configuration_Files/2D_Gates/ID_dEdeg_Z1.txt");
    
    for (i = 0; i < 8; i++){
        for(int j=0; j<6; j++){
       if(IsData(file)) file >>dEdeggnum >> X_dEdeg[i][j]>> Y_dEdeg[i][j] ;
        }
    }
  file.close();
  
  
}
///-------------------------------------------------------------------------------------------------------
 int EventAnlProc::IsData(ifstream &f) {
        char dum;
        char dumstr[300];
        int retval = 0;

        /* 'operator >>' does not read End-of-line, therefore check if read 
            character is not EOL (10) */
        do {
            dum = f.get();
            if (dum == '#' || dum==' ')    // comment line => read whole line and throw it away
            f.getline(dumstr,300);
        }
        while ((dum == '#') || ((int)dum == 10)|| dum==' '); 

        f.unget();   // go one character back
        retval = 1;
        return retval;
    }
///-------------------------------------------------------------------------------------------------------

void EventAnlProc::checkTAMEXorVME(){

  std::ifstream PL_FILE("Configuration_Files/DESPEC_General_Setup/TAMEX_or_VME.txt");

  std::string line;

  if(PL_FILE.fail()){
    std::cerr << "Could not find Configuration_Files/DESPEC_General_Setup/TAMEX_or_VME.txt file" << std::endl;
    exit(1);
  }
  bool T_or_V_bPlas = false;
  bool T_or_V_Fatima = false;
  bool T_and_V_Fatima = false;
  while(std::getline(PL_FILE,line)){
    if(line[0] == '#') continue;

    if(line == "VME_bPlas") T_or_V_bPlas = true;
    if(line == "TAMEX_bPlas") T_or_V_bPlas = false;

    if(line == "VME_Fatima") T_or_V_Fatima = true;
    if(line == "TAMEX_Fatima") T_or_V_Fatima = false;

    if(line == "VME_AND_TAMEX_Fatima") T_or_V_Fatima = false;
    if(line == "VME_AND_TAMEX_Fatima") T_and_V_Fatima = true;

    if(line == "VME_Fatima") T_and_V_Fatima = false;
    if(line == "TAMEX_Fatima") T_and_V_Fatima = false;


//     if(line != "VME_bPlas" && line != "TAMEX_bPlas"){
//       std::cerr << line << " module of PLASTIC not known!" <<std::endl;
//       exit(1);
//     }
  }

  VME_TAMEX_bPlas = T_or_V_bPlas;
  VME_TAMEX_Fatima = T_or_V_Fatima;
  VME_AND_TAMEX_Fatima = T_and_V_Fatima;

}
///-----------------------------------------------------------------------------------------------------------------///
   void EventAnlProc::get_used_systems(){
    for(int i = 0;i < 7;i++) Used_Systems[i] = false;

    ifstream data("Configuration_Files/DESPEC_General_Setup/Used_Systems.txt");
    if(data.fail()){
        cerr << "Could not find Used_Systems config file!" << endl;
        exit(0);
    }
    int i = 0;
    int id = 0;
    string line;
    char s_tmp[100];
    while(data.good()){
        getline(data,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),"%s %d",s_tmp,&id);
        Used_Systems[i] = (id == 1);
        i++;
    }
    string DET_NAME[7] = {"FRS","AIDA","PLASTIC","FATIMA_VME","FATIMA_TAMEX","GALILEO","FINGER"};

//     cout << "\n=====================================================" << endl;
//     cout << "USED SYSTEMS" << endl;
//     cout << "-----------------------------------------------------" << endl;
//     for(int j = 0;j < 6;++j){
//         if(Used_Systems[j]) cout << DET_NAME[j] << endl;
//     }
//     cout << "=====================================================" << endl;


}
//-----------------------------------------------------------------------------------------------------------------------------//
//                                                            END                                                              //
//-----------------------------------------------------------------------------------------------------------------------------//
