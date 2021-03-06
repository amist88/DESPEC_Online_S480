// $Id: EventUnpackProc.cxx 754 2011-05-18 11:04:52Z adamczew $
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
#include "EventUnpackProc.h"
#include "EventUnpackStore.h"
#include "Riostream.h"

// Root Includes //
#include "TROOT.h"
// #include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TTree.h"

#include <time.h>
#include <math.h>
#include <iomanip>

// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "Go4StatusBase/TGo4Picture.h"
#include "TGo4MbsEvent.h"

#include "TGo4MbsSubEvent.h"

// General Includes //
#include <fstream>
#include <vector>

#include "Detector_System.cxx"
#include "AIDA_Detector_System.h"
#include "FATIMA_Detector_System.h"
#include "FATIMA_TAMEX_Detector_System.h"
#include "PLASTIC_TAMEX_Detector_System.h"
//#include "PLASTIC_VME_Detector_System.h"
#include "FINGER_Detector_System.h"
#include "GALILEO_Detector_System.h"
#include "FRS_Detector_System.h"

#include "TAidaConfiguration.h"

#include "CalibParameter.h"
#include "CorrelParameter.h"

#include "White_Rabbit.h"

#include <string>


using namespace std;


//***********************************************************
EventUnpackProc::EventUnpackProc() :TGo4EventProcessor("Proc")
{
  cout << "**** EventUnpackProc: Create instance " << endl;


}
//***********************************************************
// standard factory
EventUnpackProc::EventUnpackProc(const char* name) : TGo4EventProcessor(name)
{


  cout << "**** EventUnpackProc: Create" << endl;

  //  input_data_path_old = "old";
   // WR_out.open ("WR_diff_store_270289.txt");
  WR_used = false;

  //used_systems
  get_used_systems();
  get_WR_Config();
    //read_setup_parameters();

  //  FAT_det_pos_setup();
  checkTAMEXorVME();
  
  //create White Rabbit obj
  WR = new White_Rabbit();

  fCal = (CalibParameter*) GetParameter("CalibPar"); 

 //create Detector Systems
  Detector_Systems = new Detector_System*[7];

 // all non used systems intialized as NULL
  //-> calling uninitialized system will cause an error !

  Detector_Systems[0] = !Used_Systems[0] ? nullptr : new FRS_Detector_System();
  Detector_Systems[1] = !Used_Systems[1] ? nullptr : new AIDA_Detector_System();
  
//  if(VME_TAMEX_bPlas==true) Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_VME_Detector_System();
  if(VME_TAMEX_bPlas==false) Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_TAMEX_Detector_System();

  if(VME_TAMEX_Fatima==true || VME_AND_TAMEX_Fatima==true) Detector_Systems[3] =  new FATIMA_Detector_System();


  if(VME_TAMEX_Fatima==false|| VME_AND_TAMEX_Fatima==true) Detector_Systems[4] = new FATIMA_TAMEX_Detector_System();

  Detector_Systems[5] = !Used_Systems[5] ? nullptr : new GALILEO_Detector_System();
  Detector_Systems[6] = !Used_Systems[6] ? nullptr : new FINGER_Detector_System();


  //Only create histograms if system is used
  if(Used_Systems[0]) Make_FRS_Histos();

  if(Used_Systems[1]) Make_AIDA_Histos();


  if(Used_Systems[4] && VME_TAMEX_Fatima==false && VME_AND_TAMEX_Fatima==false) Make_FATIMA_TAMEX_Histos();

  if(Used_Systems[3] && VME_TAMEX_Fatima==true && VME_AND_TAMEX_Fatima==false)  Make_FATIMA_VME_Histos();


 // if((Used_Systems[3] || Used_Systems[4]) && VME_TAMEX_Fatima==false && VME_AND_TAMEX_Fatima==true) Make_FATIMA_VME_TAMEX_Histos();

  if(Used_Systems[5]) Make_GALILEO_Histos();

  RAW = new Raw_Event();

  load_PrcID_File();

  load_FingerID_File();
  load_FatTamex_Allocationfile();
  read_setup_parameters();

    WR_count = 0;
    count = 0;
    iterator = 0;
    val_it = 0;
 /// zero FRS scalers (cumulative)
    memset(frs_scaler_value, 0, sizeof(frs_scaler_value));
  ///Clear for AIDA
  lastTime = 0;
  ID = 0;
  totalEvents = 0;
  startTime = 0;
  stopTime = 0;
  fAida.ImplantEvents.clear();
  fAida.DecayEvents.clear();
  fAida.Implants.clear();
  fAida.Decays.clear();
  /// Setup AIDA arrays
  if(Used_Systems[1])
  {
    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
    adcLastTimestamp.resize(conf->FEEs());
    adcCounts.resize(conf->FEEs());
  }

}

void EventUnpackProc::UserPostLoop()
{
  if (Used_Systems[1])
  {
    ((AIDA_Detector_System*)Detector_Systems[1])->PrintStatistics();
  }
}


//----------------------------------------------------------


EventUnpackProc::~EventUnpackProc()
{

  delete[] Detector_Systems;
  delete RAW;
  delete WR;
  cout << "**** EventUnpackProc: Delete instance" << endl;
}

//----------------------------------------------------------
Bool_t EventUnpackProc::BuildEvent(TGo4EventElement* dest)
{
  bool skip;
  static TString oldfname = "";
  TGo4MbsEvent *fMbsEvent = dynamic_cast<TGo4MbsEvent*>    (GetInputEvent("Unpack"));
// s_filhe* fileheader=fMbsEvent->GetMbsSourceHeader();
  //s_bufhe* head = GetMbsBufferHeader();

  EventUnpackStore* fOutput = (EventUnpackStore*) dest;
  TGo4MbsEvent* fInput = (TGo4MbsEvent*) GetInputEvent();

  if(fOutput==0){
    cout << "UnpackProc: no unpack output event !";
    return false;
  }
 // input_data_path = fileheader->filhe_file;

  count++;

  if (count % 100000 == 0){

    cout << "\r";
    cout << "Event " << count << " Reached!!!"<<"    Data File Number : "<<data_file_number;
    cout <<"\t\t\t\t";
    cout.flush();
  }

  Bool_t isValid=kFALSE; // validity of output event //

  if (fInput==0) // Ensures that there is data in the event //
  {
    cout << "EventUnpackProc: no input event !"<< endl;
    fOutput->SetValid(isValid);
    return isValid;
  }
  isValid=kTRUE;
  event_number=fInput->GetCount();
  fOutput-> fevent_number = event_number;
  fOutput->fTrigger = fInput->GetTrigger();

  fInput->ResetIterator();
  TGo4MbsSubEvent* psubevt(0);


  // ------------------------------------------------------ //
  // |                                                    | //
  // |               START OF EVENT ANALYSIS              | //
  // |                                                    | //
  // ------------------------------------------------------ //

 if(true){
  //if (event_number==99833){
  //cout<<"event " << event_number <<endl;
      int subevent_iter = 0;
      Int_t PrcID_Conv = 0;

      Int_t* pdata = nullptr;
      Int_t lwords = 0;
      Int_t PrcID = -1;
      Int_t sub_evt_length = 0;
      Int_t Type =-1;
      Int_t SubType =-1;
      WR_tmp = 0;
      WR_d=0;
      WR_main=0;

      while ((psubevt = fInput->NextSubEvent()) != 0) // subevent loop //
      {
        subevent_iter++;
        pdata = psubevt->GetDataField();
        lwords = psubevt->GetIntLen();
        PrcID = psubevt->GetProcid();
        Type = psubevt->GetType();
        SubType = psubevt->GetSubtype(); 
        
        PrcID_Conv = get_Conversion(PrcID);
        
	if(PrcID_Conv==-1) continue;
        fOutput -> fProcID[PrcID_Conv] = PrcID_Conv;
        sub_evt_length  = (psubevt->GetDlen() - 2) / 2;

    ///------------------------------WHITE RABBIT --------------------------------------////
        if(WHITE_RABBIT_USED){
          //sub_evt_length = sub_evt_length - 5;

            //Pulls it straight from White_Rabbit class
            WR_tmp = WR->get_White_Rabbit(pdata);
            WR_d = WR->get_Detector_id();
            pdata = WR->get_pdata();
            
            ///Temp Wr detector fix
            if(PrcID ==10|| PrcID == 30 || PrcID==20 || PrcID == 25 || PrcID==41 ||PrcID==100)WR_d=0;

           if(WR_d==0) fOutput->fFRS_WR = WR_tmp; //FRS
           if(WR_d==1) fOutput->fAIDA_WR = WR_tmp; //AIDA
           if(WR_d==2) fOutput->fbPlas_WR = WR_tmp; //bPlas (TAMEX)
           if(WR_d==3) fOutput->fFat_WR = WR_tmp; //Fatima (VME)
           if(WR_d==5) fOutput->fGal_WR = WR_tmp; //Galileo
           if(WR_d==6) fOutput->fFinger_WR = WR_tmp; //FINGER
          
            WR_main = WR_tmp;
            
            
        }

///-----------------------------------------------------------------------------------------------------------///
        //if necessary, directly print MBS for wanted Detector_System
//         if(PrcID_Conv == AIDA && false) print_MBS(pdata,lwords);
//         if(PrcID_Conv == FATIMA && false) print_MBS(pdata,lwords);
//         if(PrcID_Conv == PLASTIC && false) print_MBS(pdata,lwords);
//         if(PrcID_Conv == GALILEO && false) print_MBS(pdata,lwords);
        // if(PrcID_Conv == FINGER && false) print_MBS(pdata,lwords);

        //=================================================================
        //UNPACKING
        ///send subevent to respective unpacker

        if(Detector_Systems[PrcID_Conv] !=0){
        Detector_Systems[PrcID_Conv]->Process_MBS(psubevt);
        Detector_Systems[PrcID_Conv]->Process_MBS(pdata);
       // cout<<"Detector_Systems[PrcID_Conv] " <<Detector_Systems[PrcID_Conv] <<" PrcID_Conv " <<PrcID_Conv<<endl;
        ///get mbs stream data from unpacker (pointer copy solution)
        pdata = Detector_Systems[PrcID_Conv]->get_pdata();

        ///get data from subevent and send to RAW
        Detector_Systems[PrcID_Conv]->get_Event_data(RAW);
        }


        //=================================================================
        //HISTOGRAM FILLING (only singles)
        FILL_HISTOGRAMS(PrcID_Conv,PrcID,SubType);
        //=================================================================

        pdata = nullptr;

        ///--------------------------------------------------------------------------------------------///
                                /** Unpack Tree for each detector subsystem**/
        ///--------------------------------------------------------------------------------------------///
                                                /** Output FRS **/
        ///--------------------------------------------------------------------------------------------///

   if (Used_Systems[0] && PrcID_Conv==0){

        ///MUSIC
          if(PrcID==20){
           for(int i =0; i<2; ++i){
            fOutput->fFRS_Music_dE[i] = RAW->get_FRS_MusicdE(i);
            fOutput->fFRS_Music_dE_corr[i] = RAW->get_FRS_MusicdE_corr(i);
         //  cout<<"fOutput->fFRS_Music_dE[i] " <<fOutput->fFRS_Music_dE[i] << " i " << i << endl;
           }
           for(int i =0; i<8; ++i){
        fOutput->fFRS_Music_E1[i] = RAW->get_FRS_MusicE1(i);
        fOutput->fFRS_Music_E2[i] = RAW->get_FRS_MusicE2(i);
        fOutput->fFRS_Music_T1[i] = RAW->get_FRS_MusicT1(i);
        fOutput->fFRS_Music_T2[i] = RAW->get_FRS_MusicT2(i);
//         if(fOutput->fFRS_Music_E1[i]!=0)cout<<"fOutput->fFRS_Music_E1[i] " <<fOutput->fFRS_Music_E1[i] << " i " << i << endl;
           }
          }
        ///SCI
        
           for(int l=0;l<12;++l){
            if(PrcID==10){
            if(RAW->get_FRS_sci_l(l)!=0) fOutput->fFRS_sci_l[l] = RAW->get_FRS_sci_l(l);
            if(RAW->get_FRS_sci_r(l)!=0) fOutput->fFRS_sci_r[l] = RAW->get_FRS_sci_r(l);
                    }
                if(PrcID==30){
            if(RAW->get_FRS_sci_e(l)!=0) fOutput->fFRS_sci_e[l] = RAW->get_FRS_sci_e(l);
            if(RAW->get_FRS_sci_tx(l)!=0) fOutput->fFRS_sci_tx[l] = RAW->get_FRS_sci_tx(l);
            if(RAW->get_FRS_sci_x(l)!=0)  fOutput->fFRS_sci_x[l] = RAW->get_FRS_sci_x(l);
            }
        }
        
     if(PrcID==35){
        if(RAW->get_FRS_TRaw_vftx_21l()!=0)fOutput->fTRaw_vftx_21l = RAW->get_FRS_TRaw_vftx_21l();
        if(RAW->get_FRS_TRaw_vftx_21r()!=0)fOutput->fTRaw_vftx_21r = RAW->get_FRS_TRaw_vftx_21r();
        if(RAW->get_FRS_TRaw_vftx_22l()!=0)fOutput->fTRaw_vftx_22l = RAW->get_FRS_TRaw_vftx_22l();
        if(RAW->get_FRS_TRaw_vftx_22r()!=0)fOutput->fTRaw_vftx_22r = RAW->get_FRS_TRaw_vftx_22r();
        if(RAW->get_FRS_TRaw_vftx_41l()!=0)fOutput->fTRaw_vftx_41l = RAW->get_FRS_TRaw_vftx_41l();
        if(RAW->get_FRS_TRaw_vftx_41r()!=0)fOutput->fTRaw_vftx_41r = RAW->get_FRS_TRaw_vftx_41r();
        if(RAW->get_FRS_TRaw_vftx_42l()!=0)fOutput->fTRaw_vftx_42l = RAW->get_FRS_TRaw_vftx_42l();
        if(RAW->get_FRS_TRaw_vftx_42r()!=0)fOutput->fTRaw_vftx_42r = RAW->get_FRS_TRaw_vftx_42r();
             }
           
            ///SCI TOF
//         fOutput->fFRS_sci_tofll2 = RAW->get_FRS_tofll2();
//         fOutput->fFRS_sci_tofll3 = RAW->get_FRS_tofll3();
//         fOutput->fFRS_sci_tof2 = RAW->get_FRS_tof2();
//         fOutput->fFRS_sci_tofrr2 = RAW->get_FRS_tofrr2();
//         fOutput->fFRS_sci_tofrr3 = RAW->get_FRS_tofrr3();
//         fOutput->fFRS_sci_tof3 = RAW->get_FRS_tof3();
        ///ID 2 4
   if(RAW->get_FRS_x2()!=0)      fOutput->fFRS_ID_x2 = RAW->get_FRS_x2();      
   if(RAW->get_FRS_y2()!=0)      fOutput->fFRS_ID_y2 = RAW->get_FRS_y2();  
   if(RAW->get_FRS_a2()!=0)      fOutput->fFRS_ID_a2 = RAW->get_FRS_a2();
   if(RAW->get_FRS_b2()!=0)      fOutput->fFRS_ID_b2 = RAW->get_FRS_b2();     
   
   if(RAW->get_FRS_x4()!=0)      fOutput->fFRS_ID_x4 = RAW->get_FRS_x4();
   if(RAW->get_FRS_y4()!=0)      fOutput->fFRS_ID_y4 = RAW->get_FRS_y4();
   if(RAW->get_FRS_a4()!=0)      fOutput->fFRS_ID_a4 = RAW->get_FRS_a4();
   if(RAW->get_FRS_b4()!=0)      fOutput->fFRS_ID_b4 = RAW->get_FRS_b4();
            ///SCI dT
//         fOutput->fFRS_sci_dt_21l_21r = RAW->get_FRS_dt_21l_21r();
//         fOutput->fFRS_sci_dt_41l_41r = RAW->get_FRS_dt_41l_41r();
//         fOutput->fFRS_sci_dt_42l_42r = RAW->get_FRS_dt_42l_42r();
//         fOutput->fFRS_sci_dt_43l_43r = RAW->get_FRS_dt_43l_43r();
// 
//         fOutput->fFRS_sci_dt_21l_41l = RAW->get_FRS_dt_21l_41l();
//         fOutput->fFRS_sci_dt_21r_41r = RAW->get_FRS_dt_21r_41r();
// 
//         fOutput->fFRS_sci_dt_21l_42l = RAW->get_FRS_dt_21l_42l();
//         fOutput->fFRS_sci_dt_21r_42r = RAW->get_FRS_dt_21r_42r();
            ///ID B Rho
         for(int i =0; i<2; ++i)
            fOutput->fFRS_ID_brho[i] = RAW->get_FRS_brho(i);
//             fOutput->fFRS_ID_rho[i] = RAW->get_FRS_rho(i);
//         }
    
      ///Using TAC
//        
//         fOutput->fFRS_beta3 = RAW->get_FRS_beta3();
//        fOutput->fFRS_gamma = RAW->get_FRS_gamma();
//        fOutput->fFRS_beta = RAW->get_FRS_beta();
//         fOutput->fFRS_AoQ = RAW->get_FRS_AoQ();
//         fOutput->fFRS_AoQ_corr = RAW->get_FRS_AoQ_corr();
//         fOutput->fFRS_z = RAW->get_FRS_z();
//         fOutput->fFRS_z2 = RAW->get_FRS_z2();
//         fOutput->fFRS_znocorr =  RAW->get_FRS_znocorr();
//         fOutput->fFRS_z2nocorr =  RAW->get_FRS_z2nocorr();
//         fOutput->fFRS_dEdeg = RAW->get_FRS_dEdeg();
//         fOutput->fFRS_dEdegoQ = RAW->get_FRS_dEdegoQ();
        
        ///Using MHTDC
         if(RAW->get_FRS_id_mhtdc_aoq()>0) fOutput->fFRS_AoQ = RAW->get_FRS_id_mhtdc_aoq();
         if(RAW->get_FRS_id_mhtdc_aoq_corr()>0)fOutput->fFRS_AoQ_corr = RAW->get_FRS_id_mhtdc_aoq_corr();
         if(RAW->get_FRS_id_mhtdc_z1()>0) fOutput->fFRS_z = RAW->get_FRS_id_mhtdc_z1();
         if(RAW->get_FRS_id_mhtdc_z2()>0) fOutput->fFRS_z2 = RAW->get_FRS_id_mhtdc_z2();
        fOutput->fFRS_znocorr =  RAW->get_FRS_znocorr();
        fOutput->fFRS_z2nocorr =  RAW->get_FRS_z2nocorr();
        if(RAW->get_FRS_id_mhtdc_dEdeg()>0)fOutput->fFRS_dEdeg = RAW->get_FRS_id_mhtdc_dEdeg();
        if(RAW->get_FRS_id_mhtdc_dEdegoQ()>0)fOutput->fFRS_dEdegoQ = RAW->get_FRS_id_mhtdc_dEdegoQ();
        if(RAW->get_FRS_id_mhtdc_beta()>0.0 &&  RAW->get_FRS_id_mhtdc_beta()<1.0) fOutput->fFRS_beta = RAW->get_FRS_id_mhtdc_beta();
        
      
        
      
     
        for (int i = 0; i < 7; i++){
            if(PrcID!=35) {
            fOutput->fFRS_TPC_x[i]=0;
            fOutput->fFRS_TPC_y[i]=0;
        }
             if(PrcID==35 && SubType==3600){
        
           fOutput->fFRS_TPC_x[i] = RAW->get_FRS_tpcX(i);
           fOutput->fFRS_TPC_y[i] = RAW->get_FRS_tpcY(i);
       
        }
       
        }
       
   
        
        for (int i = 0; i < 64; i++)
        {
          fOutput->fFRS_scaler[i] = frs_scaler_value[i];
          fOutput->fFRS_scaler_delta[i] = increase_scaler_temp[i];
        }
        //fOutput->fFRS_z3 = RAW->get_FRS_z3();
            ///ID Timestamp
//         fOutput->fFRS_timestamp = RAW->get_FRS_timestamp();
//         fOutput->fFRS_ts = RAW->get_FRS_ts();
//         fOutput->fFRS_ts2 = RAW->get_FRS_ts2();
         }
         ///--------------------------------------------------------------------------------------------///
                                            /** Output AIDA **/
        ///--------------------------------------------------------------------------------------------///

        if (Used_Systems[1] && PrcID_Conv==1){
          TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
          AIDA_Hits = RAW->get_AIDA_HITS();

          AidaEvent evt;
          fOutput->fAIDAHits += AIDA_Hits;

          for(int i = 0; i<AIDA_Hits; i++){

            AIDA_Energy[i] = RAW->get_AIDA_Energy(i);
            AIDA_FEE[i] = RAW-> get_AIDA_FEE_ID(i);
            AIDA_ChID[i] = RAW-> get_AIDA_CHA_ID(i);
            AIDA_Time[i] = RAW-> get_AIDA_WR(i);
            AIDA_HighE_veto[i] = RAW-> get_AIDA_HighE_VETO(i);
            AIDA_Side[i] = RAW-> get_AIDA_SIDE(i);
            AIDA_Strip[i] = RAW-> get_AIDA_STRIP(i);
            AIDA_evtID[i] = RAW-> get_AIDA_EVTID(i);

            evt.Channel = AIDA_ChID[i];
            evt.Module = AIDA_FEE[i];
            evt.Time = AIDA_Time[i];
            evt.HighEnergy =  AIDA_HighE_veto[i];
            evt.DSSD = conf->FEE(AIDA_FEE[i]).DSSD;
            evt.Side = AIDA_Side[i];
            evt.Strip = AIDA_Strip[i];
            evt.ID = AIDA_evtID[i];
            evt.Energy = AIDA_Energy[i];
            evt.FastTime = RAW->get_AIDA_FastTime(i);
            //if(evt.HighEnergy>0){
       // cout<<"Event AIDA " << event_number << " evt.Energy " << evt.Energy <<" evt.HighEnergy " << evt.HighEnergy <<endl;}
            if (!startTime) startTime = evt.Time;
            stopTime = evt.Time;
            /// Build events from everything until there's a gap of 2000 �s (event window)

            /// If lastTime is 0 it's the first event
            /// New event for timewarps
            if (lastTime > 0 && (evt.Time - lastTime > conf->EventWindow()))
            {
              // if event happened too late, redo the event again with a new out_event
              lastTime = 0;
              ResetMultiplexer();

              totalEvents++;

              fOutput->Aida.push_back(fAida);
              fAida.ImplantEvents.clear();
              fAida.DecayEvents.clear();
              fAida.AIDATime = 0;
            }

            lastTime = evt.Time;
            CorrectTimeForMultiplexer(evt);


            if (evt.HighEnergy==1)
            {
              fAida.ImplantEvents.push_back(evt);
            }
            else
            {
              fAida.DecayEvents.push_back(evt);
            }

            if (fAida.AIDATime == 0)
            {
              fAida.AIDATime = evt.Time;
            }
          }

          if (conf->ucesb())
          {
            lastTime = 0;
            ResetMultiplexer();
            totalEvents++;

            fOutput->Aida.push_back(fAida);
            fAida.ImplantEvents.clear();
            fAida.DecayEvents.clear();
            fAida.AIDATime = 0;
          }
        }
       
     ///--------------------------------------------------------------------------------------------///
                                                /**Output bPLASTIC TAMEX and FATIMA **/
        ///--------------------------------------------------------------------------------------------///
        int Fatfired[4];
        int bPlasfired[2];
        int Phys_Channel_Lead_Fat[4][256];
        int Phys_Channel_Trail_Fat[4][256];
        int Phys_Channel_Lead_bPlas[2][256];
        int Phys_Channel_Trail_bPlas[2][256];
     if (Used_Systems[2]&& PrcID_Conv==2 && VME_TAMEX_bPlas == false){
          

          for (int i=0; i<RAW->get_PLASTIC_tamex_hits(); i++){///Loop over tamex ID's
         ///--------------------------------------------------------------------------------------------///
                                                /**Output FATIMA TAMEX **/
        ///--------------------------------------------------------------------------------------------///
               
            if(RAW->get_PLASTIC_TAMEX_ID(i) >1){
               
                Fatfired[i] = RAW->get_PLASTIC_am_Fired(i);

            for(int j = 0;j < Fatfired[i];j++){

              if(RAW->get_PLASTIC_CH_ID(i,j) % 2 == 1){ //Lead j
	
                Phys_Channel_Lead_Fat[i][j] =TAMEX_bPlasFat_ID[i][RAW->get_PLASTIC_physical_channel(i, j)]; //From allocation file in future
                int chan_fat = Phys_Channel_Lead_Fat[i][j];
                fOutput->fFat_chan = chan_fat;
	      
               int N1 = fOutput->fFat_PMT_Lead_N[chan_fat]++;
            fOutput->fFat_Lead_PMT[chan_fat][N1] = RAW->get_PLASTIC_lead_T(i,j);

              }
              else{ //Trail even j
                Phys_Channel_Trail_Fat[i][j] = TAMEX_bPlasFat_ID[i][RAW->get_PLASTIC_physical_channel(i, j)];
		
                int chan_fat = Phys_Channel_Trail_Fat[i][j];
            
                // PMT allocation succeeded
		if(chan_fat>-1){

              int N1 = fOutput->fFat_PMT_Trail_N[chan_fat]++;
	   
                fOutput->fFat_Trail_PMT[chan_fat][N1] = RAW->get_PLASTIC_trail_T(i,j);
		
                                }
                            }       
                        }
                     }
            
            
        
        ///--------------------------------------------------------------------------------------------///
                                                /**Output bPLASTIC TAMEX  **/
        ///--------------------------------------------------------------------------------------------///
      
        
      if(RAW->get_PLASTIC_TAMEX_ID(i)<2){
                                 
            int chan=-1;

            bPlasfired[i] = RAW->get_PLASTIC_am_Fired(i); ///Iterator
            int bPlasdetnum=i+1;
            for(int j = 0;j < bPlasfired[i];j++){

              if(RAW->get_PLASTIC_CH_ID(i,j) % 2 == 1){ ///Lead odd j
                  
                Phys_Channel_Lead_bPlas[i][j] = RAW->get_PLASTIC_physical_channel(i, j); //This fixes a problem with the first two channels
                    if(RAW->get_PLASTIC_CH_ID(i,j)==1)Phys_Channel_Lead_bPlas[0][j]=0;
                   
                chan = (Phys_Channel_Lead_bPlas[i][j])-i*16;
                fOutput->fbPlaschan=  chan; 
		

                // PMT allocation succeeded
                int N1 = fOutput->fbPlas_PMT_Lead_N[bPlasdetnum][chan]++;
            ///Set some limits to avoid strange behaviour
                if(N1>-1 && N1<10){
             //   cout<<"1 "<< " bPlasdetnum " << bPlasdetnum << " chan " << chan << " N1 " << N1 << endl;
                
                fOutput->fbPlas_Lead_PMT[bPlasdetnum][chan][N1] = RAW->get_PLASTIC_lead_T(i,j);
            
               // cout<<"2 fOutput->fbPlas_Lead_PMT[bPlasdetnum][chan][N1] " << fOutput->fbPlas_Lead_PMT[bPlasdetnum][chan][N1] << " bPlasdetnum " << bPlasdetnum << " chan " << chan << " N1 " << N1 << endl;
               }
              }
              
               if(RAW->get_PLASTIC_CH_ID(i,j) % 2 == 0){ ///Trail even j
                  
                Phys_Channel_Trail_bPlas[i][j] = RAW->get_PLASTIC_physical_channel(i,j);
                chan = (Phys_Channel_Trail_bPlas[i][j])-i*16;
		 
                // PMT allocation succeeded
                
                  
                int N1 = fOutput->fbPlas_PMT_Trail_N[bPlasdetnum][chan]++;
                ///Set some limits to avoid strange behaviour
                if(N1>-1 && N1<10){
         
             fOutput->fbPlas_Trail_PMT[bPlasdetnum][chan][N1] = RAW->get_PLASTIC_trail_T(i,j);  
           
                }            
              }
            }
          }
        }
       }
     
           ///--------------------------------------------------------------------------------------------///
                                                /**Output FATIMA VME **/
        ///--------------------------------------------------------------------------------------------///

        int Fat_QDC_ID;
        int Fat_TDC_ID_sing;
        int Fat_TDC_ID[48];
        int Fat_TDC_multi[50];
        bool TimID[50];
        bool EnID[50];
        int counter = 0;
        Double_t dummy_qdc_E[50];
        Double_t dummy_qdc_E_raw[50];

        Double_t dummy_tdc_t[50];
        Double_t dummy_tdc_t_raw[50];


		Long64_t dummy_qdc_t_coarse[50];
        Double_t dummy_qdc_t_fine[50];

		Long64_t dum_qdc_t_coarse[50];
        Double_t dum_qdc_t_fine[50];
        
        int dummy_qdc_id[50];
        int dummy_tdc_id[50];
        
        Double_t dum_qdc_E[50];
        Double_t dum_qdc_E_raw[50];

        Double_t dum_tdc_t[50];
        Double_t dum_tdc_t_raw[50];

        int dum_qdc_id[50];
        int dum_tdc_id[50];
        
        int dummytdcmult = 0;
        int dummyqdcmult = 0;

        int matchedmult = 0;
        
        int sc40count = 0;
        int sc41count = 0;
        
        int fatbplascount=0;
        
        int singlestdcmult = 0;
        int singlesqdcmult = 0;
                
        bool tdc_multi_hit_exclude[50];
        bool qdc_multi_hit_exclude[50];
        
        for (int i = 0; i<50; i++){
          Fat_TDC_multi[i] = 0;
          dummy_qdc_id[i] = -1;
          dummy_tdc_id[i] = -2;
          dummy_qdc_E[i] = 0;
          dummy_qdc_E_raw[i] = 0;
          dummy_qdc_t_coarse[i] = 0;
          dummy_qdc_t_fine[i] = 0;

          dummy_tdc_t[i] = 0;
          dummy_tdc_t_raw[i] = 0;
          dum_qdc_E[i] = 0;
          dum_qdc_E_raw[i] = 0;

          dum_tdc_t[i] = 0;
          dum_tdc_t_raw[i] = 0;

          dum_qdc_id[i] = -1;
          dum_tdc_id[i] = -2;
          
          dummy_qdc_t_coarse[i] = 0;
          dummy_qdc_t_fine[i] = 0;       
          
          qdc_multi_hit_exclude[i] = 0;
          tdc_multi_hit_exclude[i] = 0;
          TimID[i] = 0;
          EnID[i] = 0;

        }


        
                if (Used_Systems[3]&& PrcID_Conv==3){
            
        for (int i=0; i< RAW->get_FAT_TDCs_fired(); i++){
            for(int j = 0; j < RAW->get_FAT_TDCs_fired(); j++){
                    if(RAW->get_FAT_TDC_id(i) == RAW->get_FAT_TDC_id(j) && i!=j){
                        
                        tdc_multi_hit_exclude[RAW->get_FAT_TDC_id(i)] = true;

                    }
                }//end j loop
            }//end i loop
            
        for (int i=0; i< RAW->get_FAT_TDCs_fired(); i++){


         if(RAW->get_FAT_TDC_id(i) > -1 && RAW->get_FAT_TDC_id(i) < 36){
            fOutput->fFat_TDC_Singles_t[singlestdcmult] = RAW->get_FAT_TDC_timestamp(i);
            fOutput->fFat_TDC_Singles_t_Raw[singlestdcmult] = RAW->get_FAT_TDC_timestamp_raw(i);
            fOutput->fFat_TDC_Singles_ID[singlestdcmult] = RAW->get_FAT_TDC_id(i);

            singlestdcmult++;
          
                //cout << "Event no: " << event_number << " tdc time: " << dummy_tdc_t[i] << " tdc id: " << dummy_tdc_id[i] << endl;

                if(RAW->get_FAT_TDC_timestamp(i) != 0. && tdc_multi_hit_exclude[RAW->get_FAT_TDC_id(i)] == 0){
                    dummy_tdc_t[dummytdcmult] = RAW->get_FAT_TDC_timestamp(i);
                    dummy_tdc_t_raw[dummytdcmult] = RAW->get_FAT_TDC_timestamp_raw(i);
                    dummy_tdc_id[dummytdcmult] =  RAW->get_FAT_TDC_id(i);  
                    dummytdcmult++;
                    //cout << "Event no: " << event_number << " tdc time: " << dummy_tdc_t[i] << " tdc id: " << dummy_tdc_id[i] << endl;
                }//end if good event
                //cout << dummy_tdc_id[i] << endl;
			}//End if correct ID range
            
            if(RAW->get_FAT_TDC_id(i) == 40 && RAW->get_FAT_TDC_timestamp(i) != 0. ){
                fOutput->fSC40[sc40count] =  RAW->get_FAT_TDC_timestamp(i);
                sc40count++;
               
            }// end if statement to check for SC40 hits
            
            if(RAW->get_FAT_TDC_id(i) == 41 && RAW->get_FAT_TDC_timestamp(i) !=0.){
                fOutput->fSC41[sc41count] =  RAW->get_FAT_TDC_timestamp(i);
                sc41count++;
            }// end if statement to check for SC41 hits
            
            ///Get the bPlas coincident signal
//              if(RAW->get_FAT_TDC_id(i) == 38 && RAW->get_FAT_TDC_timestamp(i) !=0.){
//                 fOutput->fFat_bPlas[fatbplascount] =  RAW->get_FAT_TDC_timestamp(i);
//                 fatbplascount++;
//                 cout<<"RAW->get_FAT_TDC_timestamp(i) " <<RAW->get_FAT_TDC_timestamp(i) << endl;
//             }// end if statement to check for bPlas hits
            
        }// end dummy TDC for loop
        
        
        
        for (int i=0; i< RAW->get_FAT_QDCs_fired(); i++){
            
            for(int j = 0; j < RAW->get_FAT_QDCs_fired(); j++){
                if(RAW->get_FAT_QDC_id(i) == RAW->get_FAT_QDC_id(j) && i!=j){
                    qdc_multi_hit_exclude[RAW->get_FAT_QDC_id(i)] = true;
                }
            }//end j loop
        }
            
        for (int i=0; i< RAW->get_FAT_QDCs_fired(); i++){
		
          if(RAW->get_FAT_QDC_id(i) > -1 && RAW->get_FAT_QDC_id(i) < 36){
            fOutput->fFat_QDC_Singles_E[singlesqdcmult] = RAW->get_FAT_QLong(i);
            fOutput->fFat_QDC_Singles_ID[singlesqdcmult] = RAW->get_FAT_QDC_id(i);
            fOutput->fFat_QDC_Singles_t_coarse[singlesqdcmult] = RAW->get_FAT_QDC_t_Coarse(i);
            fOutput->fFat_QDC_Singles_t_fine[singlesqdcmult] = RAW->get_FAT_QDC_t_Fine(i);
            fOutput->fFat_QDC_Singles_E_Raw[singlesqdcmult] = RAW->get_FAT_QLong_Raw(i);


            singlesqdcmult++;
          
		
            if(RAW->get_FAT_QLong(i) > 10. && qdc_multi_hit_exclude[RAW->get_FAT_QDC_id(i)] == 0 ){
                dummy_qdc_E[dummyqdcmult] = RAW->get_FAT_QLong(i);
                dummy_qdc_E_raw[dummyqdcmult] = RAW->get_FAT_QLong_Raw(i);
                dummy_qdc_t_coarse[dummyqdcmult] = RAW->get_FAT_QDC_t_Coarse(i);
                dummy_qdc_t_fine[dummyqdcmult] = RAW->get_FAT_QDC_t_Fine(i);

                dummy_qdc_id[i] =  RAW->get_FAT_QDC_id(i);  
                dummyqdcmult++;
                //cout << "Event no: " << event_number << " qdc e: " << dummy_qdc_E[i] << " qdc id: " << dummy_qdc_id[i] << endl;

                //cout << RAW->get_FAT_QDC_id(i) << endl;
                //cout << dummy_qdc_id[i] << endl;

            }//end if to ensure real energies
		}//End if right IDs
	
        }// end dummy QDC for loop
        
        
        if(dummyqdcmult < dummytdcmult) counter = dummytdcmult;
        else if(dummytdcmult < dummyqdcmult) counter = dummyqdcmult;
        else if(dummyqdcmult == dummytdcmult) counter = dummytdcmult;


//cout << "event no: " << event_number << endl;

    
          for (int i=0; i< dummyqdcmult; i++){
                //cout << "i: " << i << " QDCID: " << dummy_qdc_id[i] << " TDCID " << dummy_tdc_id[i] << " energy " << RAW->get_FAT_QLong(i) << " timestamp" << RAW->get_FAT_TDC_timestamp(i) << endl;
              
              //cout << "i " << i << endl;
                if(dummy_qdc_id[i] == dummy_tdc_id[i] && EnID[i] == 0 && TimID[i] == 0 ){
                    EnID[i] = true;
                    TimID[i] = true;
                    
                                        
                    //Fat_TDC_ID_sing = RAW->get_FAT_TDC_id(i);         //come back to this



                    if(dummy_qdc_E[i] == 0) cout << "i = i zero energy" << endl;

                    dum_qdc_id[matchedmult] =  dummy_qdc_id[i];
                    dum_qdc_E[matchedmult] = dummy_qdc_E[i];
                    dum_qdc_E_raw[matchedmult] = dummy_qdc_E_raw[i];
                    dum_qdc_t_coarse[matchedmult] = dummy_qdc_t_coarse[i];
                    dum_qdc_t_fine[matchedmult] = dummy_qdc_t_fine[i];
                    
                    dum_tdc_id[matchedmult] =  dummy_tdc_id[i];
                    dum_tdc_t[matchedmult] = dummy_tdc_t[i];
                    dum_tdc_t_raw[matchedmult] = dummy_tdc_t_raw[i];

                    matchedmult++;
                    
                    //cout << "matched mult: " << matchedmult << " i = i, i: " << i << " i " << i << " QDCID: " << dummy_qdc_id[matchedmult] << " TDCID " << dummy_tdc_id[matchedmult] << " energy " << dummy_qdc_E[matchedmult] << " timestamp " << dummy_tdc_t[matchedmult] << endl;
                    
                    if( dummy_qdc_id[i] != dummy_tdc_id[i]) cout << "error --------------------------------------------------" << endl;

                    if(dummy_qdc_id[i] == -1 && dummy_tdc_id[i] == -1){
                    
                        cout << "there was a -1 match" << endl;
                        matchedmult--;
                    }

                    }//end if check           
              
            
                    else {
                    for (int j=0; j < dummytdcmult; j++){
                        
                    if(dummy_qdc_id[i] == dummy_tdc_id[j] && EnID[i] == 0 && TimID[j] == 0){
                        //cout << "qdc id " << dummy_qdc_id[i] << " tdc id " << dummy_tdc_id[j] << endl;

                    EnID[i] = true;
                    TimID[j] = true;
                    
                    //Fat_TDC_ID_sing = RAW->get_FAT_TDC_id(j);

                    //Fat_TDC_multi[Fat_TDC_ID_sing]++;


                    dum_qdc_id[matchedmult] =  dummy_qdc_id[i];
                    dum_qdc_E[matchedmult] = dummy_qdc_E[i];
                    dum_qdc_E_raw[matchedmult] = dummy_qdc_E_raw[i];
                    dum_qdc_t_coarse[matchedmult] = dummy_qdc_t_coarse[i];
                    dum_qdc_t_fine[matchedmult] = dummy_qdc_t_fine[i];
                    
                    dum_tdc_id[matchedmult] =  dummy_tdc_id[j];
                    dum_tdc_t[matchedmult] = dummy_tdc_t[j];
                    dum_tdc_t_raw[matchedmult] = dummy_tdc_t_raw[j];
                    
                    
                    if(dummy_qdc_E[i] == 0) cout << "i = j zero energy" << endl;


                    matchedmult++;
                    
                    if(dummy_qdc_id[j] == -1 && dummy_tdc_id[j] == -1){
                    
                        cout << "there was a -1 match" << endl;
                        matchedmult--;
                    }
                    //cout << "matched mult: " << matchedmult << "i = j, i: " << i << " j " << j << " QDCID: " << dummy_qdc_id[i] << " TDCID " << dummy_tdc_id[j] << " energy " << dummy_qdc_E[i] << " timestamp " << dummy_tdc_t[j] << endl;
                    if( dummy_qdc_id[i] != dummy_tdc_id[j]) cout << "error --------------------------------------------------" << endl;
                        } // end if 
                    
                    
                    }//end j    
                    
                    
                }//end else
 
            
            //if(counter > 1)cout << "i: " << i << " j " << i << " TDCID: " << RAW->get_FAT_TDC_id(i) << " QDCID " << RAW->get_FAT_QDC_id(i) << " energy " << RAW->get_FAT_QLong(i) << " timestamp" << RAW->get_FAT_TDC_timestamp(i) << endl;

            }//end of i loop

            //cout << "mult in unpacker: " << matchedmult << endl;
    
                        
          //fOutput->fFat_firedQDC = matchedqdcmult;
          //fOutput->fFat_firedTDC = matchedtdcmult;
          
          fOutput->fFat_mult = matchedmult;
        
          fOutput->fFat_SC40mult =  sc40count;
          fOutput->fFat_SC41mult =  sc41count;
          
          fOutput->fFat_bPlasmult = fatbplascount;
          
          fOutput->fFat_tdcsinglescount =  singlestdcmult;
          fOutput->fFat_qdcsinglescount =  singlesqdcmult;

        for(int i = 0; i < matchedmult; i++){
            
                fOutput->fFat_QDC_ID[i] =  dum_qdc_id[i];
                fOutput->fFat_QDC_E[i] = dum_qdc_E[i];
                fOutput->fFat_QDC_E_Raw[i] = dum_qdc_E_raw[i];
                fOutput->fFat_QDC_T_coarse[i] = dum_qdc_t_coarse[i];
                fOutput->fFat_QDC_T_fine[i] = dum_qdc_t_fine[i];

                fOutput->fFat_TDC_ID[i] =  dum_tdc_id[i];
                fOutput->fFat_TDC_Time[i] = dum_tdc_t[i];
                fOutput->fFat_TDC_Time_Raw[i] = dum_tdc_t_raw[i];

                //cout << "event no: " << event_number << " mult = " << matchedmult <<  " i: " << i << " QDCID: " << dum_qdc_id[i] << " TDCID " << dum_tdc_id[i] << " energy " << dum_qdc_E[i] << " timestamp " << dum_tdc_t[i] << endl;

                if(dum_tdc_t[i] == 0.)cerr << "Something is wrong in Unpack Fatima pointer out (dum_tdc_t==0)!!" << endl;
            
        }
     
    }//end if used systems
        ///--------------------------------------------------------------------------------------------///
                                            /**Output GALILEO **/
        ///--------------------------------------------------------------------------------------------///
        if (Used_Systems[5]&& PrcID_Conv==5){
         for (int i=fOutput->fGal_fired; i<RAW->get_GALILEO_am_Fired() && i < GALILEO_MAX_HITS; i++){
                fOutput->fGal_Detector[i] =  RAW->get_GALILEO_Det_id(i);
                fOutput->fGal_Crystal[i] =  RAW->get_GALILEO_Crystal_id(i);
                fOutput->fGal_Event_T[i] = RAW->get_GALILEO_Event_T(i);
                fOutput->fGal_E[i] = RAW->get_GALILEO_Chan_E(i);
                fOutput->fGal_T[i] = RAW->get_GALILEO_Chan_T(i);
                fOutput->fGal_Pileup[i] = RAW->get_GALILEO_Pileup(i);
                fOutput->fGal_Overflow[i] = RAW->get_GALILEO_Overflow(i);
                fOutput->fGal_fired++;
   
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                        /** Output FINGER **/
  ///--------------------------------------------------------------------------------------------///
      if (Used_Systems[6]&& PrcID_Conv==6){

          int Phys_Channel_Lead[FINGER_TAMEX_MODULES][FINGER_TAMEX_HITS] = {0,0};
          int Phys_Channel_Trail[FINGER_TAMEX_MODULES][FINGER_TAMEX_HITS] = {0,0};

          int fingfired[FINGER_TAMEX_MODULES] = {0};

          for (int i=0; i<RAW->get_FINGER_tamex_hits(); i++){
            fingfired[i] = RAW->get_FINGER_am_Fired(i);
     
            for(int j = 0;j < fingfired[i];j++){
        
              if(RAW->get_FINGER_CH_ID(i,j) % 2 == 1){ //Lead odd j
                Phys_Channel_Lead[i][j] = fingID[i][RAW->get_FINGER_physical_channel(i, j)]; //From allocation file
                int chan = Phys_Channel_Lead[i][j];

                if (chan < 0)
                  continue;

                // PMT allocation succeeded
                int N1 = fOutput->fFing_PMT_Lead_N[chan]++;
                fOutput->fFing_Lead_PMT[chan][N1] = RAW->get_FINGER_lead_T(i,j);
  
                // PMT "0" is the trigger
                if (chan == 0 || chan == 1){
                    fOutput->fFing_SC41_lead[chan][N1] = RAW->get_FINGER_lead_T(i,j);
                  continue;
                }
                // chan = "PMT" number
                // this maps to two strips to fill in
                if (chan % 2 == 0) // even PMT = up pmts
                {
                  int strip1 = chan;
                  int strip2 = chan + 1;
                  int N1 = fOutput->fFing_Strip_N_LU[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_LU[strip2]++;
                  fOutput->fFing_Lead_Up[strip1][N1] = RAW->get_FINGER_lead_T(i,j);
                  fOutput->fFing_Lead_Up[strip2][N2] = RAW->get_FINGER_lead_T(i,j);
                  fOutput->fFing_Strip_N[strip1]++;
                  fOutput->fFing_Strip_N[strip2]++;
                                }
                else // odd = lower PMT
                {
                  int strip1 = chan + 1;
                  int strip2 = chan;
                  int N1 = fOutput->fFing_Strip_N_LD[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_LD[strip2]++;
                  fOutput->fFing_Lead_Down[strip1][N1] = RAW->get_FINGER_lead_T(i,j);
                  fOutput->fFing_Lead_Down[strip2][N2] = RAW->get_FINGER_lead_T(i,j);
                      }
              }
              else{ //Trail even j
                Phys_Channel_Trail[i][j] = fingID[i][RAW->get_FINGER_physical_channel(i,j)];
             
                int chan = Phys_Channel_Trail[i][j];
                if (chan < 0)
                  continue;

                // PMT allocation succeeded
                int N1 = fOutput->fFing_PMT_Trail_N[chan]++;
                fOutput->fFing_Trail_PMT[chan][N1] = RAW->get_FINGER_trail_T(i,j);
                 // PMT "0" is the trigger
                if (chan == 0 || chan == 1){
                    fOutput->fFing_SC41_trail[chan][N1] = RAW->get_FINGER_trail_T(i,j);

                  continue;
                }
                if (chan % 2 == 0) // even PMT = up pmts
                {
                  int strip1 = chan + 1;
                  int strip2 = chan;
                  int N1 = fOutput->fFing_Strip_N_TU[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_TU[strip2]++;
                  fOutput->fFing_Trail_Up[strip1][N1] = RAW->get_FINGER_trail_T(i,j);
                  fOutput->fFing_Trail_Up[strip2][N2] = RAW->get_FINGER_trail_T(i,j);
                 }
                else // odd = lower PMT
                {
                  int strip1 = chan + 1;
                  int strip2 = chan;
                  int N1 = fOutput->fFing_Strip_N_TD[strip1]++;
                  int N2 = fOutput->fFing_Strip_N_TD[strip2]++;
                  fOutput->fFing_Trail_Down[strip1][N1] = RAW->get_FINGER_trail_T(i,j);
                  fOutput->fFing_Trail_Down[strip2][N2] = RAW->get_FINGER_trail_T(i,j);
                }
              }
            }
          }
        }      
        ///--------------------------------------------------------------------------------------------///

      } //End of subevent loop


      fOutput->SetValid(isValid);

      pdata = nullptr;
    //} //End of Skip
  }
  //
  return isValid;

}

void EventUnpackProc::FILL_HISTOGRAMS(int PrcID_Conv, int PrcID, int SubType){
    
 // switch(PrcID_Conv){
  //  case 0:
    ///WARNING CHECK IF USED SYSTEMS WORKS!!!
   
  if(PrcID_Conv==0 && Used_Systems[0])  Fill_FRS_Histos(PrcID,Type,SubType);

  if(PrcID_Conv==1 && Used_Systems[1])  Fill_AIDA_Histos();
 
  if(PrcID_Conv==3 && Used_Systems[3]) Fill_FATIMA_VME_Histos();

  //if(VME_TAMEX_Fatima==false && VME_AND_TAMEX_Fatima==false && PrcID_Conv==4)  Fill_FATIMA_TAMEX_Histos();
  
  if(VME_AND_TAMEX_Fatima==true && (PrcID_Conv==4)) Fill_FATIMA_TAMEX_Histos();

  if(PrcID_Conv==5 && Used_Systems[5]) Fill_GALILEO_Histos();
}


//-----------------------------------------------------------------------------------------------------------------------------//
void EventUnpackProc::ResetMultiplexer()
{
  for (int i = 0; i < 12; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      adcLastTimestamp[i][j] = 0;
      adcCounts[i][j] = 0;
    }
  }
}



void EventUnpackProc::CorrectTimeForMultiplexer(AidaEvent &evt)
{
  int fee = evt.Module;
  int adc = evt.Channel / 16;
  int64_t time = evt.Time;

  if ((time - adcLastTimestamp[fee][adc] > 2500) && adcLastTimestamp[fee][adc] != 0)
  adcCounts[fee][adc] = 0;

  adcLastTimestamp[fee][adc] = time;

  evt.Time = time - (2000 * adcCounts[fee][adc]++);
  if (evt.HighEnergy) evt.FastTime = evt.Time;
}

//-----------------------------------------------------------------------------------------------------------------------------//


void EventUnpackProc::load_PrcID_File(){
  ifstream data("Configuration_Files/DESPEC_General_Setup/PrcID_to_Det_Sys.txt");
  if(data.fail()){
    cerr << "Could not find PrcID config file!" << endl;
    exit(0);
  }
  int id[9] = {0,0,0,0,0,0,0,0,0};
  int i = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d %d %d %d %d %d %d %d %d",s_tmp,&id[0],&id[1],&id[2],&id[3],&id[4],&id[5],&id[6],&id[7],&id[8]);
    for(int j = 0; j < 8; ++j){ PrcID_Array[i][j] = id[j];
  
    }
    i++;
   
  }
}
//---------------------------------------------------------------------------------------------------
void EventUnpackProc::load_FatTamex_Allocationfile(){

  const char* format = "%d %d %d";
  ifstream data("Configuration_Files/FATIMA/FATIMA_TAMEX_allocation.txt");
  if(data.fail()){
    cerr << "Could not find Fatima_TAMEX_allocation config file!" << endl;
    exit(0);
  }
  //     int id[5] = {0,0,0,0,0};
  //int i = 0;
   
  for(int i=0; i<5; i++){
      for(int j=0; j<80; j++){
   TAMEX_bPlasFat_ID[i][j]=0;
      }
  }
  int TamID = 0;
  int TamCh = 0;
  int Sys_ch =0;
  string line;
  //char s_tmp[100];
  while(data.good()){

    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&TamID,&TamCh,&Sys_ch);
    TAMEX_bPlasFat_ID[TamID][TamCh] = Sys_ch;
  }
}
//---------------------------------------------------------------------------------------------------
void EventUnpackProc::load_FingerID_File(){

  const char* format = "%d %d %d";
  ifstream data("Configuration_Files/Finger/Finger_allocation.txt");
  if(data.fail()){
    cerr << "Could not find Finger_allocation config file!" << endl;
    exit(0);
  }
  //     int id[5] = {0,0,0,0,0};
  //int i = 0;
  int tamid = 0;
  int tamch = 0;
  int fingid = 0;
  string line;
  //char s_tmp[100];
  while(data.good()){

    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&tamid,&tamch,&fingid);
    fingID[tamid][tamch] = fingid;
  }
}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::read_setup_parameters(){

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
  file >> WHITE_RABBIT_USED;//dummy_var;


  cout<<endl;
  cout<<endl;
  cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Setup Parameters List Unpack Proc: "<<endl;
  if(WHITE_RABBIT_USED) cout<<"White Rabbit: Enabled"<<endl;
  else if(!WHITE_RABBIT_USED) cout<<"White Rabbit: Disabled"<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
  cout<<endl;
  cout<<endl;
 
}

//-----------------------------------------------------------------------------------------------------------------------------//
Int_t EventUnpackProc::get_Conversion(Int_t PrcID){

  for(int i = 0;i < 7;++i){
    for(int j = 0;j < 8;++j){
        ///Fix for FRS 
          if (PrcID==100) {return -1; }
      if(PrcID == PrcID_Array[i][j]) {
          
          return i;
      }
    }
  }
  cerr << "ProcID " << PrcID << " not known!" << endl;
  exit(0);
}

void EventUnpackProc::get_used_systems(){
    for(int i = 0;i < 6;i++) Used_Systems[i] = false;

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

    cout << "\n=====================================================" << endl;
    cout << "USED SYSTEMS" << endl;
    cout << "-----------------------------------------------------" << endl;
    for(int j = 0;j < 6;++j){
        if(Used_Systems[j]) cout << DET_NAME[j] << endl;
    }
    cout << "=====================================================" << endl;


}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::get_WR_Config(){
  ifstream data("Configuration_Files/DESPEC_General_Setup/White_Rabbit.txt");
  if(data.fail()){
    cerr << "Could not find White_Rabbit config file!" << endl;
    exit(0);
  }

  int id = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d",s_tmp,&id);
    WR_used = (id == 1);
  }
}

  //-----------------------------------------------------------------------------------------------------------------------------//
  // ################################################################## //
  // ################################################################## //
  // ################# Raw Histogram Filling Section ################## //
  // ################################################################## //
  // ################################################################## //
  
  
  
  /**----------------------------------------------------------------------------------------------**/
  /**---------------------------------------------  FRS  ------------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/

  void EventUnpackProc::Make_FRS_Histos(){
 char fname[50], name[50], title[60];//, title2[60];

  const char *count_title1[12]={"(0:1)", "(1:1)", "(2:1)",
                "(2:2)", "(3:1)", "(4:1)",
                "(4:2)", "(4:3)", "(6:1)",
                "(6:2)", "(8:1)"};
  const char *fext1[12]={"0", "1", "2", "2", "3", "4", "4", "4", "6", "6", "8", "8"};
  const char *fext2[12]={"01", "11", "21", "22","31", "41",
             "42", "43", "61",
             "62", "81", "82"};
             
  ///FRS Scalars           
  bool scaler_enable_hist[64];
  char scaler_name[64][256];
  scaler_ch_1kHz=39; //ch7 of 2nd scaler
  scaler_ch_spillstart=8; //ch8 of 1st scaler 
  scaler_check_first_event=1;      
  for(int ii=0; ii<64; ii++){
    sprintf(scaler_name[ii],"scaler_ch%d",ii);//default name
    scaler_enable_hist[ii]=false;
  }
  sprintf(scaler_name[0],"IC01curr-old"); 
  sprintf(scaler_name[1],"SEETRAM-old");
  sprintf(scaler_name[2],"SEETRAM-new");
  sprintf(scaler_name[3],"IC01curr-new");
  sprintf(scaler_name[4],"IC01 count");
  sprintf(scaler_name[5],"SCI00");
  sprintf(scaler_name[6],"SCI01");
  sprintf(scaler_name[7],"SCI02");
  sprintf(scaler_name[8],"Start Extr");
  sprintf(scaler_name[9],"Stop Extr");
  sprintf(scaler_name[10],"Beam Transformer");
  
  sprintf(scaler_name[32],"Free Trigger");
  sprintf(scaler_name[33],"Accept Trigger");
  sprintf(scaler_name[34],"Spill Counter");
  sprintf(scaler_name[35],"1 Hz clock");
  sprintf(scaler_name[36],"10 Hz clock");
  sprintf(scaler_name[37],"100 kHz X veto dead-time");
  sprintf(scaler_name[38],"100 kHz clock");
  sprintf(scaler_name[39],"1 kHz clock");
  
  sprintf(scaler_name[48],"SCI21L");
  sprintf(scaler_name[49],"SCI41L");
  sprintf(scaler_name[50],"SCI42L");
  sprintf(scaler_name[51],"SCI43L");
  sprintf(scaler_name[52],"SCI81L");
  sprintf(scaler_name[53],"SCI21R");
  sprintf(scaler_name[54],"SCI41R");
  sprintf(scaler_name[55],"SCI42R");
  sprintf(scaler_name[56],"SCI43R");
  sprintf(scaler_name[57],"SCI81R");
  sprintf(scaler_name[58],"SCI31L");
  sprintf(scaler_name[59],"SCI31R");
  sprintf(scaler_name[60],"SCI11");
  sprintf(scaler_name[61],"SCI51");
  

  for(int ii=0; ii<64; ii++){
    hScaler_per_s[ii]     = MakeH1I("FRS/Scaler/Scaler_per_1s",Form("%s_per_1s",scaler_name[ii]),1000,0,1000,"Time (s)", 2,5, "Count per second");
    hScaler_per_100ms[ii] = MakeH1I("FRS/Scaler/Scaler_per_0.1s",Form("%s_per_0.1s",scaler_name[ii]),4000,0,400,"Time (s)", 2,5, "Count per 0.1 second");
    hScaler_per_spill[ii] = MakeH1I("FRS/Scaler/Scaler_per_spill",Form("%s_per_spill",scaler_name[ii]),1000,0,1000,"Spill", 2,5, "Count per spill");
  }

  for (int cnt = 0; cnt<7; cnt++) //changed from 3 to 6 04.07.2018
    {
      int index = 0;
      switch(cnt)
    {
        case 0: index = 2; break;
        case 1: index = 3; break;
        case 2: index = 5; break;
        case 3: index = 6; break;
        case 4: index = 7; break;
        case 5: index = 10; break;
        case 6: index = 8; break;
    }
      sprintf(fname,"FRS/SCI/SCI%s/SCI%s",fext1[index],fext2[index]);
      sprintf(name, "SCI%s_L", count_title1[index]);
      sprintf(title, "Sc%s L dE [ch]", count_title1[index]);
      hSCI_L[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_R", count_title1[index]);
      sprintf(title, "Sc%s R dE [ch]", count_title1[index]);
      hSCI_R[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_E", count_title1[index]);
      sprintf(title, "Sc%s Energy [ch]", count_title1[index]);
      hSCI_E[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_Tx", count_title1[index]);
      sprintf(title, "Sc%s t_lr [ch] TAC", count_title1[index]);
      hSCI_Tx[index] = MakeH1I(fname,name,4096,0,4096,title,2,3);

      sprintf(name, "SCI%s_X", count_title1[index]);
      sprintf(title, "Sc%s x-pos [mm]", count_title1[index]);
      hSCI_X[index] = MakeH1I(fname,name,240,-120,120,title,2,3);

    }
      hSCI_dE24 = MakeH2I("FRS/SCI","SCI_dE21-41", 100,10,4000,100,10,4000,"SC21 dE","SC41 dE",2);
     // ToF SC21-SC41
        sprintf(fname,"FRS/SCI/TOF/TOF(%d)",2);
        sprintf(name,"SCI_21_41_TofLL");
        hSCI_TofLL2 = MakeH1I(fname,name,1500,0,100000,"TAC SC41L-SC21L [ps]",2,3);

        sprintf(name,"SCI_21_41_TofRR");
        hSCI_TofRR2 = MakeH1I(fname,name,1500,0,100000,"TAC SC41R-SC21R [ps]",2,3);

//         hSCIdE41_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE41_TPC42X", 1024,0,4096, 400,-100.,100, "SC41 dE", "TPC42 X[mm]", 2);
//         hSCIdE41L_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE41L_TPC42X", 1024,0,4096, 400,-100.,100, "SC41L dE", "TPC42 X[mm]", 2);
//         hSCIdE41L_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE41L_TPC41X", 1024,0,4096, 400,-100.,100, "SC41L dE", "TPC41 X[mm]", 2);
//         hSCIdE41R_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE41R_TPC42X", 1024,0,4096, 400,-100.,100, "SC41R dE", "TPC42 X[mm]", 2);
//         hSCIdE41R_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE41R_TPC41X", 1024,0,4096, 400,-100.,100, "SC41R dE", "TPC41 X[mm]", 2);
// 
//         hSCIdE21_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE21_TPC42X", 1024,0,4096, 400,-100.,100, "SC41 dE", "TPC42 X[mm]", 2);
//         hSCIdE21L_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE21L_TPC42X", 1024,0,4096, 400,-100.,100, "SC21L dE", "TPC42 X[mm]", 2);
//         hSCIdE21L_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE21L_TPC41X", 1024,0,4096, 400,-100.,100, "SC21L dE", "TPC41 X[mm]", 2);
//         hSCIdE21R_TPC42X= MakeH2I("FRS/SCI_TPC/","SCIdE21R_TPC42X", 1024,0,4096, 400,-100.,100, "SC21R dE", "TPC42 X[mm]", 2);
//         hSCIdE21R_TPC41X= MakeH2I("FRS/SCI_TPC/","SCIdE21R_TPC41X", 1024,0,4096, 400,-100.,100, "SC21R dE", "TPC41 X[mm]", 2);

        // ToF SC21-SC42 changed on 03.07.2018 SB
        sprintf(fname,"FRS/SCI/TOF/TOF(%d)",3);
        sprintf(name,"SCI_21_42_TofLL");
        hSCI_TofLL3 = MakeH1I(fname,name,1500,0,100000,"TAC SC42L-SC21L [ps]",2,3);

        sprintf(name,"SCI_21_42_TofRR");
        hSCI_TofRR3 = MakeH1I(fname,name,1500,0,100000,"TAC SC42R-SC21R [ps]",2,3);

        sprintf(name,"SCI_21_42_Tof3");
        hSCI_Tof3 = MakeH1I(fname,name,1000,0,100000,"TAC SC42-SC21 [ps] (pos.corr.)",2,3);
        
        
        

        sprintf(fname,"FRS/SCI/TOF/TOF(%d)",5);
        sprintf(name,"SCI_22_41_TofLL");
        hSCI_TofLL5 = MakeH1I(fname,name,1500,0,100000,"TAC SC41L-SC22L [ps]",2,3);

        sprintf(name,"SCI_22_41_TofRR");
        hSCI_TofRR5 = MakeH1I(fname,name,1500,0,100000,"TAC SC41R-SC22R [ps]",2,3);

        sprintf(name,"SCI_22_41_Tof5");
        hSCI_Tof5 = MakeH1I(fname,name,1000,0,100000,"TAC SC41-SC22 [ps] (pos.corr.)",2,3);
        
        
        

    hSCI_dT_21l_41l = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21l_41l","hSCI_dT_21l_41l",5001,0,5000); //from Multihit TDCS
    hSCI_dT_21r_41r = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21r_41r","hSCI_dT_21r_41r",5001,0,5000);

    hSCI_dT_21l_42l = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21l_42l","hSCI_dT_21l_42l",5001,0,5000);
    hSCI_dT_21r_42r = MakeTH1('D',"FRS/SCI/dT/SCI_dt_21r_42r","hSCI_dT_21r_42r",5001,0,5000);

    //ID
    hID_AoQ = MakeH1I("FRS/ID","ID_AoQ",2000,1.4,5.0,"A/Q S2-S4",2,6);
    hID_AoQ_corr = MakeH1I("FRS/ID","ID_AoQ_corr",2000,1.4,3.0,"A/Q S2-S4",2,6);
  //   hID_Z = MakeH1I("ID",Form("ID_Z, gain=%f",music->e1_gain[0]),1000,10,93,"Z s2-s4",2,6);
    hID_Z = MakeH1I("FRS/ID","ID_Z",1000,0,93,"Z s2-s4",2,6);
    hID_Z2 = MakeH1I("FRS/ID","ID_Z2",1000,0,93,"Z2 s2-s4",2,6);
   // hID_Z3 = MakeH1I("FRS/ID","ID_Z3",1000,10,93,"Z3 s2-s4",2,6);
    ////////////////////////////////////////////////////////////

    hID_Z_dE2 = MakeH2I("FRS/ID","ID_Z_dE2", 250,1,30, 250,0.,4000.,
              "Z", "MUSIC2_dE", 2);

    hID_Z_Sc21E = MakeH2I("FRS/ID","ID_Z_Sc21E", 300,0,25.,400,0,4000.,
            "Z s2-s4", "sqrt(Sc21_L*sC21_R)", 2);
    hID_x2z = MakeH2I("FRS/ID","ID_x2z", 1500,30.,60., 300,-150.,150., "Z s2-s4", "X at S2 [mm]", 2);
    hID_x4z = MakeH2I("FRS/ID","ID_x4z", 1500,30.,60., 300,-150.,150., "Z s2-s4", "X at S4 [mm]", 2);
    hID_E_Xs4 = MakeH2I("FRS/ID","ID_E_Xs4", 200,-100.,100., 400,0.,4000., "X s4 [mm]", "Delta E", 2);
    hID_E_Xs2 = MakeH2I("FRS/ID","ID_E_Xs2", 200,-100.,100., 400,0.,4000., "X s2 [mm]", "Delta E", 2);
    hID_x2a2 = MakeH2I("FRS/ID", "ID_x2_a2", 200, -100., 100., 200, -100., 100., "X s2 [mm]", "AngleX s2 [mrad]", 2);
    hID_y2b2 = MakeH2I("FRS/ID", "ID_y2_b2", 200, -100., 100., 200, -100., 100., "Y s2 [mm]", "AngleY s2 [mrad]", 2);
    hID_x4a4 = MakeH2I("FRS/ID", "ID_x4_a4", 200, -100., 100., 200, -100., 100., "X s4 [mm]", "AngleX s4 [mrad]", 2);
    hID_y4b4 = MakeH2I("FRS/ID", "ID_y4_b4", 200, -100., 100., 200, -100., 100., "Y s4 [mm]", "AngleY s4 [mrad]", 2);
    hID_x2x4 = MakeH2I("FRS/ID","ID_x2_x4",200,-100,100,200,-100,100,"x2 mm","x4 mm",2);
    hID_SC41dE_AoQ = MakeH2I("FRS/ID","ID_SC41dE_AoQ", 300,1.2,3.0, 800,0.,4000.,"A/Q s2-s4", "SC41 dE", 2);
   
    hID_TACvsMHTDC_AoQ= MakeH2I("FRS/ID","ID_TACvsMHTDC_AoQ", 1500,1.9,2.7,  1500,1.3,2.7,"A/Q MHTDC", "AoQ TAC", 2);
    
     hID_TACvsMHTDC_Z= MakeH2I("FRS/ID","ID_TACvsMHTDC_Z", 2000,10.,60.,  2000,10.,60.,"Z MHTDC", "AoQ TAC", 2);
     
     hID_TACvsMHTDC_Beta= MakeH2I("FRS/ID","ID_TACvsMHTDC_Beta", 2500,0.,1.,  2500,0.,1.,"Beta MHTDC", "Beta TAC", 2);
   
     hID_ZTACvsAoQMHTDC= MakeH2I("FRS/ID","ID_ZTACvsAoQMHTDC", 1500,1.9,2.7,  1500,10.,60.,"A/Q MHTDC", "Z TAC", 2);
     
     hID_ZMHTDCvsAoQTAC= MakeH2I("FRS/ID","ID_ZMHTDCvsAoQTAC", 1500,1.9,2.7,  1500,10.,60.,"A/Q TAC", "Z MHTDC", 2);
   
    ///////////////////////////////////////////////////////////////////////////////////////////////
for(int i=0;i<7;i++)
    {
      char fname[100];
      char name[100];
      sprintf(fname,"FRS/TPC/%s/",tpc_folder_ext1[i]);

      hTPC_X[i]=MakeH1I_TPC(fname,"X",i,800,-100.,100.,
                "x[mm]",2,3);
      hTPC_Y[i]=MakeH1I_TPC(fname,"Y",i,800,-100.,100.,
                "y[mm]",2,3);


      sprintf(name,"%s%s",tpc_name_ext1[i],"XY");
      hcTPC_XY[i]=MakeH2I(fname,name, 120,-120.,120., 120,-120.,120.,
              "X [mm] ","Y [mm] ", 2);

      sprintf(name,"%s%s",tpc_name_ext1[i],"LTRT");
      hTPC_LTRT[i]=MakeH2I(fname,name, 2048,0,4095, 2048,0,4095,
               "LT [ch]","RT[ch] ", 2);
      hTPC_DELTAX[i]=MakeH1I_TPC(fname,"x0-x1",i,100,-10.,10.,
                 "x0-x1[mm]",2,3);

    }
    hID_x2 = MakeTH1('D',"FRS/TPC/S2_X","ID_x2",200,-100,100);
    hID_y2 = MakeTH1('D',"FRS/TPC/S2_Y","ID_y2",200,-100,100);
    hID_a2 = MakeTH1('D',"FRS/TPC/S2_angA","ID_a2",200,-100,100);
    hID_b2 = MakeTH1('D',"FRS/TPC/S2_angB","ID_b2",200,-100,100);

    hID_x4 = MakeTH1('D',"FRS/TPC/S4_X","ID_x4",800,-100,100);
    hID_y4 = MakeTH1('D',"FRS/TPC/S4_Y","ID_y4",200,-100,100);
    hID_a4 = MakeTH1('D',"FRS/TPC/S4_angA","ID_a4",200,-100,100);
    hID_b4 = MakeTH1('D',"FRS/TPC/S4_angA","ID_b4",200,-100,100);

    htpc_X2 = MakeTH1('D',"FRS/TPC/S2_TPCX","tpc_x S21",800,-100,100);
    htpc_Y2 = MakeTH1('D',"FRS/TPC/S2_TPCY","tpc_y S21",800,-100,100);
    htpc_X4 = MakeTH1('D',"FRS/TPC/S4_TPCX","tpc_x S41",800,-100,100);
    htpc_Y4 = MakeTH1('D',"FRS/TPC/S4_TPCY","tpc_y S41",800,-100,100);

    hID_beta = MakeH1I("FRS/ID","ID_beta",1000,0,1000,"id.beta(2)*1000",2,6);

    hID_dEToF = MakeH2I("FRS/ID","ID_dEToF", 2000, 00000.,70000.,400,0,4000, "tof S2-S4 Sci.Tof(2)", "Music_dE(1)", 2);
    hID_BRho[0] = MakeH1I("FRS/ID","ID_BRho0",5000,2.5,14.5,"BRho of 1. Stage [Tm]",2,6);
    hID_BRho[1] = MakeH1I("FRS/ID","ID_BRho1",5000,2.5,14.5,"BRho of 2. Stage [Tm]",2,6);

  // char name[80], xtitle[80];
   for(int i=0;i<8;i++)
     {
       hMUSIC1_E[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(1)/Energy/EnergyM1%2d",i), Form("Music 1 E%2d",i), 4096,0,4096);
       hMUSIC2_E[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(2)/Energy/EnergyM2%2d",i), Form("Music 2 E%2d",i), 4096,0,4096);
       hMUSIC1_T[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(1)/Time/TimeM1%2d",i), Form("Music 1 T%2d",i), 4096,0,4096);
       hMUSIC2_T[i] = MakeTH1('D', Form("FRS/MUSIC/MUSIC(2)/Time/TimeM2%2d",i), Form("Music 2 T%2d",i), 4096,0,4096);
     }
 
  //  hMUSIC1_dE1dE2 = MakeTH2('D',"FRS/MUSIC/MUSIC(1)/E1E2","E1_E2", 1024,0,4096,1024,0,4096);
   
  //  hMUSIC1_MUSIC2 = MakeTH2('D',"FRS/MUSIC/MUSIC1_MUSIC2","dE1_dE2", 1024,0,4096,1024,0,4096);  
   
    for(int i=0;i<32;i++)
     {
       hvftx_TRaw[i] = MakeTH1('D', Form("FRS/VFTXSCI/TRaw/TRaw%2d",i), Form("T Raw %2d",i), 2000,1E6,1E8);
     }
   
   hvftx_Sci21PosRaw = MakeTH1('D',"FRS/VFTXSCI/PosRaw21","Sci21_PosRaw",1000,-20000.,20000,"PosRaw = Left - Right [10ps/bin]");
   hvftx_Sci22PosRaw = MakeTH1('D',"FRS/VFTXSCI/PosRaw22","Sci22_PosRaw",1000,-20000.,20000,"PosRaw = Left - Right [10ps/bin]");
  
   //hvftx_Sci21PosRaw_TRAW = MakeTH1('D',"FRS/VFTXSCI/PosRaw21_TRAW","Sci21_PosRaw_TRAW",1000,-20000.,20000,"PosRaw = Left - Right [10ps/bin]");
  
   hvftx_Sci41PosRaw = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/PosRaw41","Sci41_PosRaw",1000,-20000.,20000,"PosRaw = Left - Right [10ps/bin]"); 
   hvftx_Sci42PosRaw = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/PosRaw42","Sci42_PosRaw",1000,-20000.,20000,"PosRaw = Left - Right [10ps/bin]"); 
  
   hvftx_ToFraw_2141 = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/ToFraw_2141","ToFraw_41-21_mult1",1000,-20000.,20000,"ToFraw [10ps] ");
   
   hvftx_ToFraw_2141LL = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/ToFraw_2141LL","ToFraw_41-21_mult1_LL",1000,-20000.,20000,"ToFraw [10ps] ");
   
   hvftx_ToFraw_2141RR = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/ToFraw_2141RR","ToFraw_41-21_mult1_RR",1000,-20000.,20000,"ToFraw [10ps] ");
   
   hvftx_ToFraw_2142LL = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/ToFraw_2142LL","ToFraw_42-21_mult1_LL",1000,-20000.,20000,"ToFraw [10ps] ");
   
   hvftx_ToFraw_2142RR = (TH1D*)MakeTH1('D',"FRS/VFTXSCI/ToFraw_2142RR","ToFraw_42-21_mult1_RR",1000,-20000.,20000,"ToFraw [10ps] ");
  

    htimestamp = MakeTH1('D',"FRS/timestamp","timestamp",30,0.,300.);
    hts = MakeTH1('D',"FRS/ts","ts",30,0.,300.);
    hts2 = MakeTH1('D',"FRS/ts2","ts2",30,0.,300.);

    
  }
  //-----------------------------------------------------------------------------------------------------------------------------//
  void EventUnpackProc::Fill_FRS_Histos(int PrcID, int Type, int SubType){

     time_in_ms = 0;
     spill_count = 0;
     ibin_for_s = 0;
     ibin_for_100ms = 0;
     ibin_for_spill = 0;


    for(int i =0; i<8; i++){
        Music_E1[i] = 0;
        Music_E2[i] = 0;
        Music_T1[i] = 0;
        Music_T2[i] = 0;
    }
    
    TRaw_vftx_21l =0;
    TRaw_vftx_21r =0;
    TRaw_vftx_22l =0;
    TRaw_vftx_22r =0;
    TRaw_vftx_41l =0;
    TRaw_vftx_41r =0;
    TRaw_vftx_42l =0;
    TRaw_vftx_42r =0;
   
 
    for(int i =0; i<3; i++){
        Music_dE[i] = RAW->get_FRS_MusicdE(i);
        Music_dE_corr[i] = RAW->get_FRS_MusicdE_corr(i);
       
    }
    for(int i=0; i<8; i++){
        Music_E1[i] = RAW->get_FRS_MusicE1(i);
        Music_E2[i] = RAW->get_FRS_MusicE2(i);
        Music_T1[i] = RAW->get_FRS_MusicT1(i);
        Music_T2[i] = RAW->get_FRS_MusicT2(i);
 
    }

    for(int l=0;l<12;++l){
        sci_l[l] = RAW->get_FRS_sci_l(l);
        sci_r[l] = RAW->get_FRS_sci_r(l);
        sci_e[l] = RAW->get_FRS_sci_e(l);
        sci_tx[l] = RAW->get_FRS_sci_tx(l);
        sci_x[l] = RAW->get_FRS_sci_x(l);
        
      
    }
    sci_tofll2 = RAW->get_FRS_tofll2();
    sci_tofll3 = RAW->get_FRS_tofll3();
    sci_tof2 = RAW->get_FRS_tof2();
    sci_tofrr2 = RAW->get_FRS_tofrr2();
    sci_tofrr3 = RAW->get_FRS_tofrr3();
    sci_tof3 = RAW->get_FRS_tof3();
    sci_tofll5 = RAW->get_FRS_tofll5();
    sci_tofrr5 = RAW->get_FRS_tofrr5();
    sci_tof5 = RAW->get_FRS_tof5();

    ID_x2 = RAW->get_FRS_x2();
    ID_y2 = RAW->get_FRS_y2();
    ID_a2 = RAW->get_FRS_a2();
    ID_b2 = RAW->get_FRS_b2();

    ID_x4 = RAW->get_FRS_x4();
    ID_y4 = RAW->get_FRS_y4();
    ID_a4 = RAW->get_FRS_a4();
    ID_b4 = RAW->get_FRS_b4();
   

    for(int i =0; i<7; i++){
    TPC_X[i] = RAW-> get_FRS_tpcX(i);
    TPC_Y[i] = RAW-> get_FRS_tpcY(i);

    
    for(int j=0; j<2; j++){
    TPC_LT[i][j] = RAW->get_FRS_tpclt(i,j);
    TPC_RT[i][j] = RAW->get_FRS_tpcrt(i,j);
        }

    }
    TPC_X0 = RAW->get_FRS_tpcx0();
    TPC_X1 = RAW->get_FRS_tpcx0();

    sci_dt_21l_21r = RAW->get_FRS_dt_21l_21r();
    sci_dt_41l_41r = RAW->get_FRS_dt_41l_41r();
    sci_dt_42l_42r = RAW->get_FRS_dt_42l_42r();
    sci_dt_43l_43r = RAW->get_FRS_dt_43l_43r();

    sci_dt_21l_41l = RAW->get_FRS_dt_21l_41l();
    sci_dt_21r_41r = RAW->get_FRS_dt_21r_41r();

    sci_dt_21l_42l = RAW->get_FRS_dt_21l_42l();
    sci_dt_21r_42r = RAW->get_FRS_dt_21r_42r();
    

    for(int k =0; k<2; ++k){
        ID_brho[k] = RAW->get_FRS_brho(k);
        ID_rho = RAW->get_FRS_rho(k);
  
    }
    
    hID_TACvsMHTDC_AoQ->Fill(RAW->get_FRS_id_mhtdc_aoq(), RAW->get_FRS_AoQ());
    
    hID_TACvsMHTDC_Z->Fill(RAW->get_FRS_id_mhtdc_z1(), RAW->get_FRS_z());
    
    hID_TACvsMHTDC_Beta->Fill(RAW->get_FRS_beta(),RAW->get_FRS_id_mhtdc_beta());
    
    hID_ZTACvsAoQMHTDC->Fill(RAW->get_FRS_id_mhtdc_aoq(),RAW->get_FRS_z());
    
    hID_ZMHTDCvsAoQTAC->Fill(RAW->get_FRS_id_mhtdc_z1(),RAW->get_FRS_AoQ());
    ///Using TAC
//      beta = RAW->get_FRS_beta();
//    
//   
//     gamma = RAW->get_FRS_gamma();
//     AoQ = RAW->get_FRS_AoQ();
//     AoQ_corr = RAW->get_FRS_AoQ_corr();
// 
//     ID_z = RAW->get_FRS_z();
//     ID_z2 = RAW->get_FRS_z2();
//     ID_z3 = RAW->get_FRS_z3();

 // beta3 = RAW->get_FRS_beta3();
    
    ///Using MHTDC
     AoQ = RAW->get_FRS_id_mhtdc_aoq();
     AoQ_corr = RAW->get_FRS_id_mhtdc_aoq_corr();

     ID_z = RAW->get_FRS_id_mhtdc_z1();
     ID_z2 = RAW->get_FRS_id_mhtdc_z2();
     beta = RAW->get_FRS_id_mhtdc_beta();
  // ID_z3 = RAW->get_FRS_z3();
 
    TRaw_vftx_21l = RAW->get_FRS_TRaw_vftx_21l();
    TRaw_vftx_21r = RAW->get_FRS_TRaw_vftx_21r();
    TRaw_vftx_22l = RAW->get_FRS_TRaw_vftx_22l();
    TRaw_vftx_22r = RAW->get_FRS_TRaw_vftx_22r();
    TRaw_vftx_41l = RAW->get_FRS_TRaw_vftx_41l();
    TRaw_vftx_41r = RAW->get_FRS_TRaw_vftx_41r();
    TRaw_vftx_42l = RAW->get_FRS_TRaw_vftx_42l();
    TRaw_vftx_42r = RAW->get_FRS_TRaw_vftx_42r();

    for(int i=0; i<32; i++){
   TRaw_vftx[i] = RAW->get_FRS_TRaw_vftx(i);
    }
    timestamp = RAW->get_FRS_timestamp();
    ts = RAW->get_FRS_ts(); //Spill time structrue
    ts2 = RAW->get_FRS_ts2();
    
     /// --------FRS SCALARS-------------------------------- //

    
    time_in_ms           = RAW->get_FRS_time_in_ms();
    spill_count          = RAW->get_FRS_spill_count();
    ibin_for_s           = RAW->get_FRS_ibin_for_s();
    ibin_for_100ms       = RAW->get_FRS_ibin_for_100ms();
    ibin_for_spill       = RAW->get_FRS_ibin_for_spill();
    ibin_clean_for_s     = RAW->get_FRS_ibin_clean_for_s();
    ibin_clean_for_100ms = RAW->get_FRS_ibin_clean_for_100ms();
    ibin_clean_for_spill = RAW->get_FRS_ibin_clean_for_spill();
    increase_scaler_temp = RAW->get_FRS_increase_scaler_temp();
    static bool scalers_done = false;
    

    /// ------------MUSIC---------------------------- //

    
     //MUSIC 1 is TUM MUSIC (8 anodes). MUSIC 3 not required
    /*for(int i=0; i<3; i++){
   //  hMUSIC1_MUSIC2->Fill(Music_dE[0],Music_dE[1]);
   //  hMUSIC1_dE1dE2->Fill(Music_E1[0],Music_E1[1]);
    }*/
    if(PrcID==20){
    for(int i=0; i<8; i++){
   if(Music_E1[i]!=0) hMUSIC1_E[i]->Fill(Music_E1[i]);
   if(Music_E2[i]!=0) hMUSIC2_E[i]->Fill(Music_E2[i]);
   if(Music_T1[i]!=0) hMUSIC1_T[i]->Fill(Music_T1[i]);
   if(Music_T2[i]!=0) hMUSIC2_T[i]->Fill(Music_T2[i]);
 }
    }
    //SCI
     /*-------------------------------------------------------------------------*/
    /* focus index: detector number                  tof index  tof path       */
    /*       0:     Sc01                                0:     TA - S1         */
    /*       1:     Sc11                                1:     S1 - S2         */
    /*       2:     Sc21                                2:     S2 - S41        */
    /*       3:     Sc21                                3:     S2 - S42        */
    /*       4:     Sc31                                4:     S2 - 81         */
    /*       5:     Sc41                                5:     S2 - E1         */
    /*                                                                         */
    /*       6:     Sc42                              tof index not used up to */
    /*       7:     Sc43 (previously Sc51)             now, only separate      */
    /*       8:     Sc61                              variables for S2-S41 and */
    /*       9:     ScE1 (ESR)                                S2-S42           */
    /*      10:     Sc81                                                       */
    /*      11:     Sc82                                                       */
    /*-------------------------------------------------------------------------*/
    for (int cnt=0;cnt<6;cnt++) //
     {
       int idx = 0 ;
       //int mw_idx = 0;
       //Float_t mwx = 0;
       switch(cnt)
     {
     case 0:        /* SC21 */
       idx = 2;
       //mw_idx = 2;
       //mwx = clb.sc21_x;
       break;
     case 1:        /* SC21 delayed */
       idx = 3;
       //mw_idx = 2;
       //mwx = clb.sc21_x;
       break;
     case 2:        /* SC41 */
       idx = 5;
       //mw_idx = 5;
       //mwx = clb.tpc_sc41_x;
       break;
     case 3:        /* SC42 */
           idx = 6;
       break;
     case 4:
       idx = 7;     /* SC43 */
       break;
     case 5:
       idx = 10;    /* SC81 */
       break;
     default: idx = 2;
     }
   
     
    if(PrcID==10){
     
     
      scalers_done = false;
      if(sci_l[idx]!=0)   hSCI_L[idx]->Fill(sci_l[idx]); 
      if(sci_r[idx]!=0)   hSCI_R[idx]->Fill(sci_r[idx]);
      if(sci_e[idx]!=0)   hSCI_E[idx]->Fill(sci_e[idx]);
        if(sci_e[2]!=0 && sci_e[5]!=0)    hSCI_dE24->Fill(sci_e[2],sci_e[5]);
            }
      
    if(PrcID==30){
      if(sci_tx[idx]!=0)  hSCI_Tx[idx]->Fill(sci_tx[idx]);
      if(sci_x[idx]!=0)   hSCI_X[idx]->Fill(sci_x[idx]);
       }
     }
     
//      if(sci_e[5]!=0 && TPC_X[5]!=0) hSCIdE41_TPC42X->Fill(sci_e[5],TPC_X[5]); //dE_SCI_41 vs TPC_42X
//      if(sci_l[5]!=0 &&TPC_X[5]!=0) hSCIdE41L_TPC42X->Fill(sci_l[5],TPC_X[5]); //dE_SCI_41L vs TPC_42X
//      if(sci_l[5]!=0 &&TPC_X[4]!=0) hSCIdE41L_TPC41X->Fill(sci_l[5],TPC_X[4]); //dE_SCI_41L vs TPC_41X
//      if(sci_r[5]!=0 &&TPC_X[5]!=0) hSCIdE41R_TPC42X->Fill(sci_r[5],TPC_X[5]); //dE_SCI_41R vs TPC_42X
//      if(sci_r[5]!=0 &&TPC_X[4]!=0) hSCIdE41R_TPC41X->Fill(sci_r[5],TPC_X[4]); //dE_SCI_41R vs TPC_41X
//       //cout<<"sci_e[0] " << sci_e[0] <<endl;
//      if(sci_e[0]!=0 &&TPC_X[5]!=0) hSCIdE21_TPC42X->Fill(sci_e[0],TPC_X[5]); //dE_SCI_21 vs TPC_42X
//      if(sci_l[0]!=0 &&TPC_X[5]!=0) hSCIdE21L_TPC42X->Fill(sci_l[0],TPC_X[5]); //dE_SCI_21L vs TPC_42X
//      if(sci_l[0]!=0 &&TPC_X[4]!=0) hSCIdE21L_TPC41X->Fill(sci_l[0],TPC_X[4]); //dE_SCI_21L vs TPC_41X
//      if(sci_l[0]!=0 &&TPC_X[4]!=0) hSCIdE21R_TPC42X->Fill(sci_l[0],TPC_X[4]); //dE_SCI_21R vs TPC_42X
//      if(sci_r[0]!=0 &&TPC_X[4]!=0) hSCIdE21R_TPC41X->Fill(sci_r[0],TPC_X[4]); //dE_SCI_21R vs TPC_41X

    if(PrcID==25){
     if(sci_tofll2!=0)  hSCI_TofLL2->Fill(sci_tofll2);
     if(sci_tofll3!=0) hSCI_TofLL3->Fill(sci_tofll3);
  //   hSCI_Tof2->Fill(sci_tof2);
     if(sci_tofrr2!=0)   hSCI_TofRR2->Fill(sci_tofrr2);
     if(sci_tofrr3!=0)  hSCI_TofRR3->Fill(sci_tofrr3);
     if(sci_tof3!=0)  hSCI_Tof3->Fill(sci_tof3);
     
     if(sci_tofrr5!=0)   hSCI_TofLL5->Fill(sci_tofll5);
     if(sci_tofrr5!=0)   hSCI_TofRR5->Fill(sci_tofrr5);
     if(sci_tof5!=0)  hSCI_Tof5->Fill(sci_tof5);
     
    

    if(ID_x2!=0) hID_x2->Fill(ID_x2);
    if(ID_y2!=0) hID_y2->Fill(ID_y2);
    if(ID_a2!=0) hID_a2->Fill(ID_a2);
    if(ID_b2!=0) hID_b2->Fill(ID_b2);

    if(ID_x4!=0) hID_x4->Fill(ID_x4);
    if(ID_y4!=0) hID_y4->Fill(ID_y4);
    if(ID_a4!=0) hID_a4->Fill(ID_a4);
    if(ID_b4!=0) hID_b4->Fill(ID_b4);
    

    for(int i=0; i<7; i++){
       if(TPC_X[i]!=0)   hTPC_X[i]->Fill(TPC_X[i]);
       if(TPC_Y[i]!=0)    hTPC_Y[i]->Fill(TPC_Y[i]);
       if(TPC_X[i]!=0 && TPC_Y[i]!=0)    hcTPC_XY[i]->Fill(TPC_X[i],TPC_Y[i]);
       if(TPC_LT[i][0]!=0 && TPC_RT[i][1]!=0)    hTPC_LTRT[i]->Fill(TPC_LT[i][0],TPC_RT[i][1]);
       if(TPC_X0-TPC_X1!=0)    hTPC_DELTAX[i]->Fill(TPC_X0-TPC_X1);
         // cout<<"TPC_LT[i][0] " << TPC_LT[i][0] << " i " << i << endl;
    }


//     htpc_X2->Fill(TPC_X[2]);
//     htpc_Y2->Fill(TPC_Y[2]);
//     htpc_X4->Fill(TPC_X[4]);
//     htpc_Y4->Fill(TPC_Y[4]);
    for(int i=0;i<2;i++){
      if(ID_brho[i]!=0)hID_BRho[i]->Fill(ID_brho[i]);
    }

    //SCI tx
    //if(sci_dt_21l_21r) hSCI_dT_21l_21r->Fill(sci_dt_21l_21r);
//     if(sci_dt_41l_41r) hSCI_dT_41l_41r->Fill(sci_dt_41l_41r);
//     if(sci_dt_42l_42r) hSCI_dT_42l_42r->Fill(sci_dt_42l_42r);

    if(sci_dt_21l_41l!=0) hSCI_dT_21l_41l->Fill(sci_dt_21l_41l);
    if(sci_dt_21r_41r!=0) hSCI_dT_21r_41r->Fill(sci_dt_21r_41r);

    if(sci_dt_21l_42l!=0) hSCI_dT_21l_42l->Fill(sci_dt_21l_42l);
    if(sci_dt_21r_42r!=0) hSCI_dT_21r_42r->Fill(sci_dt_21r_42r);



   
   // if(beta3) hbeta3->Fill(beta3);
}
 
if(PrcID==25){
  
   if(beta!=0) hID_beta->Fill(beta*1000);
    if(AoQ!=0) hID_AoQ->Fill(AoQ);
   
    if(AoQ_corr!=0) hID_AoQ_corr->Fill(AoQ_corr);

    /****  S4  (MUSIC 1)   */
     if(ID_z!=0)hID_Z->Fill(ID_z);
     /****  S4  (MUSIC 2)   */
     if(ID_z2!=0) hID_Z2->Fill(ID_z2);
    
     /****  S4  (MUSIC OLD)   */
     //hID_Z3->Fill(ID_z3);

//      hID_Z_Z2->Fill(ID_z,ID_z2);
     if(ID_z!=0 && Music_dE[1]!=0)hID_Z_dE2->Fill(ID_z,Music_dE[1]);
    // hID_Z_Z3->Fill(ID_z,ID_z3);
     //if(ID_z!=0 && sci_l[2]!=0 && sci_r[2]!=0)hID_Z_Sc21E->Fill(ID_z, sqrt(sci_l[2]*sci_r[2]));

     if(ID_x2!=0&&ID_x4!=0 ) hID_x2x4->Fill(ID_x2, ID_x4);
     if(AoQ!=0 && sci_e[5]!=0) hID_SC41dE_AoQ->Fill(AoQ, sci_e[5]);

     if(sci_tof2!=0 && Music_dE[0]!=0) hID_dEToF->Fill(sci_tof2, Music_dE[0]);

     if(ID_z!=0 && ID_x2!=0) hID_x2z->Fill(ID_z, ID_x2);// MUSIC1
     if(ID_z!=0 && ID_x4!=0) hID_x4z->Fill(ID_z, ID_x4);// MUSIC1

     if(ID_x4!=0 && Music_dE[0]!=0) hID_E_Xs4->Fill(ID_x4,Music_dE[0]);
     if(ID_x4!=0 && Music_dE[0]!=0)hID_E_Xs2->Fill(ID_x2,Music_dE[0]);

     if(ID_x2!=0 && ID_a2!=0)hID_x2a2->Fill(ID_x2,ID_a2);
     if(ID_y2!=0 && ID_b2!=0)hID_y2b2->Fill(ID_y2,ID_b2);
     if(ID_x4!=0 && ID_a4!=0)hID_x4a4->Fill(ID_x4,ID_a4);
     if(ID_x4!=0 && ID_b4!=0) hID_y4b4->Fill(ID_y4,ID_b4);
}
    if (PrcID == 35 && !scalers_done)
    {
    if(timestamp) htimestamp->Fill(timestamp);
    if(ts) hts->Fill(ts);
    if(ts2) hts2->Fill(ts2);
    
    for(int ii=0; ii<64; ii++){
    
   
    //printf("ch %d: this event = %lld, increase =%lld\n",ii,src.sc_long[ii],increase_scaler_temp);
    //hScaler_per_s[ii]->AddBinContent(ibin_for_s, increase_scaler_temp[ii]);
//    hScaler_per_100ms[ii]->AddBinContent(ibin_for_100ms, increase_scaler_temp[ii]);
   // hScaler_per_spill[ii]->AddBinContent(ibin_for_spill, increase_scaler_temp[ii]);    
   // if(ii=50)cout<<"ibin_clean_for_s " << ibin_clean_for_s << " increase_scaler_temp " << increase_scaler_temp<< endl;
   frs_scaler_value[ii] += increase_scaler_temp[ii];
  }
  
   for(int ii=0; ii<64; ii++){
    hScaler_per_s[ii]->SetBinContent(ibin_clean_for_s, 0);
   // hScaler_per_100ms[ii]->SetBinContent(ibin_clean_for_100ms, 0);
 //   hScaler_per_spill[ii]->SetBinContent(ibin_clean_for_spill, 0);
    }
   scalers_done = true;
   }
   
   ///VFTX histogramming
  if(PrcID==35){
      ///T_RAW  
      for(int i=0; i<32; i++){
      if(TRaw_vftx[i]!=0){
      hvftx_TRaw[i]->Fill(TRaw_vftx[i]);
     //cout<<"TRaw_vftx[i] " <<TRaw_vftx[i] << " i " << i << endl;
      }          
    }
      
      
  ///  SCI21 Position
   if(TRaw_vftx_21l!=0. && TRaw_vftx_21r!=0.){
   hvftx_Sci21PosRaw->Fill((TRaw_vftx_21l-TRaw_vftx_21r)); ///1ps

   }
    /// SCI22 Position
     if(TRaw_vftx_22l!=0. && TRaw_vftx_22r!=0.){
   hvftx_Sci22PosRaw->Fill((TRaw_vftx_22l-TRaw_vftx_22r)); ///1ps
    }
    /// SCI41 Position
      if(TRaw_vftx_41l!=0. && TRaw_vftx_41r!=0.){
   hvftx_Sci41PosRaw->Fill((TRaw_vftx_41l-TRaw_vftx_41r)); ///1ps

    }
    /// SCI42 Position
     if(TRaw_vftx_42l!=0. && TRaw_vftx_42r!=0.){
   hvftx_Sci42PosRaw->Fill((TRaw_vftx_42l-TRaw_vftx_42r)); ///1ps
    }
    
   ///SCI41 -SCI21 ToF
  if(TRaw_vftx_41l!=0.&&TRaw_vftx_41r!=0.&&TRaw_vftx_21l!=0.&&TRaw_vftx_21r!=0.){
   hvftx_ToFraw_2141->Fill((0.5*((TRaw_vftx_41l+TRaw_vftx_41r)-(TRaw_vftx_21l+TRaw_vftx_21r))));//ps
    }
    ///SCI41 -SCI21 LL
  if(TRaw_vftx_41l!=0.&&TRaw_vftx_21l!=0.){
   hvftx_ToFraw_2141LL->Fill((TRaw_vftx_41l-TRaw_vftx_21l));//ps
    }
    
    ///SCI41 -SCI21 RR
  if(TRaw_vftx_41r!=0.&& TRaw_vftx_21r!=0.){
   hvftx_ToFraw_2141RR->Fill((TRaw_vftx_41r-TRaw_vftx_21r));//ps
  //cout<<"hvftx_ToFraw_2141 " <<(0.5*(TRaw_vftx_41l+TRaw_vftx_41r-(TRaw_vftx_21l+TRaw_vftx_21r)))<< endl;
    }
    
     ///SCI42 -SCI21 LL
  if(TRaw_vftx_42l!= 0 &&TRaw_vftx_21l!= 0){
   hvftx_ToFraw_2142LL->Fill((TRaw_vftx_42l-TRaw_vftx_21l));//ps
  //cout<<"hvftx_ToFraw_2141 " <<(0.5*(TRaw_vftx_41l+TRaw_vftx_41r-(TRaw_vftx_21l+TRaw_vftx_21r)))<< endl;
    }
    
    ///SCI42 -SCI21 RR
  if(TRaw_vftx_42r!=0 && TRaw_vftx_21r!=0){
   hvftx_ToFraw_2141RR->Fill((TRaw_vftx_42r-TRaw_vftx_21r));//ps
  //cout<<"hvftx_ToFraw_2141 " <<(0.5*(TRaw_vftx_41l+TRaw_vftx_41r-(TRaw_vftx_21l+TRaw_vftx_21r)))<< endl;
        }
    }
  }

  /**----------------------------------------------------------------------------------------------**/
  /**-------------------------------------------  AIDA   ------------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/

  void EventUnpackProc::Make_AIDA_Histos(){

    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
    hAIDA_ADC.resize(conf->FEEs());

    for (int i = 0; i < conf->FEEs(); i++)
    {
      for (int j = 0; j < 64; j++)
      {
        hAIDA_ADC[i][j][0] = MakeTH1('I',
          Form("AIDA/Unpacker/FEE%d/Fee%d_L_Channel%02d", i+1, i+1, j+1),
          Form("FEE %d Channel %2d (Low Energy)", i+1, j+1),
          2000, -32768, 32767
        );
      }
    }

    for (int i = 0; i < conf->FEEs(); i++)
    {
      for (int j = 0; j < 64; j++)
      {
        hAIDA_ADC[i][j][1] = MakeTH1('I',
          Form("AIDA/Unpacker/FEE%d/Fee%d_H_Channel%02d", i+1, i+1, j+1),
          Form("FEE %d Channel %2d (High Energy)", i+1, j+1),
          2000, -32768, 32767
        );
      }
    }
}

void EventUnpackProc::Fill_AIDA_Histos() {
  AIDA_Hits = RAW->get_AIDA_HITS();

  for(int i = 0; i<AIDA_Hits; i++) {
    int fee = RAW-> get_AIDA_FEE_ID(i);
    int chan = RAW-> get_AIDA_CHA_ID(i);
    int adc = RAW->get_AIDA_ADC(i);
    int veto = RAW->get_AIDA_HighE_VETO(i) ? 1 : 0;

    hAIDA_ADC[fee][chan][veto]->Fill(adc - 32767);

    //cout<<"chan " << chan << endl;
  }
}

/**----------------------------------------------------------------------------------------------**/
/**-----------------------------------------  FATIMA TAMEX  ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
 void EventUnpackProc::Make_FATIMA_TAMEX_Histos(){


        for(int i=0; i<50; i++){
             hFATlead_Coarse[i]= MakeTH1('D', Form("FATIMA_TAMEX/Lead_Coarse/Lead-CoarseCh.%02d", i), Form("Lead Coarse %2d", i), 50000, -100000., 100000.);
             hFATlead_Fine[i]= MakeTH1('D', Form("FATIMA_TAMEX/Lead_Fine/Lead-FineCh.%02d", i), Form("Lead Fine %2d", i), 601, -1., 60.);
             hFATtrail_Coarse[i]= MakeTH1('D', Form("FATIMA_TAMEX/Trail_Coarse/Trail-CoarseCh.%02d", i), Form("Trail Coarse %2d", i), 50000, -100000., 100000.);
             hFATtrail_Fine[i]= MakeTH1('D', Form("FATIMA_TAMEX/Trail_Fine/Trail-FineCh.%02d", i), Form("Trail Fine %2d", i), 601, -1., 60.);


            }
        }
  //-----------------------------------------------------------------------------------------------------------------------------//

//   void EventUnpackProc::Fill_FATIMA_TAMEX_Histos(){
// 
//       ///TAMEX
//     //get amount of fired Tamex modules
//     int TamexHits_Fatima = RAW->get_FATIMA_tamex_hits();
// 
//     int Lead_Fatima_Coarse[4][32];
//     double Lead_Fatima_Fine[4][32];
//     int Trail_Fatima_Coarse[4][32];
//     double Trail_Fatima_Fine[4][32];
//     int Phys_Channel_Fatima[32];
//     int leadHits_Fatima = 0,leadHitsCh_Fatima = 0;
//     int trailHits_Fatima = 0,trailHitsCh_Fatima = 0;
// 
// 
//     for(int i=0; i<32; i++){
//      Phys_Channel_Fatima[i] = 0;
//     }
// 
//     for(int i =0; i<4; i++){
//         for(int j=0; j<32;j++){
//      Lead_Fatima_Coarse[i][j] = 0;
//      Lead_Fatima_Fine[i][j] = 0;
//      Trail_Fatima_Coarse[i][j] = 0;
//      Trail_Fatima_Fine[i][j] = 0;
//         }
// 
//     for(int i = 0;i < TamexHits_Fatima;++i){
// 
//       leadHits_Fatima = RAW->get_FATIMA_lead_hits(i);
//       trailHits_Fatima = RAW->get_FATIMA_trail_hits(i);
// 
// 
//         //Box diagrams for leading and trailing
//       for(int j = 0;j < RAW->get_FATIMA_am_Fired(i);j++){
//          if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 1){
//             Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
//           //  Lead_Fatima[i][j] = RAW->get_FATIMA_lead_T(i,j);
//             Lead_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_lead(i,j);
//             Lead_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_lead(i,j);
//             hFATlead_Coarse[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Coarse[i][j]);
//             hFATlead_Fine[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Fine[i][j]);
// 
//         }
//          if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 0){
//             Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
//           //  Trail_Fatima[Phys_Channel_Fatima[j]] = RAW->get_FATIMA_trail_T(i,j);
//             Trail_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_trail(i,j);
//             Trail_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_trail(i,j);
//             hFATtrail_Coarse[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Coarse[i][j]);
//             hFATtrail_Fine[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Fine[i][j]);
//                 }
//             }
//         }
//     }
// 
//   }
  /**----------------------------------------------------------------------------------------------**/
/**-----------------------------------------  FATIMA VME  ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventUnpackProc::Make_FATIMA_VME_Histos(){


  for (int i=0; i<FAT_MAX_VME_DET; i++){
    hFAT_Eraw_VME[i] = MakeTH1('D', Form("FATIMA_VME/Unpacker/Energy/Raw/E_Raw_LaBrCh. %02d", i),
    Form("LaBr%02d energy (raw)", i),2000,0,40000);

    hFAT_Traw_VME[i] = MakeTH1('D', Form("FATIMA_VME/Unpacker/Timing/Raw/Traw_LaBrCh. %02d", i),
                    Form("LaBr%02d energy", i),5000,0,4E6);

            }
    hScalar_hit_pattern = MakeTH1('D',"Scalar/HitPat","Scalar Hit pattern",32,0,32);
        }
//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::Fill_FATIMA_VME_Histos(){
  double FAT_E[50],FAT_T[50];
  double En_i;
  int detQDC, detTDC;
  detTDC = 0;
  detQDC = 0;
  for (int k=0; k<50; k++){
    FAT_T[k] = 0;
    FAT_E[k] = 0;
  }
  int Scalar_iterator = 0;
  int Scalar_Chan = 0;


  /**------------------FATIMA Energy -----------------------------------------**/
  for (int i=0; i<RAW->get_FAT_QDCs_fired(); i++){ /** Loops over only channels in the QDC **/

    detQDC = RAW->get_FAT_QDC_id(i); /**FAT ID QDC*/
    En_i = RAW->get_FAT_QLong_Raw(i); /**Raw FAT Energy*/
  if (detQDC<40){
   hFAT_Eraw_VME[detQDC]->Fill(En_i);
    //cout << " EN " << En_i<<endl;
  }

  }
  /**------------------FATIMA TIMING -----------------------------------------**/
  for (int i=0; i<RAW->get_FAT_TDCs_fired(); i++){ /** Loops over only channels in the TDC 1-4 **/

   
    detTDC = (RAW->get_FAT_TDC_id(i));
     if(detTDC<50){
        FAT_T[detTDC] = (RAW->get_FAT_TDC_timestamp(i));
        hFAT_Traw_VME[detTDC]->Fill(FAT_T[detTDC]*25); //in ps
  // if(detTDC==40)  cout<<"FAT_T[detTDC]*25 " <<FAT_T[detTDC]*25 <<" detTDC " <<detTDC << endl;
                }
         }
  /**------------------Scaler TIMING -----------------------------------------**/  
         Scalar_iterator = RAW->get_scalar_iterator();
         
          for (int g=0; g<Scalar_iterator; g++){
              if(RAW->get_scalar_data(g)>0){
                hScalar_hit_pattern->Fill(g);
           }
        }
    }

void EventUnpackProc::Fill_FATIMA_TAMEX_Histos(){

    //get amount of fired Tamex modules
    int TamexHits_Fatima = RAW->get_FATIMA_tamex_hits();

    int Lead_Fatima_Coarse[4][32];
    double Lead_Fatima_Fine[4][32];
    int Trail_Fatima_Coarse[4][32];
    double Trail_Fatima_Fine[4][32];
    int Phys_Channel_Fatima[32];
    int leadHits_Fatima = 0,leadHitsCh_Fatima = 0;
    int trailHits_Fatima = 0,trailHitsCh_Fatima = 0;

    for(int i=0; i<32; i++){
     Phys_Channel_Fatima[i] = 0;
    }

    for(int i =0; i<4; i++){
        for(int j=0; j<32;j++){
     Lead_Fatima_Coarse[i][j] = 0;
     Lead_Fatima_Fine[i][j] = 0;
     Trail_Fatima_Coarse[i][j] = 0;
     Trail_Fatima_Fine[i][j] = 0;
        }
    }
        ///QDC
    double FAT_E[50];
    double En_i;
    int detQDC;
    detQDC = 0;
    for (int k=0; k<50; k++){

        FAT_E[k] = 0;
  }
 //cout<<"TamexHits_Fatima "<<TamexHits_Fatima<<endl;
    for(int i = 0;i < TamexHits_Fatima;++i){

      leadHits_Fatima = RAW->get_FATIMA_lead_hits(i);
      trailHits_Fatima = RAW->get_FATIMA_trail_hits(i);


        //Box diagrams for leading and trailing
      for(int j = 0;j < RAW->get_FATIMA_am_Fired(i);j++){
         if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 1){
            Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
          //  Lead_Fatima[i][j] = RAW->get_FATIMA_lead_T(i,j);
            Lead_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_lead(i,j);
            Lead_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_lead(i,j);
            hFATlead_Coarse[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Coarse[i][j]);
            hFATlead_Fine[Phys_Channel_Fatima[j]]->Fill(Lead_Fatima_Fine[i][j]);

        }
         if(RAW->get_FATIMA_CH_ID(i,j) % 2 == 0){
            Phys_Channel_Fatima[j] = RAW->get_FATIMA_physical_channel(i,j);
          //  Trail_Fatima[Phys_Channel_Fatima[j]] = RAW->get_FATIMA_trail_T(i,j);
            Trail_Fatima_Coarse[i][j] = RAW->get_FATIMA_coarse_trail(i,j);
            Trail_Fatima_Fine[i][j] = RAW->get_FATIMA_fine_trail(i,j);
            hFATtrail_Coarse[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Coarse[i][j]);
            hFATtrail_Fine[Phys_Channel_Fatima[j]]->Fill(Trail_Fatima_Fine[i][j]);
                }
            }
        }
  }
/**----------------------------------------------------------------------------------------------**/
/**----------------------------------------   GALILEO   -----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/


void EventUnpackProc::Make_GALILEO_Histos(){
//   for (int j; j<GALILEO_MAX_HITS; j++){
//         hGAL_Raw_E[j] = MakeTH1('D',Form("GALILEO/Raw/GALILEO_Energy_Spectra/GALILEO_Raw_E%2d",j),
//                             Form("GALILEO Channel Energy Channel Raw %2d",j),20000,0,20000);
// 
//                     }
    hFebTime  = MakeTH1('D',"SysTime/FebexClock","Febex clock",200000,0,200000);
                }
//-----------------------------------------------------------------------------------------------------------------------------//
void EventUnpackProc::Fill_GALILEO_Histos(){

    double tmpGAL[32];
    int  GALILEO_hits, GalID;

     /**------------------GALILEO Raw Energy -----------------------------------------**/
      GALILEO_hits = RAW->get_GALILEO_am_Fired();
         for(int i=0; i<GALILEO_hits; i++){
	   if(RAW->get_GALILEO_Det_id(i)>-1){
        hFebTime->Fill(RAW->get_GALILEO_Chan_T(i)*10E-9);
        GalID = RAW->get_GALILEO_Det_id(i) * 3 + RAW->get_GALILEO_Crystal_id(i);
        tmpGAL[GalID] = RAW->get_GALILEO_Chan_E(i);
      
         }
	 }
   }

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::checkTAMEXorVME(){

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



//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::print_MBS(int* pdata,int lwords){
  cout << "---------------------\n";
  for(int i = 0;i < lwords;++i){
    cout << hex << *(pdata + i) << " ";
    if(i % 5 == 0 && i > 0) cout << endl;
  }
  cout << "\n---------------------\n";
}
//-----------------------------------------------------------------------------------------------------------------------------//

TH1I* EventUnpackProc::MakeH1I(const char* fname,
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

TH2I* EventUnpackProc::MakeH2I(const char* fname,
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
TH1I* EventUnpackProc::MakeH1I_TPC(const char* foldername, const char* name, int nameindex,
                  Int_t nbinsx, Float_t xmin, Float_t xmax,
                  const char* xtitle, Color_t linecolor, Color_t fillcolor)
{
  char fullname[100];
  if(nameindex>=0)
    sprintf(fullname,"%s%s",tpc_name_ext1[nameindex],name);
  else
    strcpy(fullname, name);
  return MakeH1I(foldername, fullname, nbinsx, xmin, xmax, xtitle,
         linecolor, fillcolor);
}
const  char* EventUnpackProc::tpc_name_ext1[7]={"TPC21_","TPC22_","TPC23_","TPC24_","TPC41_","TPC42_", "TPC31_"};
const  char* EventUnpackProc::tpc_folder_ext1[7]={"TPC21","TPC22","TPC23","TPC24","TPC41","TPC42","TPC31"};


//-----------------------------------------------------------------------------------------------------------------------------//
//                                                            END                                                              //
//-----------------------------------------------------------------------------------------------------------------------------//
