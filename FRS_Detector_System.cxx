#include "FRS_Detector_System.h"

#include "Riostream.h"
#include "TROOT.h"
#include "TH1.h"
#include "TMap.h"
#include "TObject.h"
#include "TString.h"
#include "TObjString.h"
 
#include <cstdlib>
#include <bitset>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>

// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "TGo4MbsEvent.h"
#include "TGo4StepFactory.h"
#include "TGo4Analysis.h"
#include "TGo4AnalysisObjectManager.h"
#include "TGo4WinCond.h"
#include "TGo4PolyCond.h"
#include "TGo4CondArray.h"
#include "TGo4EventProcessor.h"
#include "TGo4Parameter.h"

#include "TFRSParameter.h"

#include "DESPECAnalysis.h"
#include "TH2.h"
#include "TCutG.h"

#include "Gtypes.h"

#include <sys/stat.h>
#include <sys/types.h>
//#include "TFRSVftxSetting.h"


using namespace std;

//---------------------------------------------------------------


FRS_Detector_System::FRS_Detector_System(){
    
    
     DESPECAnalysis* an = dynamic_cast<DESPECAnalysis*> (TGo4Analysis::Instance());
    
  if (an==0) {
    cout << "!!!  Script should be run in FRS analysis" << endl;
    return;
  }
  
   frs = dynamic_cast<TFRSParameter*> (an->GetParameter("FRSPar"));
  if (frs==0) {
    cout << "!!!  Parameter FRSPar not found" << endl;
    return;
  }
  
 //  ModSetup = dynamic_cast<TModParameter*>(an->GetParameter("ModPar"));
  
   mw = dynamic_cast<TMWParameter*> (an->GetParameter("MWPar"));
  if (mw==0) {
    cout << "!!!  Parameter MWPar not found" << endl;
    return;
  }

   music = dynamic_cast<TMUSICParameter*> (an->GetParameter("MUSICPar"));
  if (music==0) {
    cout << "!!!  Parameter MUSICPar not found" << endl;
    return;
  }

  sci = dynamic_cast<TSCIParameter*> (an->GetParameter("SCIPar"));
  if (sci==0) {
    cout << "!!!  Parameter SCIPar not found" << endl;
    return;
  }

   id = dynamic_cast<TIDParameter*> (an->GetParameter("IDPar"));
  if (id==0) {
    cout << "!!!  Parameter IDPar not found" << endl;
    return;
  }
   
   tpc = dynamic_cast<TTPCParameter*> (an->GetParameter("TPCPar"));
  if (tpc==0) {
    cout << "!!!  Parameter TPCPar not found" << endl;
    return;
  }

   si = dynamic_cast<TSIParameter*> (an->GetParameter("SIPar"));
  if (si==0) {
    cout << "!!!  Parameter SIPar not found" << endl;
    return;
  }

   mrtof = dynamic_cast<TMRTOFMSParameter*> (an->GetParameter("MRTOFMSPar"));
  if (mrtof==0) {
    cout << "!!!  Parameter MR-TOF-MSPar not found" << endl;
    return;
  }
  
  //TModParameter* ElecMod = dynamic_cast<TModParameter*>(an->GetParameter("ModPar"));
  if(frs == nullptr)
    {
     std::cout<<"E> FRS parameters not set properly, it is nullptr !"<<std::endl;
     std::exit(-1);
    }
    #if CALIBRATION_VFTX
  for(int mod=0; mod<VFTX_N; mod++)
    for(int ch=0; ch<VFTX_MAX_CHN; ch++)
      Vftx_Stat[mod][ch]=0;
#endif
    
    Text_t hist[256];
    
    for (int module=0; module<VFTX_N; module++) {
    for(int channel=0; channel<VFTX_MAX_CHN; channel++){
    sprintf(hist,"FRS/VFTX_Unpack/VFTX_%i/FineTime/VFTX%i_FineTime_ch%02d",module,module,channel);
    h1_vftx_ft[module][channel] = new TH1D(hist,hist,1000,0,1000);
   
    
        }
    }
      // VFTX parameters
  m_VFTX_Bin2Ps();
//     ///for VFTX calibration
//    min_val=0;
//    max_val=1000;
//    nbins =1000;
//     #if LOAD_CALIBRATION_VFTX
// //   
//     load_VFTX_Calibration_Files();
//      #endif
  //  vme2s = new Int_t*[32];           // User TDC (V1290) 
    //vme2s_trailing = new Int_t*[32];  // User TDC (V1290) 
    
     //leading_v1290_main     =   new UInt_t*[32];             // Mtof TDC (V1290) 
     //trailing_v1290_main   = new UInt_t*[32];   // Mtof TDC (V1290) 

    nhit5 = new Int_t*[32];            // multiplicity (V1290)
    
    //#if NEW_CALIBRATION_VFTX
  for(int mod=0; mod<VFTX_N; mod++)
    for(int ch=0; ch<VFTX_MAX_CHN; ch++)
      Vftx_Stat[mod][ch]=0;
 // #endif
    
    for(int i = 0; i < 32; ++i){
    
//     vme2s[i]        = new Int_t[10];           // User TDC (V1290) 
//     vme2s_trailing[i]   = new Int_t[10];  // User TDC (V1290) 

//   leading_v1290_main[i]         = new UInt_t[10];           // Mtof TDC (V1290) 
  // trailing_v1290_main[i]    = new UInt_t[10];  // Mtof TDC (V1290) 

    nhit5[i]        = new Int_t[2];       // multiplicity (V1290)
//    nhit_v1290_main[i]    = new Int_t[2];       // multiplicity (V1290)

    }
    
    
    vme2scaler = new Int_t[32];         // User Crate Messhute
    vme3scaler = new Int_t[32];         // User Crate Messhute

    firstTS = new bool[3];
    
    previousTimeStamp = new Long64_t[3];
    
    currentTimeStamp = 0;
    
    EventFlag = 0;
    
    qlength     = 0;         /* From event header     */
    qtype   = 0;           /*                       */
    qsubtype    = 0;        /*                       */
    qdummy  = 0;          /*                       */
    qtrigger    = 0;        /*                       */
    qevent_nr   = 0;       /*                       */
    
    for(int i=0; i<16; i++){
     vme_trmu_tdc[i] = 0;  
     vme_trmu_adc[i] = 0;
    }
    
 //   Proc_iterator=0;

    /*******************************************************************/
    /***SORT STUFF***/
  
    StartOfSpilTime     = 0;
    StartOfSpilTime2    = 0; //does not reset at end of extraction
    PreviousTS      = 0;
  
    counter = 0; 

    // time stamp data  
    ts_id = 0;                                
   // ts_word = new ULong64_t[4]; //for the titris time stammping
    //tsys_word = new Int_t[3]; //for the system time
    timestamp = 0;  // absolute time stamp value
    timespill = 0;  // relative time from start of the spill 
    timespill2 = 0;  // relative time from start of the spill does not reset at end extraction
    systemtime_ms = 0;
    systemtime_s = 0;

    pattern = 0;
    trigger = 0;

    // scaler readings     
    //sc_long = new UInt_t[64]; //changed from 32 to 64 (10.07.2018)
    sc_long2 = new UInt_t[32];


    // SEETRAM calib IC energy
    //      Int_t         ic_de;           /*                          */
      
     
    // part of MW parameter
    mw_an = new Int_t[13];       /*  anode time              */
    mw_xl = new Int_t[13];       /*  Rohdaten                */
    mw_xr = new Int_t[13];       /*                          */
    mw_yu = new Int_t[13];       /*                          */ 
    mw_yd = new Int_t[13];       /*                          */

    // TPC part //(HAPOL-25/03/06) 6 TPCs each one with 2 delay lines each!!    
    //7 TPCs (4 at S2, 2 at S4 and one at S3) 03.07.2018 SB
    
    tpc_l = new Int_t*[7];
    tpc_r = new Int_t*[7];
    tpc_lt = new Int_t*[7];
    tpc_rt = new Int_t*[7];
  //  tpc_timeref = new Int_t*[7];
    // [index][anode_no]
    tpc_dt = new Int_t*[7];
    tpc_a = new Int_t*[7];
    tpc_timeref = new Int_t[4];
    
    
    for(int i=0; i < 7; ++i){
    
    tpc_l[i] = new Int_t[2];
    tpc_r[i] = new Int_t[2];
    tpc_lt[i] = new Int_t[2];
    tpc_rt[i] = new Int_t[2];
    
    tpc_dt[i] = new Int_t[4];
    tpc_a[i] = new Int_t[4];
    
    }

      
    // SCI part
    de_21l  = 0;         /* dE SCI21 left            */  
    de_21r  = 0;         /* de SCI21 right           */ 
    de_41l  = 0;         /* dE SCI41 left            */
    de_41r  = 0;         /* dE SCI41 right           */
    de_41u  = 0;         /* dE SCI41 up              */
    de_41d  = 0;         /* dE SCI41 down            */
    de_42l  = 0;         /* de SCI42 left            */
    de_42r  = 0;         /* de SCI42 right           */
    de_43l  = 0;         /* de SCI43 left            */
    de_43r  = 0;         /* de SCI43 right           */
    de_81l  = 0;         /* de SCI81 left            */
    de_81r  = 0;         /* de SCI81 right           */  
//     de_21ld = 0;         /* dE SCI21 left delayed    */  
//     de_21rd = 0;         /* de SCI21 right delayed   */ 
//     de_v1l  = 0;         /* dE veto1 left            */  
//     de_v1r  = 0;         /* de veto1 right           */ 
//     de_v2l  = 0;         /* dE veto2 left            */
//     de_v2r  = 0;         /* dE veto2 right           */
//     de_v3   = 0;         /* dE veto3                 */


    dt_21l_21r  = 0;     /*                          */ 
    dt_41l_41r  = 0;     /*                          */ 
    dt_21l_41l  = 0;     /*                          */
    dt_21r_41r  = 0;     /*                          */
    dt_42l_42r  = 0;     /*                          */
    dt_43l_43r  = 0;     /*                          */
    dt_42l_21l  = 0;     /*                          */
    dt_42r_21r  = 0;     /*                          */
    dt_41u_41d  = 0;     /*                          */
    dt_81l_81r  = 0;     /*                          */
    dt_21l_81l  = 0;     /*                          */
    dt_21r_81r  = 0;     /*                          */
    dt_22l_22r  = 0;      /*                          */ 
    dt_22l_41l  = 0;
    dt_22r_41r  = 0;

    // User multihit TDC
        
    tdc_sc41l = new Int_t[10];
    tdc_sc41r = new Int_t[10];
    tdc_sc21l = new Int_t[10];
    tdc_sc21r = new Int_t[10];
    tdc_sc42l = new Int_t[10];
    tdc_sc42r = new Int_t[10];
    tdc_sc43l = new Int_t[10];
    tdc_sc43r = new Int_t[10];
    tdc_sc81l = new Int_t[10];
    tdc_sc81r = new Int_t[10];
    tdc_sc31l = new Int_t[10];
    tdc_sc31r = new Int_t[10];
    tdc_sc11 = new Int_t[10];
    // MUSIC1 part
    music_e1 = new Int_t[8];     /* Raw energy signals       */
    music_t1 = new Int_t[8];     /* Raw drift time           */
   // music_pres = new Int_t[3];   /* Music Druck              */
   // music_temp = new Int_t[3];   /* Music Temperatur         */
    
    // MUSIC2 part
    music_e2 = new Int_t[8];     /* Raw energy signals       */
    music_t2 = new Int_t[8];     /* Raw drift time           */
    
    // MUSIC3 part (OLD MUSIC)
    music_e3 = new Int_t[8];     /* Raw energy signals       */
    music_t3 = new Int_t[8];     /* Raw drift times          */

    
    /******************************************************************/
    /**CALIBRATION Parameters**/
    
    // MON data declarations
    
    check_first_event        = new int[2];
    scaler_channel_10khz = 0;//YT 20Jun
    scaler_time_count        = new Long64_t[2]; 
    scaler_spill_count       = new Long64_t[2]; //UInt_t
    scaler_time_check_last   = new Long64_t[2];//UInt_t
    scaler_spill_check_last      = new Long64_t[2];//UInt_t 
    check_increase_time      = new Long64_t[64];//UInt_t 
    check_increase_spill     = new Long64_t[64];//UInt_t
    scaler_increase_event    = new Long64_t[64];//UInt_t
    scaler_last_event        = new Long64_t[64];
    
    scaler_ch_1kHz=39; //ch7 of 2nd scaler
    scaler_ch_spillstart=8; //ch8 of 1st scaler 
    scaler_check_first_event=1;
    //TGo4WinCond   *cMW_XSUM[13];
    //TGo4WinCond   *cMW_YSUM[13];  
    
    ////condtition for control sum
    //TGo4WinCond *cTPC_CSUM[7][4];
    //TGo4PolyCond *cTPC_XY[7];
    
    //**** keep values from previous event
    focx_s2m        = 0;
    focy_s2m        = 0;
    angle_x_s2m     = 0;
    angle_y_s2m     = 0;
    
    
    
    fbFirstEvent = false;
    //unsigned int  scaler_save[64];
    scaler_save  = new Long64_t[64];
    
    firstsec        = 0;
    firsttenthsec   = 0;
    firsthundrethsec    = 0;    //mik   
    firstcycle      = 0;
    firstseetram    = 0;
    
    scalercycle_Sec         = 0;
    scalercycle_TenthSec    = 0;
    scalercycle_HundrethSec     = 0;
    scalercycle_Cycle       = 0;
    scalercycle_Seetram     = 0;
    
    dtime   = 0;
    dt_last     = 0;
    dt_diff     = 0;
    
    
    check_total_sc41        = 0;
    check_total_sc21        = 0;
    check_total_seetram     = 0;
    check_total_mrtof_start     = 0;
    check_total_mrtof_stop  = 0;
    
    
    freeTrig = 0;
    acptTrig = 0;
    mon_inc  = new Long64_t[64];
    
    // MON part
//     seconds      = 0;
//     tenthsecs        = 0;
//     hundrethsecs     = 0;    //mik
   // extraction_cycle  = 0;
    coin  = new Int_t[16];
    seetram         = 0;
    
    // MW part
    mw_xsum = new Float_t[13];     /*                          */
    mw_ysum = new Float_t[13];     /*                          */
    
    mw_x    = new Float_t[13];     /*                          */
    mw_y    = new Float_t[13];     /*                          */
    mw_wire = new Float_t[13];     /* special for Helmut       */
    
    z_x_s2 = 0;          
    z_y_s2 = 0;          
    z_x_s4 = 0;          
    z_y_s4 = 0;          
    z_x_s8 = 0;          
    z_y_s8 = 0;          
    
    b_mw_xsum = new Bool_t[13];   /*  wc on sum               */
    b_mw_ysum = new Bool_t[13];   /*                          */
    
    focx_s2 = 0;         /*  FRS foci                */  
    focy_s2 = 0;         /*                          */  
    focx_s4 = 0;         /*                          */  
    focy_s4 = 0;         /*                          */  
    focx_s8 = 0;         /*  FRS foci                */  
    focy_s8 = 0;         /*                          */
    angle_x_s2 = 0;      /*                          */
    angle_y_s2 = 0;      /*                          */
    angle_x_s4 = 0;      /*                          */
    angle_y_s4 = 0;      /*                          */
    angle_x_s8 = 0;      /*                          */
    angle_y_s8 = 0;      /*                          */
    
    tpc_x_s2_foc_21_22=-999.;
  tpc_y_s2_foc_21_22=-999.;
  tpc_x_s2_foc_23_24=-999.;
  tpc_y_s2_foc_23_24=-999.;
  tpc_x_s2_foc_22_24=-999.;
  tpc_y_s2_foc_22_24=-999.; 
  tpc_x_s4=-999.;
  tpc_y_s4=-999.;
  tpc_angle_x_s2_foc_21_22=-999.;
  tpc_angle_y_s2_foc_21_22=-999.;
  tpc_angle_x_s2_foc_23_24=-999.;
  tpc_angle_y_s2_foc_23_24=-999.;
  tpc_angle_x_s2_foc_22_24=-999.;
  tpc_angle_y_s2_foc_22_24=-999.; 
  tpc_angle_x_s4=-999.;
  tpc_angle_y_s4=-999.;
  
  tpc21_22_sc21_x = -999.;      /* SC21                     */
  tpc23_24_sc21_x = -999.;      /* SC21                     */
  tpc22_24_sc21_x = -999.;      /* SC21                     */
  tpc21_22_sc21_y = -999.;      /* SC21                     */
  tpc23_24_sc21_y = -999.;      /* SC21                     */
  tpc22_24_sc21_y = -999.;      /* SC21                     */
  
  tpc23_24_sc22_x = -999;
  tpc21_22_sc22_x = -999;
  
  tpc21_22_s2target_x = -999.;      /* S2TARGET                     */
  tpc23_24_s2target_x = -999.;      /* S2TARGET                     */
  tpc22_24_s2target_x = -999.;      /* S2TARGET                     */
  tpc21_22_s2target_y = -999.;      /* S2TARGET                     */
  tpc23_24_s2target_y = -999.;      /* S2TARGET                     */
  tpc22_24_s2target_y = -999.;      /* S2TARGET                     */
  
  tpc_sc41_x = -999.;      /* SC41                     */
  tpc_sc41_y = -999.;      /* SC41 Y                   */
  tpc_sc42_x = -999.;      /* SC42                     */
  tpc_sc42_y = -999.;      /* tracked SC42 Y pos       */
  tpc_sc43_x = -999.;      /* SC43                     */
  tpc_sc43_y = -999.;      /* SC43 Y                   */

  tpc_music41_x = -999.;      /* MUSIC41 x                     */
  tpc_music41_y = -999.;      /* tracked MUSIC41 Y pos       */
  tpc_music42_x = -999.;      /* MUSIC42 x                     */
  tpc_music42_y = -999.;      /* tracked MUSIC42 Y pos       */
  tpc_music43_x = -999.;      /* MUSIC43 x                     */
  tpc_music43_y = -999.;      /* tracked MUSIC43 Y pos       */
  tpc_s4target_x = -999.;     /* S4 target */
  tpc_s4target_y = -999.;     /* S4 target */

    
    // TPC part
    tpc_x    = new Float_t[7];
    tpc_y    = new Float_t[7];
    b_tpc_xy = new Bool_t[7];
    tpc_de = new Float_t[7];
    b_tpc_de = new Bool_t[7];

    tpc_csum = new Int_t*[7];
    b_tpc_csum = new Bool_t*[7];
    
    
    
    for(int i=0; i < 7; ++i){
    
    tpc_csum[i] = new Int_t[4];
    b_tpc_csum[i] = new Bool_t[4];
    
    }
    x0 = 0;
    x1 = 0;
    
    //TPCs 21 & 22 @ S2 focus
    tpc_x_s2_foc_21_22      = 0;
    tpc_y_s2_foc_21_22      = 0;
    tpc_angle_x_s2_foc_21_22    = 0;
    tpc_angle_y_s2_foc_21_22    = 0;
    
    //TPCs 23 & 24 @ S2 focus
    tpc_x_s2_foc_23_24      = 0;
    tpc_y_s2_foc_23_24      = 0;
    tpc_angle_x_s2_foc_23_24    = 0;
    tpc_angle_y_s2_foc_23_24    = 0;
    
    //TPCs 22 & 24 @ S2 focus
    tpc_x_s2_foc_22_24      = 0;
    tpc_y_s2_foc_22_24      = 0;
    tpc_angle_x_s2_foc_22_24    = 0;
    tpc_angle_y_s2_foc_22_24    = 0;
    
    //TPCs 3 & 4 @ S2 focus (old tracking)
    tpc_x_s2_foc        = 0;
    tpc_y_s2_foc        = 0;
    tpc_angle_x_s2_foc      = 0;
    tpc_angle_y_s2_foc      = 0;
    
    //TPCs 41 & 42 @ s4 focus
    tpc_x_s4            = 0;
    tpc_y_s4            = 0;
    tpc_angle_x_s4      = 0;
    tpc_angle_y_s4      = 0;
    
    //TPCs 5 & 6 @ S4 OTPC (exp s388)
    tpc_x_s4_target2        = 0;
    tpc_y_s4_target2        = 0;
    
    
    //TPCs 3 & 4 @ S2 first Si tracking detector (exp s388)
    tpc_x_s2_target1        = 0;
    tpc_y_s2_target1        = 0;
    //      Float_t tpc_angle_x_s4_target2;
    //      Float_t tpc_angle_y_s4_target2;
    
    
   // tpc_mw_sc21_x = 0;      /* SC21 x                   */
    tpc_sc41_x = 0;      /* SC41 x                   */
  //  tpc_sc21_y = 0;      /* SC21 y                   */
    tpc_sc41_y = 0;      /* SC41 y                   */
    
    
    mw_sc21_x = 0;          /* SC21                     */
    mw_sc22_x = 0;          /* SC22 (LAND finger)       */
    sc41_x = 0;          /* SC41                     */
    itag_stopper_x = 0;  /* ITAG Stopper             */
    sc42_x = 0;          /* SC42                     */
    sc43_x = 0;          /* SC43                     */
    sc21_y = 0;          /* tracked SC21 Y pos       */
    sc41_y = 0;          /* SC41 Y                   */
    sc42_y = 0;          /* SC42 Y                   */
    sc43_y = 0;          /* SC43 Y                   */
    sc81_x = 0;          /* SC81          positions  */
    
    music1_x1 = 0;       /* parameters to calibrate  */
    music1_x2 = 0;       /* parameters to calibrate  */
    music1_x3 = 0;       /* parameters to calibrate  */
    music1_x4 = 0;       /* parameters to calibrate  */
    music1_y1 = 0;       /* parameters to calibrate  */
    music1_y2 = 0;       /* parameters to calibrate  */
    music1_y3 = 0;       /* parameters to calibrate  */
    music1_y4 = 0;       /* parameters to calibrate  */
    music2_x  = 0;       /* MUSIC1,2 positions       */
    
    
    //Si detectors
    si_e1 = 0;
    si_e2 = 0;
    si_e3 = 0;
    si_e4 = 0;
    si_e5 = 0;
    
 for(int module=0; module<VFTX_N; module++){
    for(int channel=0; channel<VFTX_MAX_CHN; channel++) {
      for(int hit=0; hit<VFTX_MAX_HITS; hit++) { 
	vftx_cc[module][channel][hit] = 0;
	vftx_ft[module][channel][hit] = 0;
      }      
      vftx_mult[module][channel] = 0;
    }
 }
    
    
//     dssd_e = new Float_t[32]; //[3][2][16]
//     dssd_maxenergy = 0; //[3] 
//     dssd_maxindex  = 0;  //[3]
//     
//     dssd_xmaxenergy = 0;
//     dssd_xmaxindex  = 0;
//     dssd_ymaxenergy = 0;
//     dssd_ymaxindex  = 0;
    
    /******************************************************************/
    /*** ANALYSIS Parameters ***/
    
    /*** MUSIC Conditions ***/
    
        // MUSIC part
        music1_anodes_cnt = 0;
        music2_anodes_cnt = 0;
        music3_anodes_cnt = 0;
        
        music1_x_mean = -999.;
        music2_x_mean = -999.;
        music3_x_mean = -999.;
        
        music_b_e1 = new Bool_t[8];
        music_b_t1 = new Bool_t[8];
        music_b_e2 = new Bool_t[8];
        music_b_t2 = new Bool_t[8];
        music_b_e3 = new Bool_t[4];
        music_b_t3 = new Bool_t[4];
        
        for(int i = 0; i < 8; ++i){
        
        music_b_e1[i] = false;
        music_b_t1[i] = false;
        music_b_e2[i] = false;
        music_b_t2[i] = false;
        if(i < 4) music_b_e3[i] = false;
        if(i < 4) music_b_t3[i] = false;

        }
        
        b_de1 = false;
        b_de2 = false;
        b_de3 = false;
        
        de = new Float_t[3];
        de_cor = new Float_t[3];
        
        for(int i = 0; i < 3; ++i){

            de[i] = 0;
        de_cor[i] = 0;
        }
        
        b_dt3 = false;
        
        x1_mean = 0;
        
      //  b_decor = false;
    
    //SCI part
    sci_l = new Float_t[12];  
    sci_r = new Float_t[12];  
    sci_e = new Float_t[12];  
    sci_tx = new Float_t[12];  
    sci_x = new Float_t[12];
    itag_42_e = 0;
    itag_43_e = 0;
    sci_u5 = 0;
    sci_d5 = 0;
    sci_ty5 = 0;  
    sci_y5 = 0;  
    sci_tofll2 = 0;
    sci_tofrr2 = 0;
    sci_tof2 = 0;
    sci_tofll3 = 0;
    sci_tofrr3 = 0;
    sci_tof3 = 0;
    sci_tofll4 = 0;  
    sci_tofrr4 = 0;
    sci_tof4 = 0;
    sci_tof2 = 0; sci_tof3 = 0; sci_tof4 = 0;
    
    sci_veto_l = new Float_t[3];
    sci_veto_r = new Float_t[3];
    sci_veto_e = new Float_t[3];
    
    sci_b_l = new Bool_t[12];  
    sci_b_r = new Bool_t[12];  
    sci_b_e = new Bool_t[12];  
    sci_b_tx = new Bool_t[12];  
    sci_b_x = new Bool_t[12];  
    sci_b_u5 = false;
    sci_b_d5 = false;
    sci_b_ty5 = false;
    sci_b_y5 = false;
    sci_b_tofll2 = false;
    sci_b_tofrr2 = false;
    sci_b_tofll3 = false;
    sci_b_tofrr3 = false;
    sci_b_tofll4 = false;
    sci_b_tofrr4 = false;  
    sci_b_detof = false;
    sci_b_veto_l = new Bool_t[3];  
    sci_b_veto_r = new Bool_t[4];  
    sci_b_veto_e = new Bool_t[3];
    
    for(int i = 0; i < 12; ++i){
        
        sci_l[i] = 0;  
        sci_r[i] = 0;  
        sci_e[i] = 0;  
        sci_tx[i] = 0;  
        sci_x[i] = 0;
           
        sci_b_l[i] = false;  
        sci_b_r[i] = false;  
        sci_b_e[i] = false;  
        sci_b_tx[i] = false;  
        sci_b_x[i] = false;  
           
        if(i < 3){
        
        sci_veto_l[i] = 0;
        sci_veto_r[i] = 0;
        sci_veto_e[i] = 0;
            
        sci_b_veto_l[i] = false;  
        sci_b_veto_e[i] = false;
        
        }
        
        if(i < 4) sci_b_veto_r[i] = false;          
    }

    // ID part
    
    id_x2 = 999;
    id_y2 = 999;
    id_a2 = 999;
    id_b2 = 999;
    id_x4 = 999;
    id_y4 = 999;
    id_a4 = 999;
    id_b4 = 999;
    id_x8 = 999;
    id_y8 = 999;
    id_a8 = 999;
    id_b8 = 999;
    for(int i=0; i<7; ++i){
    tpc_x[i] = 999;
    tpc_y[i] = 999;
    }
    
    tpc_sc41_x = -999.;      /* SC41                     */
  tpc_sc41_y = -999.;      /* SC41 Y                   */
  tpc_sc42_x = -999.;      /* SC42                     */
  tpc_sc42_y = -999.;      /* tracked SC42 Y pos       */
  tpc_sc43_x = -999.;      /* SC43                     */
  tpc_sc43_y = -999.;      /* SC43 Y                   */

  tpc_music41_x = -999.;      /* MUSIC41 x                     */
  tpc_music41_y = -999.;      /* tracked MUSIC41 Y pos       */
  tpc_music42_x = -999.;      /* MUSIC42 x                     */
  tpc_music42_y = -999.;      /* tracked MUSIC42 Y pos       */
  tpc_music43_x = -999.;      /* MUSIC43 x                     */
  tpc_music43_y = -999.;      /* tracked MUSIC43 Y pos       */
  tpc_s4target_x = -999.;     /* S4 target */
  tpc_s4target_y = -999.;     /* S4 target */
  
   // MultiHitTDC
   mhtdc_tof8121 = -999.9;
   mhtdc_tof4121 = -999.9;
   mhtdc_tof4221 = -999.9;
   mhtdc_tof4321 = -999.9;
   mhtdc_tof3121 = -999.9;

   mhtdc_sc21lr_dt = -999.9;
   mhtdc_sc31lr_dt = -999.9;
   mhtdc_sc41lr_dt = -999.9;
   mhtdc_sc42lr_dt = -999.9;
   mhtdc_sc43lr_dt = -999.9;
   mhtdc_sc81lr_dt = -999.9;

   mhtdc_sc21lr_x = -999.9;
   mhtdc_sc31lr_x = -999.9;
   mhtdc_sc41lr_x = -999.9;
   mhtdc_sc42lr_x = -999.9;
   mhtdc_sc43lr_x = -999.9;
   mhtdc_sc81lr_x = -999.9;
  
    id_b_x2 = false;
    id_b_x4 = false;
   // id_b_x8 = false;
    id_b_detof2 = false;
    
    id_brho = new Float_t[2];      /* position-corr. BRho      */
    id_rho = new Float_t[2];       /* Position-corrected Rho   */
    
    id_brho[0] = 0;
    id_brho[1] = 0;
    id_rho[0] = 0;
    id_rho[1] = 0;
    
    id_beta = 0;        /* Beta from TOF            */
    id_beta3 = 0;        /* Beta from TOF            */
    id_gamma = 0;       /* Gamma from TOF           */
    id_AoQ = 0;
    id_AoQ_corr = 0;
    
    id_v_cor = 0;       /* Velocity correction  */
    id_v_cor2 = 0;      /* Velocity correction TUM 2 */
    id_v_cor3 = 0;      /* Velocity correction Old Music */
    id_z = 0;
    id_z2 = 0;
    id_z3 = 0;
    id_znocorr = 0;
    id_z2nocorr = 0;
    id_energy_geL = 0;
    id_tac_41_geL = 0;
    id_stopper_x = 0;
    id_energy_geL_raw = 0;
    id_tac_41_geL_raw = 0;
    id_trigger = 0;
    id_scal_seconds = 0;
    id_scal_geL = 0;
    id_scal_sc21 = 0;
    id_scal_sc41 = 0;
    id_scal_sc42 = 0;
    id_scal_sc43 = 0;
    id_scal_sc81 = 0;
    id_de_s2tpc =0.0;
    id_b_de_s2tpc = kFALSE;
    id_z_sc81 = 0.0;
    id_v_cor_sc81 = 0.0;
    id_b_z_sc81 = kFALSE;
    id_z_s2tpc =0.0;
    id_v_cor_s2tpc =0.0;
    id_b_z_s2tpc = kFALSE;
  
    id_b_AoQ = false;
    id_b_z = false;
    id_b_AoQ_s2s8    = kFALSE;
    
    id_b_z2 = false;
    id_b_z3 = false;
    //id_b_x2AoQ = false;
    
    id_b_x2AoQ = new Bool_t[5];
    id_b_x4AoQ = new Bool_t[5];
    id_b_x4AoQ_Z = new Bool_t[5]; 
    id_b_z_AoQ = new Bool_t[5];
    id_b_music_z = new Bool_t[5];
    
    id_mhtdc_beta_s2s8  = -999;
  id_mhtdc_gamma_s2s8 = -999;
  id_mhtdc_delta_s2s8 = -999;
  id_mhtdc_aoq_s2s8   = -999;
  id_mhtdc_z_s2tpc    = -999;
  id_mhtdc_zcor_s2tpc = -999;
  id_mhtdc_z_sc81     = -999;
  id_mhtdc_zcor_sc81  = -999;
  id_mhtdc_v_cor_s2tpc = 0.;
  id_mhtdc_v_cor_sc81 = 0.;
  
  id_mhtdc_beta_s2s4  = -999;
  id_mhtdc_gamma_s2s4 = -999;
  id_mhtdc_delta_s2s4 = -999;
  id_mhtdc_aoq_s2s4   = -999;
  id_mhtdc_z_music41  = -999;
  id_mhtdc_zcor_music41  =-999;
  id_mhtdc_v_cor_music41 =0;
  id_mhtdc_z_music42 =-999;
  id_mhtdc_zcor_music42 =-999;
  id_mhtdc_v_cor_music42 =0;
    
    for(int i = 0; i < 5; i++){
        id_b_x2AoQ[i] = false;
        id_b_x4AoQ[i] = false;
        id_b_x4AoQ_Z[i] = false; 
        id_b_z_AoQ[i] = false;
        id_b_music_z[i] = false;
        
    }
    
    firsttimestamp = 0;
    firstevent = false;
    
    ts = 0;  // relative time from start of the spill 
    ts2 = 0;  // relative time from start of the spill does not reset at end extraction
    
    time_in_ms =0;
    spill_count =0;
    ibin_for_s =0;
    ibin_for_100ms =0;
    ibin_for_spill =0;
    ibin_clean_for_s=0;
    ibin_clean_for_100ms=0;
    ibin_clean_for_spill=0;
    // MRTOF part :
    
    mrtof_tof = 0;
    
    mrtof_si_e1 = 0;
    mrtof_si_e2 = 0;
    mrtof_si_e3 = 0;
    mrtof_si_e4 = 0;
    mrtof_si_e5 = 0;
    
    Setup_Conditions();   
}

//---------------------------------------------------------------

FRS_Detector_System::~FRS_Detector_System(){
    

    
    for(int i = 0; i < 32; ++i){
    
   // delete[] vme2s[i];           // User TDC (V1290) 
    //delete[] vme2s_trailing[i];  // User TDC (V1290) 

    //delete[] leading_v1290_main[i];           // Mtof TDC (V1290) 
    //delete[] trailing_v1290_main[i];  // Mtof TDC (V1290) 

   // delete[] nhit5[i];       // multiplicity (V1290)
    //delete[] nhit_v1290_main[i];       // multiplicity (V1290)

    }
    /*delete[] vme2s;           // User TDC (V1290)   
    delete[] vme2s_trailing;  // User TDC (V1290) */  
        
    //delete[] leading_v1290_main;           // Mtof TDC (V1290)    
    //delete[] trailing_v1290_main;  // Mtof TDC (V1290)    
        
   // delete[] nhit5;       // multiplicity (V1290)   
   // delete[] nhit_v1290_main;       // multiplicity (V1290)
    
    delete[] vme2scaler;         // User Crate Messhure
    delete[] vme3scaler;         // User Crate Messhure

    delete[] firstTS;
    
    delete[] previousTimeStamp;
    
    /*******************************************************************/
    /***SORT STUFF***/
  
   // delete[] ts_word; //for the titris time stammping
    //delete[] tsys_word; //for the system time

    // scaler readings     
    //delete[] sc_long; //changed from 32 to 64 (10.07.2018)
    delete[] sc_long2;


    // part of MW parameter
    delete[]  mw_an;       /*  anode time              */
    delete[] mw_xl;       /*  Rohdaten                */
    delete[] mw_xr;       /*                          */
    delete[] mw_yu;       /*                          */ 
    delete[] mw_yd;       /*                          */

    delete[] tpc_timeref;
    // TPC part //(HAPOL-25/03/06) 6 TPCs each one with 2 delay lines each!!    
    //7 TPCs (4 at S2, 2 at S4 and one at S3) 03.07.2018 SB
    
    
    for(int i=0; i < 7; ++i){
    
    delete[] tpc_l[i];
    delete[] tpc_r[i];
    delete[] tpc_lt[i];
    delete[] tpc_rt[i];
    delete[] tpc_dt[i];
    delete[] tpc_a[i];
    
    }
     delete[] tpc_l;    
     delete[] tpc_r; 
     delete[] tpc_lt;    
     delete[] tpc_rt;    
     delete[] tpc_dt;    
     delete[] tpc_a;

      
    // User multihit TDC
    
    delete[] tdc_sc41l;
    delete[] tdc_sc41r;
    delete[] tdc_sc21l;
    delete[] tdc_sc21r;
    delete[] tdc_sc42l;
    delete[] tdc_sc42r;
    delete[] tdc_sc43l;
    delete[] tdc_sc43r;
    delete[] tdc_sc81l;
    delete[] tdc_sc81r;
    delete[] tdc_sc31l;
    delete[] tdc_sc31r;
    delete[] tdc_sc11;

    // MUSIC1 part
    delete[]  music_e1;     /* Raw energy signals       */
    delete[] music_t1;     /* Raw drift time           */
 //   delete[] music_pres;   /* Music Druck              */
  //  delete[]  music_temp;   /* Music Temperatur         */
    
    // MUSIC2 part
    delete[] music_e2;     /* Raw energy signals       */
    delete[] music_t2;     /* Raw drift time           */
    
    // MUSIC3 part (OLD MUSIC)
    delete[] music_e3;     /* Raw energy signals       */
    delete[] music_t3;     /* Raw drift times          */

    delete[] dssd_adc;
    
    id_AoQ=0;
    //Proc_iterator=0;

    /******************************************************************/
    /**CALIBRATION Parameters**/
    
    // MON data declarations
    
    delete[] check_first_event;
    delete[] scaler_time_count; 
    delete[] scaler_spill_count; //UInt_t
    delete[] scaler_time_check_last;//UInt_t
    delete[] scaler_spill_check_last;//UInt_t 
    delete[] check_increase_time;//UInt_t 
    delete[] check_increase_spill;//UInt_t
    delete[] scaler_increase_event;//UInt_t
    delete[] scaler_last_event;
    
    
    delete[] scaler_save;
    
    delete[] mon_inc;
    
    delete[] coin;
    
    // MW part
   delete[]  mw_xsum;     /*                          */
   delete[]  mw_ysum;     /*                          */
    
   delete[]  mw_x;        /*                          */
   delete[]  mw_y;        /*                          */
   delete[]  mw_wire;     /* special for Helmut       */
    
   delete[]  b_mw_xsum;   /*  wc on sum               */
   delete[]  b_mw_ysum;   /*                          */
    
    // TPC part
    delete[] tpc_x;
    delete[] tpc_y;
    delete[] b_tpc_xy;
    delete[] tpc_de;
    delete[] b_tpc_de;

    for(int i=0; i < 7; ++i){
    
    delete[] tpc_csum[i];
    delete[] b_tpc_csum[i];
    
    }
    delete[] tpc_csum;  
    delete[] b_tpc_csum;
  
//    delete[] dssd_e; //[3][2][16]
    
    /******************************************************************/
    /*** ANALYSIS Parameters ***/
    
    /*** MUSIC Conditions ***/
    
    // MUSIC part
    
    delete[] music_b_e1;
    delete[] music_b_t1;
    delete[] music_b_e2;
    delete[] music_b_t2;
    delete[] music_b_e3;
    delete[] music_b_t3;
    
    delete[] de;
    delete[] de_cor;
    
    //SCI part
    delete[] sci_l;  
    delete[] sci_r;  
    delete[]  sci_e;  
    delete[]  sci_tx;  
    delete[] sci_x;
    
    delete[]  sci_veto_l;
    delete[]  sci_veto_r;
    delete[]  sci_veto_e;
    
    delete[]  sci_b_l;  
    delete[]  sci_b_r;  
    delete[]  sci_b_e;  
    delete[]  sci_b_tx;  
    delete[]  sci_b_x;  
    delete[]  sci_b_veto_l;  
    delete[]  sci_b_veto_r;  
    delete[]  sci_b_veto_e;
    
    // ID part
    
    delete[] id_brho;      /* position-corr. BRho      */
    delete[] id_rho;       /* Position-corrected Rho   */
    
    delete[] id_b_x2AoQ;
    delete[] id_b_x4AoQ;
    
    delete[] id_b_x4AoQ_Z; 
    delete[] id_b_z_AoQ;
    delete[] id_b_music_z;
    
}
// void FRS_Detector_System::load_VFTX_Calibration_Files(){
//   
//     //Load all wired Calibration files specified by MAP
//     char filename[1000];
//     ifstream file;
// 
//     int b_iter = 0;
//     double bin,val;
// 
//     //Cal_arr = new double**[10];
// 
//     const char* format = "%lf %lf";
// 
//     bool first = true;
//       string line;
//     for(int i = 0;i < 1;++i){
//        for(int j = 0;j < 32;++j){
//      
//         b_iter = 0;
// 
// 
//         sprintf(filename,"Configuration_Files/VFTX_Calib/VFTX_%05d_Bin2Ps_ch%02d.dat",i,j);
//         file.open(filename);
//         
//         if(file.fail()){
//             cerr << "Could not find VFTX Calibration file " << i << " " << j << endl;
//             exit(0);
//         }
//         while(file.good()){
//          //   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
//             getline(file,line,'\n');
//             if(line[0] == '#') continue;
//             sscanf(line.c_str(),format,&bin,&val);
//             if(first) bins_x_arr[b_iter] = bin;
//             VFTX_Cal_arr[i][j][b_iter] = val;
//             b_iter++;
//            
//         }
//  
//         first = false;
// 
//         file.close();
//         file.clear();
//         }
//     }
// }
//---------------------------------------------------------------

void FRS_Detector_System::Process_MBS(TGo4MbsSubEvent* psubevt){
    
 skip=false;
     // cout<<"frs->dist_focS4 !!!!!!!!" << frs->dist_focS4<< endl;
   // cout<<"frs->bfield[0] " << frs->bfield[0]<< " frs->bfield[1]  " << frs->bfield[1] <<" frs->bfield[2] " << frs->bfield[2]<< " frs->bfield[3]  " << frs->bfield[3]<< endl;
   // Setup_Parameters();

    FRS_Unpack(psubevt);
   
   // if(skip==false){
    FRS_Sort();
    FRS_Calib();
    FRS_Anal();
    //}
    
}

//-----------------------------------------------------------------//
//--------------------- UNPACKER STUFF ----------------------------//
//-----------------------------------------------------------------//

void FRS_Detector_System::FRS_Unpack(TGo4MbsSubEvent* psubevent){
 
  psubevt = psubevent;
  pdata = psubevt->GetDataField();
   len = 0;      
   //vme_chn;
   lenMax = (psubevt->GetDlen()-2)/2;
   Int_t event_flag = *pdata++;  len++; // 0x200 or 0x000 etc. 20200123MERGE
  
  if( (psubevt->GetType() == 36) && (psubevt->GetSubtype() == 3600) ) //tpat
   {
     Int_t *pdata = psubevt->GetDataField();
     int tpat_len = 0xFF & (*pdata) ;
     pdata++;
     for(int ii=0; (ii<tpat_len && ii<32) ; ii++){
       //tgt->tpat_main[ii] = *pdata;
       //  printf("write tpat data %d = 0x%08x\n",ii,tgt->tpat_main[ii]);
       pdata++;
     }
     return;
   }
  switch(psubevt->GetProcid())
    {

      
    //===========
    case 30: ///FRS Crate
    
   
    { // start of v830 (do not remove this bracket)
            Int_t vme_geo = getbits(*pdata,2,12,5);
            Int_t vme_type = getbits(*pdata,2,9,3);
            Int_t vme_nlw =  getbits(*pdata,2,3,6);
            pdata++; len++;
            if(vme_type!=4){   std::cout<<"E> Scaler type missed match ! GEO"<<vme_geo<<" "<<" type 4 =/="<<vme_type<<std::endl; }
            for(int i_ch=0; i_ch<vme_nlw; i_ch++){
              scaler_frs[i_ch] = *pdata;
             // cout<<"11scaler_frs[i_ch] " <<scaler_frs[i_ch] <<" i_ch " <<i_ch<< endl;
              ////               printf("scaler_frs[ch=%d] = %d\n",i_ch,*pdata);
              pdata++; len++;
            }
        } //end of V830
        
        //----- v775 TDC/V785 ADC -----//
        while (len < lenMax){          
          /* read the header longword and extract slot, type & length  */
          Int_t vme_chn = 0; 
          Int_t vme_geo = getbits(*pdata,2,12,5);
          Int_t vme_type = getbits(*pdata,2,9,3);
          Int_t vme_nlw =  getbits(*pdata,1,9,6); 
    
    
    
          pdata++; len++;
     
          /* read the data */
          if ((vme_type == 2) && (vme_nlw > 0)) {
            for(int i=0;i<vme_nlw;i++) {  
              vme_geo = getbits(*pdata,2,12,5);
              vme_type = getbits(*pdata,2,9,3);
              vme_chn = getbits(*pdata,2,1,5);
              vme_frs[vme_geo][vme_chn] = getbits(*pdata,1,1,16);
              ////        printf("vme_frs[geo=%d][ch=%d] = %d\n",vme_geo,vme_chn,getbits(*pdata,1,1,16));
              pdata++; len++;
            }
          }
          // read and ignore the expected "end-of-block" longword   
           pdata++; len++;
           }  /* end of the while... loop  */
           break;
      
      
    //==========  
    case 10:  /// Main crate
        
      if((psubevt->GetType() == 10) && (psubevt->GetSubtype() == 1) && lenMax < 15) // in this case, time stamp data.
      {
    int ii =0;
    while(len < lenMax){
       timestamp_main[ii] = *pdata;
      pdata++; len++; ii++;
    }
    break;
      }
      
      ///1st babababa
      if(*pdata != 0xbabababa){ 
         // std::cout<<"1E> ProcID 10 : Barrier-1 missed !" << std::endl;
          //break;
          skip=true;
    }
  
      pdata++; len++;
      //----------
//       /// there are 2 words (time stamp) but these are redundant
//       /// skip 2 words;
//         pdata++; len++;
//         pdata++; len++;
// 
// 
//    /// additional babababa request by JoseLuis
//     if(*pdata == 0xbabababa){
//       pdata++; len++;
//     }
//          else{ 
//               //std::cout<<"E> ProcID 10 : Barrier-2 missed !" << *pdata  << std::endl; 
//             //  break;
//         }
//           pdata++; len++;
    

     // if(ModSetup->Nb_TDC>0)
    if (len < lenMax){
      ///v1290 TDC              
      Int_t vme_geo = getbits(*pdata,1,1,5);
      Int_t vme_type = getbits(*pdata,2,12,5);
      pdata++; len++;
      Int_t multihit = 0;//, counter = 0;
          for(int i=0;i<32;i++){
            for(int j=0;j<2;j++){
                nhit_v1290_main[i][j]=0;
     }
 }
     // if(vme_geo==0) break;
//       const auto it_Module = it_Crate->second.find(vme_geo);
//       int IdMod = it_Module->second;
//       if(it_Module == it_Crate->second.end())
//       std::cout<<"E> Crate Mapping does not have this module (vmeGEO**) "<<vme_geo<<" in Crate :"<<psubevt->GetProcid()<<std::endl;

#ifdef DEBUG
      std::cout<<"mTDC geo = "<<vme_geo<<", type = "<<vme_type<<std::endl;
#endif  

      if (vme_type == 8)
        {
          while (len < lenMax) 
        {
#ifdef DEBUG
          std::cout<<"word : "<<std::bitset<32>(*pdata)<<" ";
#endif
          vme_type = getbits(*pdata,2,12,5);
          if(vme_type==1) // headerTDC
            { 
              pdata++; len++;
            }
          //multihit = 0;
#ifdef DEBUG
          std::cout<<"reading :"<<vme_type<<std::endl;
#endif
          vme_type = getbits(*pdata,2,12,5);
#ifdef DEBUG
          std::cout<<"word : "<<std::bitset<32>(*pdata)<<" type:"<<vme_type<<" ";
#endif
          if(vme_type == 0)
            {
              /// this indicates a TDC measurement
              Int_t vme_chn = getbits(*pdata,2,6,5);
              Int_t LeadingOrTrailing = getbits(*pdata,2,11,1);
              Int_t value = getbits(*pdata,1,1,21);

              multihit =  nhit_v1290_main[vme_chn][LeadingOrTrailing];
               // std::cout << "tdc vme_chn = " << vme_chn<<" multihit: " << multihit << " value: " << value << " LeadingOrTrailing: " << LeadingOrTrailing <<endl;
#ifdef DEBUG
              std::cout << "     tdc vme_chn = " << vme_chn;
              std::cout << " multihit: " << multihit << " ";
#endif

              if (multihit >= 10)
            {
              pdata++; len++;
              continue;
            }
              if(LeadingOrTrailing == 0)
            {
#ifdef DEBUG
              std::cout << " +-> tdc L value = " << value << std::endl;
#endif
              
              if (value > 0)
                {
                   leading_v1290_main[vme_chn][multihit] = value;
           //     cout<<"1111leading_v1290_main[3][multihit] "<<leading_v1290_main[3][multihit]<<" psubevt->GetProcid() " <<psubevt->GetProcid() <<" psubevt->GetSubtype() " <<psubevt->GetSubtype() <<  endl;
                }
            }
          
              else
            {
#ifdef DEBUG
              std::cout << " +-> tdc T value = " << value << std::endl;
#endif
              if (value > 0)
                 trailing_v1290_main[vme_chn][multihit] = value;
            }
               nhit_v1290_main[vme_chn][LeadingOrTrailing]++;
        //cout<<"nhit_v1290_main " << nhit_v1290_main[vme_chn][LeadingOrTrailing] << " chn " << vme_chn << " LeadingOrTrailing "<<LeadingOrTrailing<< endl;;
              pdata++; len++;
            }
          else
            {
              // TDC trailer vme_type == 3 
#ifdef DEBUG
              std::cout<<"\n";
#endif
              if(vme_type != 3 && vme_type !=16)
            std::cout<<"E> MTDC strange type :"<<vme_type<<std::endl;
              if(vme_type==16)
            {
              Int_t vme_geoEnd = getbits(*pdata,1,1,5);
              if(vme_geo!=vme_geoEnd)
                std::cout<<"E> MTDC strange end buffer header :"<<vme_type<<" "<<vme_geo<<" != "<<vme_geoEnd<<std::endl;

              pdata++; len++;
              break;
            }
              pdata++; len++;
            }
           } ///while len loop
         }///if type==8
       }///end if len more due to strange logic // end of V1290
    
      

  if(*pdata != 0xbabababa) {
    //std::cout<<"E> ProcID 10 : Barrier-3 missed !" << std::endl;
  
     // break;
}
      pdata++; len++;


      
   
      //----- Next is V830 Scaler-----
     // if(ModSetup->Nb_Scaler > 0)
      {// start of v830 (do not remove this bracket)
      /** \note FRS SCALER module data (1 longword per channel)   
       *  This module has sequential readout therefore no channel
       *  number contained in the data longwords. 
       */
      // read the header longword and extract slot, type & length  
      Int_t vme_geo = getbits(*pdata,2,12,5);
      Int_t vme_type = getbits(*pdata,2,9,3);
      Int_t vme_nlw =  getbits(*pdata,2,3,6);

          pdata++; len++;
if(vme_type!=4){  
   // std::cout<<"E> Scaler type missed match ! GEO"<<vme_geo<<" "<<" type 4 =/="<<vme_type<<std::endl; 
    
}
        for(int i_ch=0; i_ch<vme_nlw; i_ch++){
          scaler_main[i_ch] = *pdata;
          //////         printf("scaler_main[ch=%d] = %d\n",i_ch,*pdata);
          pdata++; len++;
        }
          } ///end of V830

     ///--------the rest of the unpacking...V792--------    

        
      while (len < (psubevt->GetDlen()-2)/2)
    {
             
      Int_t vme_chn = 0;
      Int_t vme_geo = getbits(*pdata,2,12,5);
      Int_t vme_type = getbits(*pdata,2,9,3);
      Int_t vme_nlw =  getbits(*pdata,1,9,6); 
      
      pdata++; len++;
      #ifdef DEBUG
          std::cout<<"data "<<vme_geo<<" "<<vme_type<<" "<<vme_nlw<<" idmod:"<<IdMod<<std::endl;
    #endif
     
      
      /// read the data
     if(vme_type == 6){
    
              // not valid data !
          }
         
          if ((vme_type == 2) && (vme_nlw > 0)){
            for(int i=0;i<vme_nlw;i++){  
              vme_geo = getbits(*pdata,2,12,5);
              vme_type = getbits(*pdata,2,9,3);
              vme_chn = getbits(*pdata,2,1,5);
              vme_main[vme_geo][vme_chn] = getbits(*pdata,1,1,16);
            
               pdata++; len++;
            }
             /// read and ignore the expected "end-of-block" longword 
            pdata++; len++;
          }
        }  /// end of the while... loop  
    
       
          break; /// proc ID 10

      ///====================================       
    case 20:
    case 25:
 
      ///=============================== (Added on 07.07.2018)
      /* v775 TDC/V785 ADC */
      
      while (len < lenMax){
    
    /* read the header longword and extract slot, type & length  */
    Int_t vme_chn = 0;
    Int_t vme_geo = getbits(*pdata,2,12,5);
    Int_t vme_type = getbits(*pdata,2,9,3);
    Int_t vme_nlw =  getbits(*pdata,1,9,6); 
    // printf("pdata=0x%08x\n",*pdata); fflush(stdout);
    //std::cout<< vme_geo << std::endl;
    pdata++; len++;
    
    /* read the data */
    if ((vme_type == 2) && (vme_nlw > 0)) {
      for(int i=0;i<vme_nlw;i++) {  
        vme_geo = getbits(*pdata,2,12,5);
        vme_type = getbits(*pdata,2,9,3);
        vme_chn = getbits(*pdata,2,1,5);
        //printf("pdata=0x%08x\n",*pdata); fflush(stdout);
        if(20 == psubevt->GetProcid()||25 == psubevt->GetProcid())  vme_tpc[vme_geo][vme_chn] = getbits(*pdata,1,1,16);
       // if(25 == psubevt->GetProcid())  vme_tpc[vme_geo][vme_chn] = getbits(*pdata,1,1,16);
        pdata++; len++;
      }   
    }

    /* read and ignore the expected "end-of-block" longword */
    pdata++; len++;
      }  /* end of the while... loop  */
      
       ///===============================      
      ///All the ADCs are read with the code just above..(07.07.2018)
        
      /* rest of the unpacking */
      while (len < (psubevt->GetDlen()-2)/2)
    {
      
      /* read the header longword and extract slot, type & length  */
      Int_t vme_chn = 0;
      Int_t vme_geo = getbits(*pdata,2,12,5);
      Int_t vme_type = getbits(*pdata,2,9,3);
      Int_t vme_nlw =  getbits(*pdata,1,9,6);
      pdata++; len++;
#ifdef DEBUG
      std::cout<<"data "<<vme_geo<<" "<<vme_type<<" "<<vme_nlw;
#endif      
      
    }  /// end of the while... loop  
      break;
      ///===============================================================================
      ///================================================================ =============    
      
      case 35: // --- vftx at S2 SOFIA --- (2020/Mar/29)
            uint32_t p32_tmp, marker;
            uint16_t cnt, channel;
            int      word;
      if((psubevt->GetType() == 12) && (psubevt->GetSubtype() == 1)) // in this case, time stamp data.
	  {
           for(int module=0; module<VFTX_N; module++)
            for(int channel=0; channel<VFTX_MAX_CHN; channel++) {
             for(int hit=0; hit<VFTX_MAX_HITS; hit++) { 
	       vftx_cc[module][channel][hit] = 0;
	       vftx_ft[module][channel][hit] = 0;
             }      
             vftx_mult[module][channel] = 0;
            }
           pdata++;
           for(int module=0; module<VFTX_N; module++){      


  // first 32-bit word is a marker
  p32_tmp = *pdata++;
  marker  = (uint32_t)(p32_tmp & 0xff000000);
  cnt     = (uint16_t)((p32_tmp & 0x0003fe00)>>9);

 // printf("1: 0x%x    %d \n",p32_tmp, cnt);

  if(marker==0xab000000) {
      
    // second 32-bit word is the trigger window
    p32_tmp = *pdata++;
    //printf("2: 0x%x\n",p32_tmp);
    // third 32-bit word is a header we don't care of
    p32_tmp = *pdata++;
  //printf("3: 0x%x\n",p32_tmp);
    // DATA BLOCK
    for(word=0; word<cnt;word++) {
      p32_tmp = *pdata++;
      channel = (uint16_t)((p32_tmp & 0x3e000000)>>25);
      
     // printf("3: 0x%x  %d  \n",p32_tmp,channel);
      //printf("error in UNPACK TDC/VFTX_%02d,  %02d: %02d \n",module,channel, vftx_mult[module][channel]);
      
      if(vftx_mult[module][channel]<VFTX_MAX_HITS){
          ///Coarse clock
        vftx_cc[module][channel][vftx_mult[module][channel]] = (uint16_t)((p32_tmp & 0x00fff800)>>11);
        ///Fine time
	    vftx_ft[module][channel][vftx_mult[module][channel]] = (uint16_t)(p32_tmp & 0x07ff);
        h1_vftx_ft[module][channel]->Fill(vftx_ft[module][channel][0]);
        
//  cout<<"vftx_ft " <<vftx_ft[module][channel][vftx_mult[module][channel]] << " CC " << vftx_cc[module][channel][vftx_mult[module][channel]]<< " channel " << channel<< " vftx_mult[module][channel] " <<vftx_mult[module][channel] << endl;
      }
      vftx_mult[module][channel]++;
    }
   
     #if CALIBRATION_VFTX
      for(channel=0; channel<VFTX_MAX_CHN; channel++) {
      Vftx_Stat[module][channel]++;
      cout<<"Vftx_Stat[module][channel] " << Vftx_Stat[module][channel]<<" module " << module << " channel " << channel <<"\r"<< flush;
      if(Vftx_Stat[module][channel]==25000) {
          VFTX_Calibration(module,channel);
        }
      }
#endif
  } 
//   else {
//     printf("error in UNPACK TDC/VFTX_%02d marker: 0x%x\n",module,p32_tmp);
//   }
	    
//	    break;
        }
      }

      break; 
      
      
   ///===============================================================================
      ///================================================================ =============    
      
      
      
         case 40: /// Travel music crate (2020/Jan/23)
           
           { ///-----MDPP module----- (do not remove this bracket)
 // header
         Int_t header = *pdata;
         Int_t nword_mdpp = (0x3FF & header);
         pdata++;
           
         // main data (data or time information or end counter)
         for(int ii=0; ii<nword_mdpp; ii++){
           int tmp_data = *pdata;
           if( 1 == (0xF & (tmp_data>>28))   ){
             //      printf("real data 0x%08x\n",tmp_data);
             int ch = 0x3F & (tmp_data >> 16);
             if(0<=ch && ch<=15){
               int adc_data = 0xFFFF & tmp_data;
               //    printf("ADC data !!! ch = %2d, ADC = %d \n",ch, adc_data);
               if( (vme_trmu_adc[ch]) <= 0 ){ //first-come-first-served, for detailed multi-hit analysis, investigation needed.
              vme_trmu_adc[ch] = adc_data;
             //// printf(" vme_trmu_adc[%d] = 0x%08x; \n",ch,adc_data);
               }
             }else if(16<=ch && ch<=31){
               int tdc_data = 0xFFFF & tmp_data;
               //printf("TDC data !!! ch = %2d, TDC = %d \n",ch-16, tdc_data);
               if( ( vme_trmu_tdc[ch-16]) <= 0 ){//first-come-first-served, for detailed multi-hit analysis, investigation needed.
              vme_trmu_tdc[ch-16] = tdc_data;
             //////          printf("vme_trmu_tdc[%d] = 0x%08x; \n",ch-16,tdc_data);
               }
             }         
           }else if( 2 == (0xF & (tmp_data>>28))){
             // printf("ext time stamp 0x%08x\n",tmp_data);
           }else if(0x0 == tmp_data ){
             // printf("dummy 0x%08x\n",tmp_data);
           }else if( 3 == (0x3 & (tmp_data>>30))){
             // printf("end counter 0x%08x\n",tmp_data);
           }else{
             // printf("unknown data0x%08x\n",tmp_data);
           } 
           pdata++;
         }
           }//---end of MDPP module ---
           
           break;
         
         //=========
         default :
           break;
           
      } // end switch prodID
    
//      return kTRUE;
    
    }     



//-----------------------------------------------------------------//
//----------------------- SORT STUFF ------------------------------//
//-----------------------------------------------------------------//

void FRS_Detector_System::FRS_Sort(){
    ///CLEAR arrays to avoid multicounting for each proc ID 
//         if(psubevt->GetProcid() == 25){
//             
//         }
    tdc_sc41l[0] = 0;
    tdc_sc41r[0] = 0;
    tdc_sc21l[0] = 0;
    tdc_sc21r[0] = 0;
    tdc_sc42l[0] = 0;
    tdc_sc42r[0] = 0;
    tdc_sc43l[0] = 0;
    tdc_sc43r[0] = 0;
    tdc_sc81l[0] = 0;
    tdc_sc81r[0] = 0;
    tdc_sc31l[0] = 0;
    tdc_sc31r[0] = 0;
    tdc_sc11[0]  = 0;
      ///Clear for main crate 
      if(psubevt->GetProcid() == 10){   
          for (int n=0; n<32; n++){
            scaler_frs[n]=0;   
           }
           for(int o =0; o<64; o++){
            sc_long[o]=0;   
            }
           if(skip==true){
               tdc_sc41l[0] = 0;
                tdc_sc41r[0] = 0;
                tdc_sc21l[0] = 0;
                tdc_sc21r[0] = 0;
                tdc_sc42l[0] = 0;
                tdc_sc42r[0] = 0;
                tdc_sc43l[0] = 0;
                tdc_sc43r[0] = 0;
                tdc_sc81l[0] = 0;
                tdc_sc81r[0] = 0;
                tdc_sc31l[0] = 0;
                tdc_sc31r[0] = 0;
                tdc_sc11[0]  = 0;
            for(int i=0; i<21; i++){ 
                for (int j=0; j<32;j++){ 
//                    //vme_tpc[i][j]=0; 
//                   // vme_frs[i][j]=0; 
                   vme_main[i][j]=0; 
                  
                   
                }    
           }
           for(int i=0; i<32; i++){ 
                for (int j=0; j<10;j++){ 
                   leading_v1290_main[i][j]=0;
                   
                }    
           }
           for(int g=0; g<4; g++){
                    timestamp_main[g]=0;
                    ts_word[g]=0;
            }
           
           
           }///end skip
           
           for(int i=0; i<21; i++){ 
               for (int j=0; j<32;j++){ 
                   vme_tpc[i][j]=0; 
                   vme_frs[i][j]=0; 
                  // vme_main[i][j]=0; 
               }    
           }   
        id_x2=999;
        id_x4 =999;        
        id_AoQ=0;   
        id_AoQ_corr=0;     
 
       for(int i=0; i<7; i++){ 
             for(int j=0; j<2; j++){    
           tpc_a[i][j]=0;   
           tpc_l[i][j]=0;   
           tpc_r[i][j]=0;   
           tpc_lt[i][j]=0; 
         //  cout<<"CLEAR " <<" tpc_lt[i][j] " <<  tpc_lt[i][j] <<" i " << i  << " j " << j << endl;
           tpc_rt[i][j]=0;  
             }  
            for(int k=0; k<4; k++){     
                tpc_a[i][k]=0;  
                tpc_dt[i][k]=0; 
                tpc_csum[i][k] =0;  
                
                }   
            }   
           /* for(int i=0; i<10;i++){ 
        tdc_sc41l[i]= 0;    
        tdc_sc41r[i] = 0;   
        tdc_sc21l[i] = 0;   
        tdc_sc21r[i] = 0;   
        tdc_sc42l[i] = 0;   
        tdc_sc42r[i] = 0;   
        tdc_sc43l[i]= 0;    
        tdc_sc43r[i]= 0;    
        tdc_sc81l[i] = 0;   
        tdc_sc81r[i] = 0;   
            }  */ 
                
             for(int i=0;i<8;i++){  
            
        music_e1[i] =0; 
        music_e2[i] =0; 
         //printf("anode value= %d\n",music_e1[i]); 
        //music_e2[i] = vme_frs[12][24+i] & 0xfff;  
            
        music_t1[i] =0; 
        music_t2[i] =0; 
            
                }   
      
        }
        
        //Clear for case 30
        if(psubevt->GetProcid()==30){
//            
           for(int i=0; i<21; i++){ 
               for (int j=0; j<32;j++){ 
                  
                   vme_main[i][j]=0; 
               }    
           }   
           
          dt_21l_21r =0;  
        dt_41l_41r =0;  
        dt_21l_41l =0;  
     //   std::cout<<"vme_frs[12][5] & 0xfff " << dt_21l_41l <<std::endl;   
        dt_21r_41r =0;  
        dt_42l_42r =0;  
        dt_42l_21l =0;  
        dt_42r_21r =0;  
        dt_43l_43r =0;  
        dt_81l_81r =0;  
        dt_21l_81l =0;  
        dt_21r_81r =0; 
        dt_22l_22r  = 0;      /*                          */ 
        dt_22l_41l  = 0;
        dt_22r_41r  = 0;
        
        
        }
   if(psubevt->GetProcid()==20||psubevt->GetProcid()==25||psubevt->GetProcid()==30){
  tdc_sc41l[0] = 0;
    tdc_sc41r[0] = 0;
    tdc_sc21l[0] = 0;
    tdc_sc21r[0] = 0;
    tdc_sc42l[0] = 0;
    tdc_sc42r[0] = 0;
    tdc_sc43l[0] = 0;
    tdc_sc43r[0] = 0;
    tdc_sc81l[0] = 0;
    tdc_sc81r[0] = 0;
    tdc_sc31l[0] = 0;
    tdc_sc31r[0] = 0;
    tdc_sc11[0]  = 0;
       for(int i=0; i<32; i++){
            for(int j=0; j<10; j++){
       leading_v1290_main[i][j]=0;
       
       }
}       
       
       for(int i=0; i<32; i++){
            for(int j=0; j<10; j++){
       leading_v1290_main[i][j]=0;
       
       }
}

    de_21l = 0;
    de_21r = 0;
 
    de_41l = 0;
    de_41r = 0;
    de_42l = 0;
    de_42r =0;
    de_43l =0;
    de_43r = 0;
  
    de_81l = 0;
    de_81r = 0;
    de_31l = 0; 
    de_31r = 0;
  //  if(psubevt->GetSubtype() ==1){
    TRaw_vftx_21l=0;
    TRaw_vftx_21r=0;
    TRaw_vftx_22l=0;
    TRaw_vftx_22r=0;
    TRaw_vftx_41l=0;
    TRaw_vftx_41r=0;
    TRaw_vftx_42l=0;
    TRaw_vftx_42r=0;
//     TRaw_vftx_sofia_l=0;
//     TRaw_vftx_sofia_r=0;
   // }
    
    
        }
         if(psubevt->GetSubtype() !=1){
    TRaw_vftx_21l=0;
    TRaw_vftx_21r=0;
    TRaw_vftx_41l=0;
    TRaw_vftx_41r=0;
//     TRaw_vftx_sofia_l=0;
//     TRaw_vftx_sofia_r=0;
    for(int i=0; i<32; i++){
       TRaw_vftx[i]=0; 
        
    }
    }
    
      if(psubevt->GetProcid()==20||psubevt->GetProcid()==25){
       for (int n=0; n<32; n++){
            scaler_frs[n]=0;   
           }
           for(int o =0; o<64; o++){
            sc_long[o]=0;   
            }
      }

//     if(EventFlag == 0x100){
//     for(int i = 0; i < 4; ++i) ts_word[i] = vme_frs[20][i];
//     }
//     else if(EventFlag == 0x200){
//     for(int i = 0; i < 4; ++i) ts_word[i] = vme[3][20][i];
//     }
//     else if(EventFlag == 0x300){
//     for(int i = 0; i < 4; ++i) ts_word[i] = vme[4][20][i];
//     }

  ts_word[0] = timestamp_main[0];
  ts_word[1] = timestamp_main[1];
  ts_word[2] = timestamp_main[2];
  ts_word[3] = timestamp_main[3];

    timestamp = Long64_t(1)*ts_word[0] + Long64_t(0x10000)*ts_word[1] + Long64_t(0x100000000)*ts_word[2] + Long64_t(0x1000000000000)*ts_word[3];
    //cout<<"timestamp " <<timestamp << endl;
    // printf("qtrigger=%d timestamp=%ld \n",qtrigger,timestamp);
    
    tsys_word[0] = Long64_t(1)*ts_word[0] + Long64_t(0x10000)*ts_word[1] ; //s time low word
    tsys_word[1] = Long64_t(0x100000000)*ts_word[2] + Long64_t(0x1000000000000)*ts_word[3] ; //s time high worid... we do not use it

    if(PreviousTS < 0) tsys_word[2] = 0;
 
    else tsys_word[2] = (timestamp - PreviousTS)*1.e-5 ; //ms time since the previous s time (ftime routine)

    systemtime_s = tsys_word[2]*1e-3; //tsys_word[0] ; 
    systemtime_ms= tsys_word[2]; 


    if (qtrigger == 12){
    
    StartOfSpilTime = timestamp;
      
    StartOfSpilTime2 = timestamp;
   
    }   
    else if (qtrigger == 13){
      
    StartOfSpilTime = -1;
    }
    // rest are interesting only if trigger == 1
  
    // calculate time from spill start in sec
    if (StartOfSpilTime >= 0){
     
    timespill = (timestamp - StartOfSpilTime) * 1e-2;// microsec // 50000000.;
    }
  
    timespill2 = (timestamp - StartOfSpilTime2) * 1e-2; //microsec  // 50000000.;

    /* ### Pattern ###*/
    pattern = vme_frs[5][0];
    trigger = qtrigger;
 

  /* ### scalers:  */
  /* these are treated as 32-bit integers directly  */
    if (psubevt->GetProcid() == 30 && psubevt->GetType() == 12)
    {
      for(int i = 0 ;i < 32; ++i){
     
        sc_long[i]  = scaler_frs[i]; //frs crate
        sc_long[32+i] = scaler_main[i]; //main crate
     }
    }

  /* ### TA Ionization Chamber dE:  */
  
    
  /* ### MW anodes:  */
  for(int i = 0; i < 13; ++i)
    mw_an[i] = vme_frs[8][i] & 0xfff;
    
  /* ### MW cathodes:  */
  // from MW11 -> MW31
  for(int i = 0; i < 4; ++i)
    {  
      mw_xr[i] = vme_frs[8][16+i*4] & 0xfff;
      mw_xl[i] = vme_frs[8][17+i*4] & 0xfff;
      mw_yu[i] = vme_frs[8][18+i*4] & 0xfff;
      mw_yd[i] = vme_frs[8][19+i*4] & 0xfff;
    }
  // from MW41 -> MW51
  for(int i = 0 ;i < 3; ++i)
    {
      mw_xr[i+4] = vme_frs[9][0+i*4] & 0xfff;
      mw_xl[i+4] = vme_frs[9][1+i*4] & 0xfff;
      mw_yu[i+4] = vme_frs[9][2+i*4] & 0xfff;
      mw_yd[i+4] = vme_frs[9][3+i*4] & 0xfff;
    }
  // skip MW61
  //putting MW61 all to zero (09.07.2018)
      mw_xr[7] = 0;
      mw_xl[7] = 0;
      mw_yu[7] = 0;
      mw_yd[7] = 0;
     
  // from MW71 -> MWB2
  for(int i = 0; i < 5; ++i)
    {
      mw_xr[i+8] = vme_frs[9][12+i*4] & 0xfff;
      mw_xl[i+8] = vme_frs[9][13+i*4] & 0xfff;
      mw_yu[i+8] = vme_frs[9][14+i*4] & 0xfff;
      mw_yd[i+8] = vme_frs[9][15+i*4] & 0xfff;
    }



   //////////////////////////////////////
  // TPC part                         //
  //                                  //
  //////////////////////////////////////
  
  //ADC

  //TPC 1 at S2 (TPC 21) in vaccuum
   tpc_a[0][0]= vme_tpc[12][0] & 0xfff;
   tpc_a[0][1]= vme_tpc[12][1] & 0xfff;
   tpc_a[0][2]= vme_tpc[12][2] & 0xfff;
   tpc_a[0][3]= vme_tpc[12][3] & 0xfff;
   tpc_l[0][0]= vme_tpc[12][4] & 0xfff;
   tpc_r[0][0]= vme_tpc[12][5] & 0xfff;
   tpc_l[0][1]= vme_tpc[12][6] & 0xfff;
   tpc_r[0][1]= vme_tpc[12][7] & 0xfff;

  //TPC 2 at S2 (TPC 22) in vaccuum
   tpc_a[1][0]= vme_tpc[12][8] & 0xfff;
   tpc_a[1][1]= vme_tpc[12][9] & 0xfff;
   tpc_a[1][2]= vme_tpc[12][10] & 0xfff;
   tpc_a[1][3]= vme_tpc[12][11] & 0xfff;
   tpc_l[1][0]= vme_tpc[12][12] & 0xfff;
   tpc_r[1][0]= vme_tpc[12][13] & 0xfff;
   tpc_l[1][1]= vme_tpc[12][14] & 0xfff;
   tpc_r[1][1]= vme_tpc[12][15] & 0xfff;

  //TPC 3 at S2 (TPC 23) in air
   tpc_a[2][0]= vme_tpc[12][16] & 0xfff;
   tpc_a[2][1]= vme_tpc[12][17] & 0xfff;
   tpc_a[2][2]= vme_tpc[12][18] & 0xfff;
   tpc_a[2][3]= vme_tpc[12][19] & 0xfff;
   tpc_l[2][0]= vme_tpc[12][20] & 0xfff;
   tpc_r[2][0]= vme_tpc[12][21] & 0xfff;
   tpc_l[2][1]= vme_tpc[12][22] & 0xfff;
   tpc_r[2][1]= vme_tpc[12][23] & 0xfff;
  
  //TPC 4 at S2 (TPC 24) in air
   tpc_a[3][0]= vme_tpc[12][24] & 0xfff;
   tpc_a[3][1]= vme_tpc[12][25] & 0xfff;
   tpc_a[3][2]= vme_tpc[12][26] & 0xfff;
   tpc_a[3][3]= vme_tpc[12][27] & 0xfff;
   tpc_l[3][0]= vme_tpc[12][28] & 0xfff;
   tpc_r[3][0]= vme_tpc[12][29] & 0xfff;
   tpc_l[3][1]= vme_tpc[12][30] & 0xfff;
   tpc_r[3][1]= vme_tpc[12][31] & 0xfff;
  
  //TPC 5  at S4 (TPC 41) in air
   tpc_a[4][0]= vme_tpc[1][0] & 0xfff;
   tpc_a[4][1]= vme_tpc[1][1] & 0xfff;
   tpc_a[4][2]= vme_tpc[1][2] & 0xfff;
   tpc_a[4][3]= vme_tpc[1][3] & 0xfff;
   tpc_l[4][0]= vme_tpc[1][4] & 0xfff;
   tpc_r[4][0]= vme_tpc[1][5] & 0xfff;
   tpc_l[4][1]= vme_tpc[1][6] & 0xfff;
   tpc_r[4][1]= vme_tpc[1][7] & 0xfff;

  //TPC 6 at S4 (TPC 42) in air
   tpc_a[5][0]= vme_tpc[1][8] & 0xfff;
   tpc_a[5][1]= vme_tpc[1][9] & 0xfff;
   tpc_a[5][2]= vme_tpc[1][10] & 0xfff;
   tpc_a[5][3]= vme_tpc[1][11] & 0xfff;
   tpc_l[5][0]= vme_tpc[1][12] & 0xfff;
   tpc_r[5][0]= vme_tpc[1][13] & 0xfff;
   tpc_l[5][1]= vme_tpc[1][14] & 0xfff;
   tpc_r[5][1]= vme_tpc[1][15] & 0xfff;

  //TPC at S3 (TPC 31) (checkSB)
   tpc_a[6][0]= vme_tpc[1][16] & 0xfff;
   tpc_a[6][1]= vme_tpc[1][17] & 0xfff;
   tpc_a[6][2]= vme_tpc[1][18] & 0xfff;
   tpc_a[6][3]= vme_tpc[1][19] & 0xfff;
   tpc_l[6][0]= vme_tpc[1][20] & 0xfff;
   tpc_r[6][0]= vme_tpc[1][21] & 0xfff;
   tpc_l[6][1]= vme_tpc[1][22] & 0xfff;
   tpc_r[6][1]= vme_tpc[1][23] & 0xfff;
  
 
  //TDC

  //TPC 1 at S2 (TPC 21) in vaccuum
   tpc_dt[0][0]= vme_tpc[8][0] & 0xfff;
   tpc_dt[0][1]= vme_tpc[8][1] & 0xfff;
   tpc_dt[0][2]= vme_tpc[8][2] & 0xfff;
   tpc_dt[0][3]= vme_tpc[8][3] & 0xfff;
   tpc_lt[0][0]= vme_tpc[8][4] & 0xfff;
   tpc_rt[0][0]= vme_tpc[8][5] & 0xfff;
   tpc_lt[0][1]= vme_tpc[8][6] & 0xfff;
   tpc_rt[0][1]= vme_tpc[8][7] & 0xfff;
   
  //TPC 2 at S2 (TPC 22) in vaccuum
   tpc_dt[1][0]= vme_tpc[8][8] & 0xfff;
   tpc_dt[1][1]= vme_tpc[8][9] & 0xfff;
   tpc_dt[1][2]= vme_tpc[8][10] & 0xfff;
   tpc_dt[1][3]= vme_tpc[8][11] & 0xfff;
   tpc_lt[1][0]= vme_tpc[8][12] & 0xfff;
   tpc_rt[1][0]= vme_tpc[8][13] & 0xfff;
   tpc_lt[1][1]= vme_tpc[8][14] & 0xfff;
   tpc_rt[1][1]= vme_tpc[8][15] & 0xfff;

  //TPC 3 at S2 (TPC 23) in air
   tpc_dt[2][0]= vme_tpc[8][16] & 0xfff;
   tpc_dt[2][1]= vme_tpc[8][17] & 0xfff;
   tpc_dt[2][2]= vme_tpc[8][18] & 0xfff;
   tpc_dt[2][3]= vme_tpc[8][19] & 0xfff;
   tpc_lt[2][0]= vme_tpc[8][20] & 0xfff;
   tpc_rt[2][0]= vme_tpc[8][21] & 0xfff;
   tpc_lt[2][1]= vme_tpc[8][22] & 0xfff;
   tpc_rt[2][1]= vme_tpc[8][23] & 0xfff;

  //TPC 4 at S2 (TPC 24) in air
   tpc_dt[3][0]= vme_tpc[8][24] & 0xfff;
   tpc_dt[3][1]= vme_tpc[8][25] & 0xfff;
   tpc_dt[3][2]= vme_tpc[8][26] & 0xfff;
   tpc_dt[3][3]= vme_tpc[8][27] & 0xfff;
   tpc_lt[3][0]= vme_tpc[8][28] & 0xfff;
   tpc_rt[3][0]= vme_tpc[8][29] & 0xfff;
   tpc_lt[3][1]= vme_tpc[8][30] & 0xfff;
   tpc_rt[3][1]= vme_tpc[8][31] & 0xfff;

  //TPC 5 at S4 (TPC 41) in air  
   tpc_dt[4][0]= vme_tpc[0][0] & 0xfff;
   tpc_dt[4][1]= vme_tpc[0][1] & 0xfff;
   tpc_dt[4][2]= vme_tpc[0][2] & 0xfff;
   tpc_dt[4][3]= vme_tpc[0][3] & 0xfff;
   tpc_lt[4][0]= vme_tpc[0][4] & 0xfff;
   tpc_rt[4][0]= vme_tpc[0][5] & 0xfff;
   tpc_lt[4][1]= vme_tpc[0][6] & 0xfff;
   tpc_rt[4][1]= vme_tpc[0][7] & 0xfff;

  //TPC 6 at S4 (TPC 42) in air
   tpc_dt[5][0]= vme_tpc[0][8] & 0xfff;
   tpc_dt[5][1]= vme_tpc[0][9] & 0xfff;
   tpc_dt[5][2]= vme_tpc[0][10] & 0xfff;
   tpc_dt[5][3]= vme_tpc[0][11] & 0xfff;
   tpc_lt[5][0]= vme_tpc[0][12] & 0xfff;
   tpc_rt[5][0]= vme_tpc[0][13] & 0xfff;
   tpc_lt[5][1]= vme_tpc[0][14] & 0xfff;
   tpc_rt[5][1]= vme_tpc[0][15] & 0xfff;
//cout<<"tpc_rt " << tpc_rt[5][1] << endl;
  //TPC at S3 (TPC 31) (checkSB)
   tpc_dt[6][0]= vme_tpc[0][16] & 0xfff;
   tpc_dt[6][1]= vme_tpc[0][17] & 0xfff;
   tpc_dt[6][2]= vme_tpc[0][18] & 0xfff;
   tpc_dt[6][3]= vme_tpc[0][19] & 0xfff;
   tpc_lt[6][0]= vme_tpc[0][20] & 0xfff;
   tpc_rt[6][0]= vme_tpc[0][21] & 0xfff;
   tpc_lt[6][1]= vme_tpc[0][22] & 0xfff;
   tpc_rt[6][1]= vme_tpc[0][23] & 0xfff;
   //for(int i=0; i<2; i++) cout<<"SORT tpc_rt[6][1] " <<tpc_rt[6][i] << " tpc_lt[6][i]  " <<tpc_lt[6][i]<<" i " << i << endl;
   tpc_timeref[0] =  vme_tpc[0][24] & 0xfff;
   tpc_timeref[1] =  vme_tpc[0][25] & 0xfff;
   tpc_timeref[2] =  vme_tpc[0][26] & 0xfff;
   tpc_timeref[3] =  vme_tpc[0][27] & 0xfff;

//  printf("timeref[0] (sc41) = %d \n", tpc_timeref[0] );
//    for(int i=0; i<7; i++){
//        for(int j=0; j<2; j++){
//  cout<<"tpc_lt[i][j] " << tpc_lt[i][j] << " i " << i << " j " << j << endl;
//        }
//    }
    /// SCI dE: 
  
      de_41l = vme_main[11][0] & 0xfff;
      de_41r = vme_main[11][1] & 0xfff;
      de_21l = vme_main[11][2] & 0xfff;
      de_21r = vme_main[11][3] & 0xfff;
      de_42l = vme_main[11][4] & 0xfff; 
      de_42r = vme_main[11][5] & 0xfff;
  
    
      de_31l = vme_main[11][6] & 0xfff;
      de_31r = vme_main[11][7] & 0xfff;
      de_43l = vme_main[11][11] & 0xfff;
      de_43r = vme_main[11][12] & 0xfff;  
      de_81l = vme_main[11][13] & 0xfff;
      de_81r = vme_main[11][14] & 0xfff;
   
 
    
//     de_v1l = 0;
//     de_v1r = 0;  
//     de_v2l = 0;
//     de_v2r = 0;
//     de_v3  = 0;

    
    ///  SCI times: 
      dt_21l_21r = vme_frs[12][0] & 0xfff;
      dt_41l_41r = vme_frs[12][1] & 0xfff;
      dt_42l_42r = vme_frs[12][2] & 0xfff;
      dt_43l_43r = vme_frs[12][3] & 0xfff; //
      dt_81l_81r = vme_frs[12][4] & 0xfff;
      dt_21l_41l = vme_frs[12][5] & 0xfff;
      dt_21r_41r = vme_frs[12][6] & 0xfff;
      dt_42l_21l = vme_frs[12][8] & 0xfff; //swapped recently 28/02/2020 
      dt_42r_21r = vme_frs[12][7] & 0xfff; //swapped recently 28/02/2020
      dt_21l_81l = vme_frs[12][9] & 0xfff;
      dt_21r_81r = vme_frs[12][10] & 0xfff;
      
      dt_22l_22r = vme_frs[12][11] & 0xfff;
      dt_22l_41l = vme_frs[12][12] & 0xfff;
      dt_22r_41r = vme_frs[12][13] & 0xfff; ///ADDED 140920 AKM 
      
      
 
   // Time raw in ps from VFTX module
      if(psubevt->GetSubtype() ==1){
      
     float r=0.;
   for(int i=0; i<32; i++){
       if((vftx_cc[SCI_MOD][i][0]!=0) && (vftx_ft[SCI_MOD][i][0]!=0)){ 
           r = (double)rand.Rndm() - 0.5 ;
   TRaw_vftx[i] = VFTX_GetTraw_ps(SCI_MOD,i,vftx_cc[SCI_MOD][i][0],vftx_ft[SCI_MOD][i][0],r);
       }
   }
   TRaw_vftx_21l = TRaw_vftx[SCI21L_CH];
   TRaw_vftx_21r = TRaw_vftx[SCI21R_CH];
   TRaw_vftx_22l = TRaw_vftx[SCI22L_CH];
   TRaw_vftx_22r = TRaw_vftx[SCI22R_CH];
   TRaw_vftx_41l = TRaw_vftx[SCI41L_CH];
   TRaw_vftx_41r = TRaw_vftx[SCI41R_CH];
   TRaw_vftx_42l = TRaw_vftx[SCI42L_CH];
   TRaw_vftx_42r = TRaw_vftx[SCI42R_CH];
//    TRaw_vftx_sofia_l = TRaw_vftx[SCISOFIAL_CH];
//    TRaw_vftx_sofia_r = TRaw_vftx[SCISOFIAR_CH];
      }
  
  if(psubevt->GetSubtype() == 3600){
     for(int i=0; i<21; i++){ 
                for (int j=0; j<32;j++){ 
                   vme_main[i][j]=0; 
                  
                   
                }
           }
            for(int i=0; i<32; i++){ 
                for (int j=0; j<10;j++){ 
                  
                   leading_v1290_main[i][j]=0;
                   
                }
           }
}

  ///AKM 12.07.20 since we only take the first hit (for now...)
  /* ### SCI Multihit TDC time:  */
//     for(int i = 0; i < 10; ++i){
// //     cout<<"leading_v1290_main[0][i] " << leading_v1290_main[0][i] << " tdc_sc41l[i] " << tdc_sc41l[i] << " i " << i << " psubevt->GetProcid() " <<psubevt->GetProcid() << endl; 
//     tdc_sc41l[i] = leading_v1290_main[0][i];
//     tdc_sc41r[i] = leading_v1290_main[1][i];
//     
//     tdc_sc21l[i] = leading_v1290_main[2][i];
//     tdc_sc21r[i] = leading_v1290_main[3][i];
//    cout<<" tdc_sc21l[i] " << tdc_sc21l[i] << "  tdc_sc21r[i] " << tdc_sc21r[i]<< " i " << i << " leading_v1290_main[3][i] " <<leading_v1290_main[3][i]<< endl;
//     tdc_sc42l[i] = leading_v1290_main[4][i];
//     tdc_sc42r[i] = leading_v1290_main[5][i];
//     tdc_sc43l[i] = leading_v1290_main[6][i];
//     tdc_sc43r[i] = leading_v1290_main[7][i];
//     tdc_sc81l[i] = leading_v1290_main[8][i];
//     tdc_sc81r[i] = leading_v1290_main[9][i];
//     tdc_sc31l[i] = leading_v1290_main[10][i];
//     tdc_sc31r[i] = leading_v1290_main[11][i];
//     tdc_sc11[i]  = leading_v1290_main[12][i];
//     
//     }

    tdc_sc41l[0] = leading_v1290_main[0][0];
    tdc_sc41r[0] = leading_v1290_main[1][0];
    tdc_sc21l[0] = leading_v1290_main[2][0];
    tdc_sc21r[0] = leading_v1290_main[3][0];
  
    tdc_sc42l[0] = leading_v1290_main[4][0];
    tdc_sc42r[0] = leading_v1290_main[5][0];
    tdc_sc43l[0] = leading_v1290_main[6][0];
    tdc_sc43r[0] = leading_v1290_main[7][0];
    tdc_sc81l[0] = leading_v1290_main[8][0];
    tdc_sc81r[0] = leading_v1290_main[9][0];
    tdc_sc31l[0] = leading_v1290_main[10][0];
    tdc_sc31r[0] = leading_v1290_main[11][0];
    tdc_sc11[0]  = leading_v1290_main[12][0];
  

 //---MUSIC configuration. 2x TUM-MUSIC from FRS crate and 1 TRavel-MUsic from TRMU crate (2020/Jan/23, YT)
       for(int i=0;i<8;i++)
        {
          music_e1[i] = (vme_frs[3][i]) & 0xfff;   // 
          music_e2[i] = (vme_frs[3][8+i]) & 0xfff; // Travel-MUSIC (from special VME crate)
          music_e3[i] = (vme_trmu_adc[i]); // 
          
          music_t1[i] = leading_v1290_main[16+i][0] & 0xfff; //TUM-MUSIC
          music_t2[i] = leading_v1290_main[24+i][0] & 0xfff; //TUM-MUSIC
          music_t3[i] = (vme_trmu_tdc[i]);                   //Travel-MUSIC (from special VME crate)
        }
        
      /* ### MUSIC temp & pressure:  */
      // music_pres[0] = 0; 
      // music_temp[0] = 0;
      // music_pres[1] = 0;
      // music_temp[1] = 0;
      // music_pres[2] = 0;
      // music_temp[2] = 0;
    
      // //Channeltron detectors (vme must be adjusted)
      // ct_time = vme3s_MT[2][0] & 0xfff;
      // ct_signal = vme2scaler[20] & 0xfff; 
      // ct_trigger_DU = vme2scaler[5]& 0xfff;
      // ct_trigger_SY = vme2scaler[6]& 0xfff;
    
      // //Electron current measurement (vme must be adjused)
      // ec_signal = vme0[10][1]& 0xfff;
      
      // // mrtof
      // mrtof_start = vme3s_MT[0][0] ;//& 0x7ffff;
      // mrtof_stopDelay = vme3s_MT[1][0];// & 0x7ffff;
      // mrtof_stop = vme3s_MT[2][0] ;//& 0x7ffff;
    
    //Electron current measurement (vme must be adjused)
    //ec_signal = vme_frs[10][1]& 0xfff;
    
    // mrtof
   // mrtof_start = leading_v1290_main[0][0] ;//& 0x7ffff;
   // mrtof_stopDelay = leading_v1290_main[1][0];// & 0x7ffff;
   // mrtof_stop = leading_v1290_main[2][0] ;//& 0x7ffff;
    
    
}


void FRS_Detector_System::FRS_Calib(){
    ///Scalars
   if (psubevt->GetProcid() == 30 && psubevt->GetType() == 12)
    {
      if(13 == trigger || 12 == trigger) return; // skip spill trigger
    
      if(1==scaler_check_first_event){
        for(int ii=0; ii<64; ii++){
           scaler_previous[ii] = sc_long[ii];
           scaler_initial[ii] = sc_long[ii];
        }
        scaler_check_first_event = 0;
      }
      //  cout<<"sc_long[scaler_ch_1kHz] " << sc_long[scaler_ch_1kHz] << endl;
       time_in_ms  = sc_long[scaler_ch_1kHz]  - scaler_initial[scaler_ch_1kHz];
       spill_count = sc_long[scaler_ch_spillstart] - scaler_initial[scaler_ch_spillstart];
  
       ibin_for_s      = ((time_in_ms / 1000) % 1000) + 1;
       ibin_for_100ms  = ((time_in_ms / 100) % 4000) + 1;
       ibin_for_spill  = (spill_count % 1000) +1; 
       
       for(int kk=0; kk<64; kk++){
       increase_scaler_temp[kk] =  sc_long[kk]  - scaler_previous[kk];
     // if(increase_scaler_temp>0) cout<<"increase_scaler_temp " <<increase_scaler_temp <<" sc_long[kk] " <<sc_long[kk]<< " scaler_previous[kk] " << scaler_previous[kk] <<  " kk " << kk << endl;
       }
      extraction_time_ms += sc_long[scaler_ch_1kHz] - scaler_previous[scaler_ch_1kHz];
      if(0 != sc_long[scaler_ch_spillstart] - scaler_previous[scaler_ch_spillstart]){
        extraction_time_ms = 0; 
      }

      
       ibin_clean_for_s      = (((time_in_ms / 1000) +20) % 1000) + 1; // 
       ibin_clean_for_100ms  = (((time_in_ms / 100) +200 ) % 4000) + 1; //
       ibin_clean_for_spill  = ((spill_count + 990 )% 20) +1; // 
      
      // put current data into _previous for the next event
      for(int ii=0; ii<64; ii++){
        scaler_previous[ii] = sc_long[ii];
      }
    }
  /**  for SeeTram calibration purposes:   **/


  //10HzClock = hMON_scaler[4]->Integral();

  /*   Naming conventions:  index     detector   focus #                     */
  /*                         1         MW11        1                         */
  /*                         2         MW21        2                         */
  /*                         3         MW22                                  */
  /*                         4         MW31        3                         */
  /*                         5         MW41                                  */
  /*                         6         MW42                                  */
  /*                         7         MW51                                  */
  /*                         8         MW61                                  */
   /*                        9         MW71                                  */
  /*                         10        MW81                                  */
  /*                         11        MW82                                  */
  /*                         12        MWB1                                  */
  /*                         13        MWB2                                  */
  /////////////////////////////////////////////////////////////////////////////
  int max_index = 13; //upto MWB2 (09.07.2018) /*  up to S3 */

  for(int i=0;i<max_index;i++)
    {   
      
      /********************************************/
      /* Calculate the sum spectra and conditions */
      /********************************************/
      
      /* better test first existence of xl, xr, an before filling of xsum */
      
      if(mw_an[i] && mw_xl[i] && mw_xr[i])
    {
      mw_xsum[i] = 1000+(mw_xl[i] - mw_an[i]) + (mw_xr[i] - mw_an[i]);
      
    }

    b_mw_xsum[i] = Check_WinCond_Multi(mw_xsum[i], lim_xsum, i);

      /* better test first existence of yu, yd, an before filling of ysum */
      
      if(mw_an[i] && mw_yu[i] && mw_yd[i])
    {
       
      //if(mw_yu[i]&&mw_yd[i]) {
      mw_ysum[i] = 1000+(mw_yu[i] - mw_an[i]) + (mw_yd[i] - mw_an[i]);    
      
      //mw_ysum[i] = (mw_yu[i]) + (mw_yd[i]); //when an doesn't work            
    }
      
        b_mw_ysum[i] = Check_WinCond_Multi(mw_ysum[i], lim_ysum, i);

    //if(mw_ysum[i] >= lim_ysum[i][0] && mw_ysum[i] <= lim_ysum[i][1]) b_mw_ysum[i] =  true; //cMW_YSUM[i]->Test(mw_ysum[i]);
    //else b_mw_ysum[i] =  false;
    
      /*******************************************************************/
      /* If the signals in x and y are valid, calculate position spectra */
      /*******************************************************************/
      
      if (b_mw_xsum[i])
    {
      //      Int_t r_x = mw_xl[i] - mw_xr[i];
      Float_t r_x = mw_xl[i] *  mw->gain_tdc[1][i] - mw_xr[i] *  mw->gain_tdc[2][i]; //14.09.05 CN+AM
      mw_x[i] = mw->x_factor[i] * r_x + mw->x_offset[i];
          
    }
    
      if (b_mw_ysum[i])
    {
      //      Int_t r_y = mw_yd[i] - mw_yu[i];
      Float_t r_y = mw_yd[i] *  mw->gain_tdc[4][i] - mw_yu[i] *  mw->gain_tdc[3][i]; //14.09.05 CN+AM
      mw_y[i] = mw->y_factor[i] * r_y + mw->y_offset[i];
    }

    
      if(b_mw_xsum[i] && b_mw_ysum[i])
    {
        
    }
      
    } // for(int i=0;i<max_index;i++)


  /*********************************************/
  /* S2 Angle and transformed position spectra */
  /*********************************************/
  
  Float_t  dist_MW21_MW22  = frs->dist_MW22  - frs->dist_MW21;
  Float_t  dist_MW22_focS2 = frs->dist_focS2 - frs->dist_MW22;
  Float_t  dist_MW22_SC21  = frs->dist_SC21  - frs->dist_MW22;
  Float_t  dist_MW22_SC22  = frs->dist_SC22  - frs->dist_MW22;
  

  if (b_mw_xsum[1] && b_mw_xsum[2])
    {
      /* MW21 and MW22 X okay */
      /*  X angle at S2  [mrad]:  */
      angle_x_s2 = (mw_x[2]-mw_x[1])/dist_MW21_MW22*1000.;

      /*  X at nominal S2 focus:  */
      focx_s2 = mw_x[2] + dist_MW22_focS2 * angle_x_s2/1000.;
      
      /*  X at SC21 position:    */
      mw_sc21_x = mw_x[2] + dist_MW22_SC21 * angle_x_s2/1000.;

      /*  X at SC22 position:    */
      mw_sc22_x = mw_x[2] + dist_MW22_SC22 * angle_x_s2/1000.;
      

      /* 'real' z-position of S2 X focus (cm) */
      Float_t rh = (angle_x_s2 - angle_x_s2m);
      if(fabs(rh)>1e-4)
    {
      z_x_s2 = ((focx_s2m - focx_s2)/rh)*100. + frs->dist_focS2/10.;  
    }
    
      /* keep values for next event */
      focx_s2m = focx_s2; 
      angle_x_s2m = angle_x_s2;
      
    
      
    }
  
  if (b_mw_ysum[1] && b_mw_ysum[2])
    {
      /* MW21 and MW22 Y okay */
      /*  Y angle at S2 [mrad]:   */
      angle_y_s2 = (mw_y[2] - mw_y[1])/dist_MW21_MW22*1000.;
      
      /*  Y at nominal S2 focus:  */
      focy_s2 = mw_y[2] + dist_MW22_focS2 * angle_y_s2/1000.;

      /* 'real' z-position of S2 Y focus (cm) */
      Float_t rh = (angle_y_s2 - angle_y_s2m);
      if(fabs(rh)>1.e-4)
    {
      z_y_s2 = ((focy_s2m - focy_s2)/rh)*100. + frs->dist_focS2/10.;  
    }
      
      /* keep values for next event */
      focy_s2m = focy_s2; 
      angle_y_s2m = angle_y_s2; 
    }
  
  /*if (b_mw_ysum[1] && b_mw_ysum[2] && b_mw_xsum[1] && b_mw_xsum[2])
    hMW_FocS2->Fill(focx_s2,focy_s2);*/ 


//////////////////////////////////////////////////////////////////
///  ***********TPC Analysis*************************** //////////

 //================================
  // Nomenclature for TPCs
  //================================
  // TPCs at S2
  // TPC 1 = TPC 21 (in vaccum) [i=0]
  // TPC 2 = TPC 22 (in vaccum) [i=1]
  // TPC 3 = TPC 23 (in air)    [i=2]
  // TPC 4 = TPC 24 (in air)    [i=3]

  // TPCs at S4
  // TPC 5 = TPC 41 (in air)    [i=4]
  // TPC 6 = TPC 42 (in air)    [i=5]

  //TPC at S3
  // TPC 7 = TPC 31             [i=6]
  

  for(int i=0;i<7;i++)
    {  
//cout<<"mw->x_factor[0] " << mw->x_factor[0] << endl;
      int count =0;
      Int_t r_y = 0;
      
      int countx =0;
      Int_t r_x0 = 0;
      Int_t r_x1 = 0;
      
      tpc_y[i] = 0;
      tpc_x[i] = 9999;

      for(int j=0;j<4;j++)
    {

      //// calculate control sums
      if(j < 2)
        tpc_csum[i][j] = (tpc_lt[i][0] + tpc_rt[i][0]- 2*tpc_dt[i][j]);
      else
        tpc_csum[i][j] = (tpc_lt[i][1] + tpc_rt[i][1]- 2*tpc_dt[i][j]);


           // if((de_42l>230&&de_42l<450)||(de_42r>540&&de_42r<750)){

      //b_tpc_csum[i][j] = cTPC_CSUM[i][j]->Test(tpc_csum[i][j]);
      
      b_tpc_csum[i][j] = Check_WinCond_Multi_Multi(tpc_csum[i][j], lim_csum, j, i);
      
      
      //if(tpc_csum[i][j] >= lim_csum[j][i][0] && tpc_csum[i][j] <= lim_csum[j][i][1]) b_tpc_csum[i][j] = true;
      //else b_tpc_csum[i][j] = false;
           // }
      
      if(tpc_lt[i][0]==0 && tpc_rt[i][0]==0 && j<2)
        b_tpc_csum[i][j]=0;
      
      if(tpc_lt[i][1]==0 && tpc_rt[i][1]==0 && j>1)
        b_tpc_csum[i][j]=0;
      
      if (b_tpc_csum[i][j])
        {
          r_y += tpc_dt[i][j];
          tpc_y[i] += tpc->y_factor[i][j]*tpc_dt[i][j] + tpc->y_offset[i][j];
          count++;
        }
       // cout<<"tpc_lt[6][j] " <<tpc_lt[6][j] << " tpc_rt[6][j] " <<tpc_rt[6][j] <<" i " << i << " j " << j << endl;
    }

      if (b_tpc_csum[i][0] || b_tpc_csum[i][1])
    {
      r_x0 =  tpc_lt[i][0]-tpc_rt[i][0];
      //          r_x0 =  tpc_lt[i][0];
     
      tpc_x[i]=tpc->x_factor[i][0]*r_x0 + tpc->x_offset[i][0];
     // cout<<"1111i " << i << " r_x0 " << r_x0 <<" tpc_x[i] " <<tpc_x[i] << endl;
      countx++;
    }
      
      if (b_tpc_csum[i][2] || b_tpc_csum[i][3])
    {
      r_x1 =  tpc_lt[i][1]-tpc_rt[i][1];
      if (countx == 0){
        tpc_x[i]=tpc->x_factor[i][1]*r_x1 + tpc->x_offset[i][1];
     //  cout<<"2222i " << i << " r_x0 " << r_x1 <<" tpc_x[i] " <<tpc_x[i] << endl;
      }
      else{
        tpc_x[i]+=tpc->x_factor[i][1]*r_x1 + tpc->x_offset[i][1];
     // cout<<"3333i " << i << " r_x1 " << r_x1 <<" tpc_x[i] " <<tpc_x[i] << endl;
      countx++;
      }
    }

      b_tpc_xy[i] = kFALSE;  // HWE

      if (countx > 0)
    {
      tpc_y[i] = tpc_y[i]/count;
      tpc_x[i] = tpc_x[i]/countx;
    // cout<<"4444i " << i << " r_x1 " << r_x1 <<" tpc_x[i] " <<tpc_x[i] << endl;
      //Int_t r_x=tpc_lt[i]-tpc_rt[i];
      //tpc_x[i]=tpc->x_factor[i]*r_x + tpc->x_offset[i];
    
      //      if(r_x0<-40)

      b_tpc_xy[i] = kTRUE;
    }
      
      if(countx>1)
    {
       x0=tpc->x_factor[i][0]*r_x0 + tpc->x_offset[i][0];
       x1=tpc->x_factor[i][1]*r_x1 + tpc->x_offset[i][1];

    }    
    }

  //*********************************************
  // *****Relative distances for tracking *******
  //*********************************************
  //Float_t dist_TPC1_TPC2 = 1.;         
  //distances at S2 for exp s388 measured 'by eye' (i.e. not very precise)! AE, 8.8.2012
  //Float_t dist_TPC4_target1 = 455.;  //8.8.12, s388: taget1 is first Si detector at S2 
  //using: TPC4->chamber= 150 mm, chamber->Si(1)= 305 mm
  //Float_t dist_TPC3_TPC4 = 610.; // exp S417+S411 Oct2014
  //Float_t dist_TPC3_focS2 = 2860 - 2013 ; //exp S411 Oct2014
  //Float_t dist_TPC3_focS2 = 0. ; //exp S417 Oct2014 focus on TPC21
  //Float_t dist_TPC5_TPC6 = 1595. - 385. ; // Oct.2014, exp s411 
  //Float_t dist_TPC6_focS4  = 2349. - 1595. ; // Oct.2014, exp s411
  //THIS SHOULD BE NOT HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Yhis was for Uranium--- see logbook page 6
  //Float_t dist_TPC6_target2  = 1015.; // 8.08.12, exp s388 target2= OTPC entrance (check distance)
  //Float_t dist_TPC4_target1 = 455.;  // 8.8.12, s388: taget1 is first Si detector at S2 
                                     // using: TPC4->chamber= 150 mm, chamber->Si(1)= 305 mm

  //================================
  // Distances of TPCs at S2 and S4
  //================================

  Float_t dist_TPC21_TPC22 = frs->dist_TPC22 - frs->dist_TPC21; 
  Float_t dist_TPC23_TPC24 = frs->dist_TPC24 - frs->dist_TPC23; 
  Float_t dist_TPC22_TPC24 = frs->dist_TPC24 - frs->dist_TPC22; 
  Float_t dist_TPC21_focS2 = frs->dist_TPC21 - frs->dist_focS2; 
  Float_t dist_TPC22_focS2 = frs->dist_TPC22 - frs->dist_focS2; 
  Float_t dist_TPC23_focS2 = frs->dist_TPC23 - frs->dist_focS2; 
  Float_t dist_TPC41_TPC42 = frs->dist_TPC42 - frs->dist_TPC41; 
  Float_t dist_TPC42_focS4 = frs->dist_focS4 - frs->dist_TPC42; 


  //Float_t dist_TPC3_TPC4 = 1782.0- 604.; //old value for old tracking 
  //Float_t dist_TPC3_focS2 = 604.-2165. ; //old value for old tracking  
/*
  //==================================
  //Distance of TPCs at S4
  //==================================
  
  Float_t dist_TPC41_TPC42 = 1525.0 - 415.0;  //rough check
  //Float_t dist_TPC41_TPC42 = 1660.0 - 360.0;  //check

  Float_t dist_TPC42_focS4 = 3300.0 - 1525.0; //check
  //Float_t dist_TPC42_focS4 = 3300.0 - 1660.0; //check*/

  
  //Float_t dist_TPC5_TPC6 = 2063.5 - 382.5 ; 
  //Float_t dist_TPC6_focS4  = 2200. - 2063.5 ;
  //Float_t dist_TPC6_focS4  = 4300. - 1660. ;//before 13Jun2016 17:10
  //Float_t dist_TPC6_focS4  = 2200. - 1660. ;//for S4 focus 2200 , changed at 13Jun2016 17:10 
  //Float_t dist_TPC6_focS4  = 2200. -1135.9 +300.0 - 1660. +100.; //short focus 11jul2016
  //Float_t dist_TPC6_focS4  = 2202 - 1660.; // tensor28 mode
  //Float_t dist_TPC6_focS4  = 3300 - 1660.; // long focus

  //=====================
  // Old tracking for S2 
  //=====================
  /*
  //Position S2 tracked with TPCs 3 and 4 for focus
  if (b_tpc_xy[2]&&b_tpc_xy[3])
    {
      tpc_angle_x_s2_foc = (tpc_x[3] - tpc_x[2])/dist_TPC3_TPC4*1000.;
      tpc_angle_y_s2_foc = (tpc_y[3] - tpc_y[2])/dist_TPC3_TPC4*1000.;
      tpc_x_s2_foc = -tpc_angle_x_s2_foc * dist_TPC3_focS2/1000. + tpc_x[2];
      tpc_y_s2_foc = -tpc_angle_y_s2_foc * dist_TPC3_focS2/1000. + tpc_y[2];

 
     // Position S2 tracked with TPCs 3 and 4 to first Si tracker in exp S388
     // tpc_x_s2_target1 = tpc_angle_x_s2_foc * dist_TPC4_target1/1000. + tpc_x[3];
     // tpc_y_s2_target1 = tpc_angle_y_s2_foc * dist_TPC4_target1/1000. + tpc_y[3];
      
    }
  else  tpc_x_s2_foc= -999 ;
  */

  //=================================
  // Tracking with TPC 21 and TPC 22 
  //=================================
  
  if (b_tpc_xy[0]&&b_tpc_xy[1])
    {
      tpc_angle_x_s2_foc_21_22 = (tpc_x[1] - tpc_x[0])/dist_TPC21_TPC22*1000.;
      tpc_angle_y_s2_foc_21_22 = (tpc_y[1] - tpc_y[0])/dist_TPC21_TPC22*1000.;
      tpc_x_s2_foc_21_22 = -tpc_angle_x_s2_foc_21_22 * dist_TPC21_focS2/1000. + tpc_x[0]; //check
      tpc_y_s2_foc_21_22 = -tpc_angle_y_s2_foc_21_22 * dist_TPC21_focS2/1000. + tpc_y[0]; //check
      
        ///This stuff is new, check why its needed 
      Float_t dist_SC21_focS2 = frs->dist_SC21 - frs->dist_focS2; 
      
      tpc21_22_sc21_x = (tpc_angle_x_s2_foc_21_22/1000.*dist_SC21_focS2)+tpc_x_s2_foc_21_22;
      tpc21_22_sc21_y = (tpc_angle_y_s2_foc_21_22/1000.*dist_SC21_focS2)+tpc_y_s2_foc_21_22; 
      
      Float_t dist_S2target_focS2 = frs->dist_S2target - frs->dist_focS2;   
      
      tpc21_22_s2target_x = (tpc_angle_x_s2_foc_21_22/1000.*dist_S2target_focS2)+tpc_x_s2_foc_21_22;
      tpc21_22_s2target_y = (tpc_angle_y_s2_foc_21_22/1000.*dist_S2target_focS2)+tpc_y_s2_foc_21_22;  
      
    }
  else  tpc_x_s2_foc_21_22= -999 ;


  //=================================
  // Tracking with TPC 23 and TPC 24 
  //=================================
  
  if (b_tpc_xy[2]&&b_tpc_xy[3])
    {
      tpc_angle_x_s2_foc_23_24 = (tpc_x[3] - tpc_x[2])/dist_TPC23_TPC24*1000.;
      tpc_angle_y_s2_foc_23_24 = (tpc_y[3] - tpc_y[2])/dist_TPC23_TPC24*1000.;
      tpc_x_s2_foc_23_24 = -tpc_angle_x_s2_foc_23_24 * dist_TPC23_focS2/1000. + tpc_x[2]; //check
      tpc_y_s2_foc_23_24 = -tpc_angle_y_s2_foc_23_24 * dist_TPC23_focS2/1000. + tpc_y[2]; //check
      
      ///This stuff is new, check why its needed 
       Float_t dist_SC21_focS2 = frs->dist_SC21 - frs->dist_focS2;
       tpc23_24_sc21_x = (tpc_angle_x_s2_foc_23_24/1000.*dist_SC21_focS2)+tpc_x_s2_foc_23_24;
       tpc23_24_sc21_y = (tpc_angle_y_s2_foc_23_24/1000.*dist_SC21_focS2)+tpc_y_s2_foc_23_24;
       Float_t dist_S2target_focS2 = frs->dist_S2target - frs->dist_focS2;   
       tpc23_24_s2target_x = (tpc_angle_x_s2_foc_23_24/1000.*dist_S2target_focS2)+tpc_x_s2_foc_23_24;
       tpc23_24_s2target_y = (tpc_angle_y_s2_foc_23_24/1000.*dist_S2target_focS2)+tpc_y_s2_foc_23_24; 

    }
  else  tpc_x_s2_foc_23_24= -999 ;


  //=================================
  // Tracking with TPC 22 and TPC 24 
  //=================================
  
  if (b_tpc_xy[1]&&b_tpc_xy[3])
    {
      tpc_angle_x_s2_foc_22_24 = (tpc_x[3] - tpc_x[1])/dist_TPC22_TPC24*1000.;
      tpc_angle_y_s2_foc_22_24 = (tpc_y[3] - tpc_y[1])/dist_TPC22_TPC24*1000.;
      tpc_x_s2_foc_22_24 = -tpc_angle_x_s2_foc_22_24 * dist_TPC22_focS2/1000. + tpc_x[1]; //check
      tpc_y_s2_foc_22_24 = -tpc_angle_y_s2_foc_22_24 * dist_TPC22_focS2/1000. + tpc_y[1]; //check
      
      ///This stuff is new, check why its needed 
        Float_t dist_SC21_focS2 = frs->dist_SC21 - frs->dist_focS2;   
      tpc22_24_sc21_x = (tpc_angle_x_s2_foc_22_24/1000.*dist_SC21_focS2)+tpc_x_s2_foc_22_24;
      tpc22_24_sc21_y = (tpc_angle_y_s2_foc_22_24/1000.*dist_SC21_focS2)+tpc_y_s2_foc_22_24;     
      Float_t dist_S2target_focS2 = frs->dist_S2target - frs->dist_focS2;   
      tpc22_24_s2target_x = (tpc_angle_x_s2_foc_22_24/1000.*dist_S2target_focS2)+tpc_x_s2_foc_22_24;
      tpc22_24_s2target_y = (tpc_angle_y_s2_foc_22_24/1000.*dist_S2target_focS2)+tpc_y_s2_foc_22_24; 
      
    }
  else  tpc_x_s2_foc_22_24= -999 ;
  
  
  //=====================================================
  // Tracking with TPC 41 and TPC 42 (TPC 5 and 6) at S4  
  //=====================================================
  if (b_tpc_xy[4]&&b_tpc_xy[5])
    {
      
      tpc_angle_x_s4 = (tpc_x[5] - tpc_x[4])/dist_TPC41_TPC42*1000.;
      tpc_angle_y_s4 = (tpc_y[5] - tpc_y[4])/dist_TPC41_TPC42*1000.;
      tpc_x_s4 = tpc_angle_x_s4 * dist_TPC42_focS4/1000. + tpc_x[5];
      tpc_y_s4 = tpc_angle_y_s4 * dist_TPC42_focS4/1000. + tpc_y[5];
    
      
      //This is new stuff, check
        Float_t dist_SC41_focS4 = frs->dist_SC41 - frs->dist_focS4;
       tpc_sc41_x = (tpc_angle_x_s4/1000.*dist_SC41_focS4)+tpc_x_s4;
       tpc_sc41_y = (tpc_angle_y_s4/1000.*dist_SC41_focS4)+tpc_y_s4;
       
       Float_t dist_SC42_focS4 = frs->dist_SC42 - frs->dist_focS4;
       tpc_sc42_x = (tpc_angle_x_s4/1000.*dist_SC42_focS4)+tpc_x_s4;
       tpc_sc42_y = (tpc_angle_y_s4/1000.*dist_SC42_focS4)+tpc_y_s4;
       
       Float_t dist_SC43_focS4 = frs->dist_SC43 - frs->dist_focS4;
       tpc_sc43_x = (tpc_angle_x_s4/1000.*dist_SC43_focS4)+tpc_x_s4;
       tpc_sc43_y = (tpc_angle_y_s4/1000.*dist_SC43_focS4)+tpc_y_s4;
       
       Float_t dist_MUSIC41_focS4 = frs->dist_MUSIC41 - frs->dist_focS4;
       tpc_music41_x = (tpc_angle_x_s4/1000.*dist_MUSIC41_focS4)+tpc_x_s4;
       tpc_music41_y = (tpc_angle_y_s4/1000.*dist_MUSIC41_focS4)+tpc_y_s4;
       
       Float_t dist_MUSIC42_focS4 = frs->dist_MUSIC42 - frs->dist_focS4;
       tpc_music42_x = (tpc_angle_x_s4/1000.*dist_MUSIC42_focS4)+tpc_x_s4;
       tpc_music42_y = (tpc_angle_y_s4/1000.*dist_MUSIC42_focS4)+tpc_y_s4;
       
       Float_t dist_MUSIC43_focS4 = frs->dist_MUSIC43 - frs->dist_focS4;
       tpc_music43_x = (tpc_angle_x_s4/1000.*dist_MUSIC43_focS4)+tpc_x_s4;
       tpc_music43_y = (tpc_angle_y_s4/1000.*dist_MUSIC43_focS4)+tpc_y_s4;
       
       Float_t dist_S4target_focS4 = frs->dist_S4target - frs->dist_focS4;
       tpc_s4target_x = (tpc_angle_x_s4/1000.* dist_S4target_focS4)+tpc_x_s4;
       tpc_s4target_y = (tpc_angle_y_s4/1000.* dist_S4target_focS4)+tpc_y_s4;
       
//       music1_x1 = tpc_angle_x_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa1)/1000. + tpc_x[5];
//       music1_x2 = tpc_angle_x_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa2)/1000. + tpc_x[5];
//       music1_x3 = tpc_angle_x_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa3)/1000. + tpc_x[5];
//       music1_x4 = tpc_angle_x_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa4)/1000. + tpc_x[5];
// 
//       music1_y1 = tpc_angle_y_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa1)/1000. + tpc_y[5];
//       music1_y2 = tpc_angle_y_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa2)/1000. + tpc_y[5];
//       music1_y3 = tpc_angle_y_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa3)/1000. + tpc_y[5];
//       music1_y4 = tpc_angle_y_s4 * (frs->dist_MUSIC3+frs->dist_MUSICa4)/1000. + tpc_y[5];

      
      //S4 entrance of OTPC tracked with TPCs 5 and 6 for s388    
      //tpc_x_s4_target2 = tpc_angle_x_s4 * dist_TPC6_target2/1000. + tpc_x[5];
      //tpc_y_s4_target2 = tpc_angle_y_s4 * dist_TPC6_target2/1000. + tpc_y[5];
    
    }
  
  
  // tracking to musics and scis (S323,410 exp.)
  //Float_t dist_TPC3_SC21  = 1375.5;
  //Float_t dist_TPC6_SC41  = 510.; //8.09.11 //check AEv

  // X position at SC21 position:
  // SC21 is before TPC3 and TPC4
  //tpc_mw_sc21_x = tpc_x[2] - dist_TPC3_SC21 * tpc_angle_x_s2_foc/1000.;
 

  // Y position at SC21
  //tpc_sc21_y = tpc_y[2] - dist_TPC3_SC21 * tpc_angle_y_s2_foc/1000.;

  
  // X position at SC41 
  // SC41 is after TPC6
 // tpc_sc41_x  = tpc_x[5] + dist_TPC6_SC41 * tpc_angle_x_s4 / 1000.;
  
  // Y position at SC41
  //tpc_sc41_y = tpc_y[5] + dist_TPC6_SC41 * tpc_angle_y_s4 /1000.;


//////////////////////////////////////////////////////////////////
///  ***********Si Analysis**************************** //////////


  // Si alpha
//   si_e1 = (si->si_factor1*si_adc1)+si->si_offset1;
//   si_e2 = (si->si_factor2*si_adc2)+si->si_offset2;
//   si_e3 = (si->si_factor3*si_adc3)+si->si_offset3;
//   si_e4 = (si->si_factor4*si_adc4)+si->si_offset4;
//   si_e5 = (si->si_factor5*si_adc5)+si->si_offset5;


  //========================
  //Active Stopper Analysis
  //========================
  
//   dssd_xmaxindex=-10;
//   dssd_xmaxenergy=1000;

  
  //for(int i=0;i<16;i++){
  //dssd_e[15-i] = (si->dssd_factor[15-i]*(dssd_adc[i]+gRandom->Rndm()-0.5))+si->dssd_offset[15-i];       }
  

//   for(int i=0;i<16;i++){
//    //dssd_e[i] = dssd_adc[i];
//    //  dssd_e[i] = (si->dssd_factor[i]*dssd_adc[i])+si->dssd_offset[i];
//   }


//   for (int i=0 ; i<16 ; i++)
//       {
//         if (dssd_xmaxenergy<dssd_e[i])
//            {
//       dssd_xmaxenergy=dssd_e[i];
//       dssd_xmaxindex= i;
//     }
//       }
  
    
//   dssd_ymaxindex=-10;
//   dssd_ymaxenergy=1000;
//   
  
  //for(int i=16;i<24;i++){
  //dssd_e[i+8] = (si->dssd_factor[i+8]*(dssd_adc[i]+gRandom->Rndm()-0.5))+si->dssd_offset[i+8];
  //}
    
  //for(int i=24;i<32;i++){
  //dssd_e[i-8] = (si->dssd_factor[i-8]*(dssd_adc[i]+gRandom->Rndm()-0.5))+si->dssd_offset[i-8];
  //}   
   

//    for(int i=16;i<32;i++){
//       // dssd_e[i] = dssd_adc[i];
//   dssd_e[i] = (si->dssd_factor[i]*dssd_adc[i])+si->dssd_offset[i];
//    }
// 
// 
//    for (int i=16 ; i<32 ; i++){
//      if(dssd_ymaxenergy<dssd_e[i])
//        {
//   dssd_ymaxenergy=dssd_e[i];
//   dssd_ymaxindex= i-16;
//        } 
// 
//    }
    
    
}



void FRS_Detector_System::FRS_Anal(){
    
    ///==================================================================================///
                    /// Start of MUSIC  analysis
///==================================================================================///
  music1_anodes_cnt = 0;  
  music2_anodes_cnt = 0;  
  music3_anodes_cnt = 0;


   // Munich MUSIC 
  
  for (int i=0;i<8;i++)
    {
      /* 8 anodes of TUM MUSIC */       
      /****** first MUSIC ***** threshold changed to 5, 9/8/2012 **/
      if ( music_e1[i] > 5)
    { 

      music_b_e1[i] = Check_WinCond_Multi(music_e1[i], cMusic1_E, i);//cMusic1_E[i]->Test(music_e1[i]);
      

      if (music_b_e1[i])
        music1_anodes_cnt++;
    }

        if (music_t1[i] > 0)
    { 
      music_b_t1[i] = Check_WinCond_Multi(music_t1[i], cMusic1_T, i); // cMusic1_T[i]->Test(music_t1[i]);
    }

       
      //hMUSIC1_dE1dE2->Fill(music_e1[0],music_e1[1]);

      //       /****** second MUSIC *****/
     
      if ( music_e2[i] > 5)
    { 
      music_b_e2[i] = Check_WinCond_Multi(music_e2[i], cMusic2_E, i);// cMusic2_E[i]->Test(music_e2[i]);
      if (music_b_e2[i]) music2_anodes_cnt++;
    }

      if (music_t2[i] > 0)
    { 
      music_b_t2[i] = Check_WinCond_Multi(music_t2[i], cMusic2_T, i);// cMusic2_T[i]->Test(music_t2[i]);
    }
    }

  for (int i=0;i<4;i++)
    {
      /* 4 anodes of MUSIC OLD */       
      /****** third MUSIC *****/
      if ( music_e3[i] > 0)
    { 

      music_b_e3[i] = Check_WinCond_Multi(music_e3[i], cMusic3_E, i);// cMusic3_E[i]->Test(music_e3[i]);
      if (music_b_e3[i])
        music3_anodes_cnt++;
    }
      if (music_t3[i] > 0)
    { 
      music_b_t3[i] = Check_WinCond_Multi(music_t3[i], cMusic3_T, i);// cMusic3_T[i]->Test(music_t3[i]);
    }
    }
  
  // calculate truncated dE from 8 anodes, Munich MUSIC 
  
  if (music1_anodes_cnt == 8)
    {
      /* "quick-n-dirty solution according to Klaus:   */
      //      Float_t r1 = (music_e1[0] - music->e1_off[0])*(music_e1[1] - music->e1_off[1]);
      //      Float_t r2 = (music_e1[2] - music->e1_off[2])*(music_e1[3] - music->e1_off[3]);
      //      Float_t r3 = (music_e1[4] - music->e1_off[4])*(music_e1[5] - music->e1_off[5]);
      //      Float_t r4 = (music_e1[6] - music->e1_off[6])*(music_e1[7] - music->e1_off[7]);

      Float_t r1 = ((music_e1[0])*music->e1_gain[0] + music->e1_off[0])*((music_e1[1])*music->e1_gain[1] + music->e1_off[1]);
      Float_t r2 = ((music_e1[2])*music->e1_gain[2] + music->e1_off[2])*((music_e1[3])*music->e1_gain[3] + music->e1_off[3]);
      Float_t r3 = ((music_e1[4])*music->e1_gain[4] + music->e1_off[4])*((music_e1[5])*music->e1_gain[5] + music->e1_off[5]);
      Float_t r4 = ((music_e1[6])*music->e1_gain[6] + music->e1_off[6])*((music_e1[7])*music->e1_gain[7] + music->e1_off[7]);

      if ( (r1 > 0) && (r2 > 0) && (r3 > 0) && (r4 > 0) )
    {
      b_de1 = kTRUE;
      de[0] = sqrt( sqrt( sqrt(r1) * sqrt(r2) ) * sqrt( sqrt(r3) * sqrt(r4) ) );
      de_cor[0] = de[0];
    }  
    }

      
  if (music2_anodes_cnt == 8)
    {
      
      Float_t r1 = ((music_e2[0])*music->e2_gain[0] + music->e2_off[0])*((music_e2[1])*music->e2_gain[1] + music->e2_off[1]);
      Float_t r2 = ((music_e2[2])*music->e2_gain[2] + music->e2_off[2])*((music_e2[3])*music->e2_gain[3] + music->e2_off[3]);
      Float_t r3 = ((music_e2[4])*music->e2_gain[4] + music->e2_off[4])*((music_e2[5])*music->e2_gain[5] + music->e2_off[5]);
      Float_t r4 = ((music_e2[6])*music->e2_gain[6] + music->e2_off[6])*((music_e2[7])*music->e2_gain[7] + music->e2_off[7]);

      
      if ( (r1 > 0) && (r2 > 0) && (r3 > 0) && (r4 > 0) )
    {
      b_de2 = kTRUE;
      de[1] = sqrt( sqrt( sqrt(r1) * sqrt(r2) ) * sqrt( sqrt(r3) * sqrt(r4) ) );
      de_cor[1] = de[1];
    }  
    }
       

  if (music3_anodes_cnt == 8)
    {         // OLD MUSIC
      
      Float_t r1 = ((music_e3[0])*music->e3_gain[0] + music->e3_off[0])*((music_e3[1])*music->e3_gain[1] + music->e3_off[1]);
      Float_t r2 = ((music_e3[2])*music->e3_gain[2] + music->e3_off[2])*((music_e3[3])*music->e3_gain[3] + music->e3_off[3]);
      //Float_t r3 = ((music_e3[4])*music->e3_gain[4] + music->e3_off[4])*((music_e3[5])*music->e3_gain[5] + music->e3_off[5]);  //added on 22.05.2018 SB
      //Float_t r4 = ((music_e3[6])*music->e3_gain[6] + music->e3_off[6])*((music_e3[7])*music->e3_gain[7] + music->e3_off[7]);  //added on 22.05.2018 SB
      
      if ( (r1 > 0) && (r2 > 0) )
    {
    //  b_de3 = kTRUE;
      de[2] = sqrt( sqrt(r1) * sqrt(r2) ) ;  // corrrected JSW 19.03.2010
      //de[2] = sqrt( sqrt( sqrt(r1) * sqrt(r2) ) * sqrt( sqrt(r3) * sqrt(r4) ) );   //changed on 22.05.2018 SB
      de_cor[2] = de[2];
    } 

//       if (music_b_t3[0] && music_b_t3[1] && music_b_t3[2] && music_b_t3[3]) 
//  b_dt3 = kTRUE;
    


    //hMUSIC1_MUSIC2->Fill(de[0],de[1]);

      /* Position (X) correction by TPC */      
      
     if(b_tpc_xy[4]&&b_tpc_xy[5]){

    music1_x_mean = tpc_music41_x;
    music2_x_mean = tpc_music42_x;
    music3_x_mean = tpc_music43_x;

    Float_t power, Corr;
    // correction for MUSIC41
    if(b_de1){
      power = 1., Corr = 0.;
      for(int i=0;i<4;i++) {
	Corr += music->pos_a1[i] * power;
	power *= music1_x_mean;
      }
      if (Corr!=0) {
	Corr = music->pos_a1[0] / Corr;
	de_cor[0] = de[0] * Corr;
      }
      
    }

    // correction for MUSIC42
    if(b_de2){
      power = 1., Corr = 0.;
      for(int i=0;i<4;i++) {
	Corr += music->pos_a2[i] * power;
	power *= music2_x_mean;
      }
      if (Corr!=0) {
	Corr = music->pos_a2[0] / Corr;
	de_cor[1] = de[1] * Corr;
      }
        }
    // correction for MUSIC43
    if(b_de3){
      power = 1., Corr = 0.;
      for(int i=0;i<4;i++) {
	Corr += music->pos_a3[i] * power;
	power *= music3_x_mean;
      }
      if (Corr!=0) {
	Corr = music->pos_a3[0] / Corr;
	de_cor[2] = de[2] * Corr;
      }
     
            }
        }    
    }
///==================================================================================///
                    /// Start of Scintillator  analysis
///==================================================================================///

   mhtdc_sc21lr_dt=0;
   mhtdc_sc31lr_dt=0;
   mhtdc_sc41lr_dt=0;
   mhtdc_sc42lr_dt=0;
   
  if(0!=tdc_sc21l[0] && 0!=tdc_sc21r[0]){
    mhtdc_sc21lr_dt = sci->mhtdc_factor_ch_to_ns*( rand3() + tdc_sc21l[0]  -  tdc_sc21r[0] );
    mhtdc_sc21lr_x  = mhtdc_sc21lr_dt * sci->mhtdc_factor_21l_21r + sci->mhtdc_offset_21l_21r;
  
    
    float sc21pos_from_tpc    = -999.9;
    if(b_tpc_xy[0]&&b_tpc_xy[1]){
      sc21pos_from_tpc =  tpc21_22_sc21_x ;
    }else if(b_tpc_xy[2]&&b_tpc_xy[3]){
      sc21pos_from_tpc =  tpc23_24_sc21_x ;
    }
  }

  if(0!=tdc_sc41l[0] && 0!=tdc_sc41r[0]){
    mhtdc_sc41lr_dt = sci->mhtdc_factor_ch_to_ns*( rand3() + tdc_sc41l[0]  -  tdc_sc41r[0] );
    mhtdc_sc41lr_x  = mhtdc_sc41lr_dt * sci->mhtdc_factor_41l_41r + sci->mhtdc_offset_41l_41r;
   }

  if(0!=tdc_sc42l[0] && 0!=tdc_sc42r[0]){
    mhtdc_sc42lr_dt = sci->mhtdc_factor_ch_to_ns*( rand3() + tdc_sc42l[0]  -  tdc_sc42r[0] );
    mhtdc_sc42lr_x  = mhtdc_sc42lr_dt * sci->mhtdc_factor_42l_42r + sci->mhtdc_offset_42l_42r;
   }

  if(0!=tdc_sc43l[0] && 0!=tdc_sc43r[0]){
    mhtdc_sc43lr_dt = sci->mhtdc_factor_ch_to_ns*( rand3() + tdc_sc43l[0]  -  tdc_sc43r[0] );
    mhtdc_sc43lr_x  = mhtdc_sc43lr_dt * sci->mhtdc_factor_43l_43r + sci->mhtdc_offset_43l_43r;
     }

  if(0!=tdc_sc31l[0] && 0!=tdc_sc31r[0]){
    mhtdc_sc31lr_dt = sci->mhtdc_factor_ch_to_ns*( rand3() + tdc_sc31l[0]  -  tdc_sc31r[0] );
    mhtdc_sc31lr_x  = mhtdc_sc31lr_dt * sci->mhtdc_factor_31l_31r + sci->mhtdc_offset_31l_31r;
   
  }

  if(0!=tdc_sc81l[0] && 0!=tdc_sc81r[0]){
    mhtdc_sc81lr_dt = sci->mhtdc_factor_ch_to_ns*( rand3() + tdc_sc81l[0]  -  tdc_sc81r[0] );
    mhtdc_sc81lr_x  = mhtdc_sc81lr_dt * sci->mhtdc_factor_81l_81r + sci->mhtdc_offset_81l_81r;
      }

  if(0!=tdc_sc21l[0] && 0!=tdc_sc21r[0] && 0!=tdc_sc41l[0] && 0!=tdc_sc41r[0]){
    mhtdc_tof4121 = sci->mhtdc_factor_ch_to_ns*( 0.5*(tdc_sc41l[0]+tdc_sc41r[0])  - 0.5*(tdc_sc21l[0]+tdc_sc21r[0]) ) + sci->mhtdc_offset_41_21;
    //if(bDrawHist) hMultiHitTDC_TOF_41_21->Fill(mhtdc_tof4121);
  }
  if(0!=tdc_sc21l[0] && 0!=tdc_sc21r[0] && 0!=tdc_sc42l[0] && 0!=tdc_sc42r[0]){
    mhtdc_tof4221 = sci->mhtdc_factor_ch_to_ns*( 0.5*(tdc_sc42l[0]+tdc_sc42r[0])  - 0.5*(tdc_sc21l[0]+tdc_sc21r[0]) ) + sci->mhtdc_offset_42_21;
   // if(bDrawHist) hMultiHitTDC_TOF_42_21->Fill(mhtdc_tof4221);
  }
  if(0!=tdc_sc21l[0] && 0!=tdc_sc21r[0] && 0!=tdc_sc43l[0] && 0!=tdc_sc43r[0]){
    mhtdc_tof4321 = sci->mhtdc_factor_ch_to_ns*( 0.5*(tdc_sc43l[0]+tdc_sc43r[0])  - 0.5*(tdc_sc21l[0]+tdc_sc21r[0]) ) + sci->mhtdc_offset_43_21;
   // if(bDrawHist) hMultiHitTDC_TOF_43_21->Fill(mhtdc_tof4321);
  }
  if(0!=tdc_sc21l[0] && 0!=tdc_sc21r[0] && 0!=tdc_sc31l[0] && 0!=tdc_sc31r[0]){
    mhtdc_tof3121 = sci->mhtdc_factor_ch_to_ns*( 0.5*(tdc_sc31l[0]+tdc_sc31r[0])  - 0.5*(tdc_sc21l[0]+tdc_sc21r[0]) ) + sci->mhtdc_offset_31_21;
    //if(bDrawHist) hMultiHitTDC_TOF_31_21->Fill(mhtdc_tof3121);
  }
  if(0!=tdc_sc21l[0] && 0!=tdc_sc21r[0] && 0!=tdc_sc81l[0] && 0!=tdc_sc81r[0]){
    mhtdc_tof8121 = sci->mhtdc_factor_ch_to_ns*( 0.5*(tdc_sc81l[0]+tdc_sc81r[0])  - 0.5*(tdc_sc21l[0]+tdc_sc21r[0]) ) + sci->mhtdc_offset_81_21;
  //  if(bDrawHist) hMultiHitTDC_TOF_81_21->Fill(mhtdc_tof8121);
  }



        
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
    
    

  /*  Raw data  */
   sci_l[2] = de_21l;  /* 21L         */
   sci_r[2] = de_21r;  /* 21R         */
   sci_tx[2] = dt_21l_21r + rand3();

   sci_l[5] = de_41l;  /* 41L         */
   sci_r[5] = de_41r;  /* 41R         */
   sci_tx[5] = dt_41l_41r + rand3();
  
   sci_l[6] = de_42l;  /* 42L         */
   
   sci_r[6] = de_42r;  /* 42R         */
   sci_tx[6] = dt_42l_42r + rand3();

   sci_l[7] = de_43l;  /* 43L         */
   sci_r[7] = de_43r;  /* 43R         */
   sci_tx[7] = dt_43l_43r + rand3();

   sci_l[10] = de_81l; /* 81L         */
   sci_r[10] = de_81r; /* 81R         */
   sci_tx[10] = dt_81l_81r + rand3();
   
   sci_tx[3] = dt_22l_22r + rand3();
   
    
    for (int cnt=0;cnt<7;cnt++) // 
     {
       int idx = 0 ;
       float posref =-999;
       //int mw_idx = 0;
       //Float_t mwx = 0;
       switch(cnt)
     {
     case 0:        /* SC21 */
       idx = 2; 
        // posref from tpc
	   if(b_tpc_xy[0]&&b_tpc_xy[1]){
	     posref =  tpc21_22_sc21_x ;
	   }else if(b_tpc_xy[2]&&b_tpc_xy[3]){
	     posref =  tpc23_24_sc21_x ;
	   }
       //mw_idx = 2;
       //mwx = mw_sc21_x;
       break;    
     case 1:        /* SC21 delayed */
       idx = 3; 
       if(b_tpc_xy[6]){ posref = tpc_x[6]; }
       //mw_idx = 2;
       //mwx = mw_sc21_x;
       break;    
     case 2:        /* SC41 */
       idx = 5; 
       if(b_tpc_xy[4]&&b_tpc_xy[5]){
	     posref =  tpc_sc41_x ;
	   }
       //mw_idx = 5;
       //mwx = tpc_sc41_x;
       break;    
     case 3:        /* SC42 */
           idx = 6;
       if(b_tpc_xy[4]&&b_tpc_xy[5]){
	     posref =  tpc_sc42_x ;
	   }
       break;
     case 4:
       idx = 7;     /* SC43 */
       if(b_tpc_xy[4]&&b_tpc_xy[5]){
	     posref =  tpc_sc43_x ;
	   }
       break;
     case 5:
       idx = 10;    /* SC81 */
       // no position reference from tpc
       break;
       
       case 6:
       idx = 3;    /* SC22 */
       if(b_tpc_xy[2]&&b_tpc_xy[3]){
         posref =  tpc23_24_sc22_x ;
       }else if(b_tpc_xy[0]&&b_tpc_xy[1]){
         posref =  tpc21_22_sc22_x ;
       }
       break;
     default: idx = 2;
     }   
     
       
       sci_b_l[idx] = Check_WinCond(sci_l[idx], cSCI_L);// cSCI_L[idx]->Test(sci_l[idx]);
       sci_b_r[idx] = Check_WinCond(sci_r[idx], cSCI_R);// cSCI_R[idx]->Test(sci_r[idx]);

       if(sci_b_l[idx] && sci_b_r[idx])
     {     
       sci_e[idx] = sqrt( (sci_l[idx] - sci->le_a[0][idx]) * sci->le_a[1][idx] 
                  * (sci_r[idx] - sci->re_a[0][idx]) * sci->re_a[1][idx]);
       
       sci_b_e[idx] = Check_WinCond(sci_e[idx], cSCI_E);// cSCI_E[idx]->Test(sci_e[idx]);
     }


       /*   Position in X direction:   */
       sci_b_tx[idx] = Check_WinCond(sci_tx[idx], cSCI_Tx);// cSCI_Tx[idx]->Test(sci_tx[idx]);
       if (sci_b_tx[idx])
     {

       /* mm-calibrated     */
       Float_t R = sci_tx[idx] ;//+ rand3(); 
       
       Float_t power = 1., sum = 0.;
       for(int i=0;i<7;i++)
         {
           sum += sci->x_a[i][idx] * power;
           power *= R;  
         }
       
       sci_x[idx] = sum;
       
       sci_b_x[idx] = Check_WinCond(sci_x[idx], cSCI_X);// cSCI_X[idx]->Test(sci_x[idx]);
     }

     } // end of loop for indices 2,3,5,6,7,10   


   /***  Scintillator Tof  spectra ***/

   /// S21 - S41 

   /*  Calibrated tof  */
   sci_tofll2 = dt_21l_41l*sci->tac_factor[2] - sci->tac_off[2] ;   /* S41L- S21L */
  
   sci_tofrr2 = dt_21r_41r*sci->tac_factor[3] - sci->tac_off[3] ;   /* S41R- S21R */
  
  
   sci_b_tofll2 = Check_WinCond(sci_tofll2, cSCI_LL2);// cSCI_TofLL2->Test(sci_tofll2);
   sci_b_tofrr2 = Check_WinCond(sci_tofrr2, cSCI_RR2);// cSCI_TofRR2->Test(sci_tofrr2);

    /* sum of Tof_LL and Tof_RR corrects for position in stop/start scint.      */
    if(sci_b_tofll2 && sci_b_tofrr2){
        /* TOF SC41 - SC21 [ps]  */
     sci_tof2        =   (sci->tof_bll2 * sci_tofll2 + sci->tof_a2 + sci->tof_brr2 * sci_tofrr2)/2.0 ;  // tof_a2  is essentially unnecessary (even confusing) = 0
     sci_tof2_calib   =  -1.0*sci_tof2 + id->id_tofoff2;
   }
   
   
    /*   S21 - S42 Calibrated tof  */
   sci_tofll3 = dt_42l_21l*sci->tac_factor[5] - sci->tac_off[5] ;   /* S42L- S21L */  // tac_off is essentially unnecessary (even confusing)
   sci_tofrr3 = dt_42r_21r*sci->tac_factor[6] - sci->tac_off[6] ;   /* S42R- S21R */  // tac_off is essentially unnecessary (even confusing)
  sci_b_tofll3 = Check_WinCond(sci_tofll3, cSCI_LL3);// cSCI_TofLL2->Test(sci_tofll2);
  sci_b_tofrr3 = Check_WinCond(sci_tofrr3, cSCI_RR3);// cSCI_TofRR2->Test(sci_tofrr2);
   if(sci_b_tofll3 && sci_b_tofrr3){
     sci_tof3        =   (sci->tof_bll3 * sci_tofll3 + sci->tof_a3 + sci->tof_brr3 * sci_tofrr3)/2.0 ;  // tof_a3  is essentially unnecessary (even confusing) = 0
     sci_tof3_calib   =  -1.0*sci_tof3 + id->id_tofoff3;
   }
  
  /*   S21 - S81 Calibrated tof  */
   sci_tofll4 = dt_21l_81l*sci->tac_factor[9] - sci->tac_off[9] ;     /* S81L- S21L */  // tac_off is essentially unnecessary (even confusing)
   sci_tofrr4 = dt_21r_81r*sci->tac_factor[10] - sci->tac_off[10] ;   /* S82R- S21R */  // tac_off is essentially unnecessary (even confusing)
 
   sci_b_tofll4 = Check_WinCond(sci_tofll4, cSCI_LL4);// cSCI_TofLL4->Test(sci_tofll4);
   sci_b_tofrr4 = Check_WinCond(sci_tofrr4, cSCI_RR4);// cSCI_TofRR4->Test(sci_tofrr4);

   // sum of Tof_LL and Tof_RR corrects for position in stop/start scint.      
   if (sci_b_tofll4 && sci_b_tofrr4)
     {      // TOF SC81 - SC21 [ps]
       sci_tof4 =  (sci->tof_bll4 * sci_tofll4 + sci->tof_a4 
            + sci->tof_brr4 * sci_tofrr4)/2.0 ;
      
     }
   
  /*   S22 - S41 Calibrated tof  */
   sci_tofll5 = dt_22l_41l*sci->tac_factor[12] - sci->tac_off[12] ;     /* S41L- S22L */  // tac_off is essentially unnecessary (even confusing)
   sci_tofrr5 = dt_22r_41r*sci->tac_factor[13] - sci->tac_off[13] ;     /* S41R- S22R */  // tac_off is essentially unnecessary (even confusing)
   
   sci_b_tofll5 = Check_WinCond(sci_tofll5, cSCI_LL5);// cSCI_TofLL4->Test(sci_tofll4);
   sci_b_tofrr5 = Check_WinCond(sci_tofrr5, cSCI_RR5);// cSCI_TofRR4->Test(sci_tofrr4);

  // sci_b_tofll5 = cSCI_TofLL5->Test(sci_tofll5);
  // sci_b_tofrr5 = cSCI_TofRR5->Test(sci_tofrr5);
   if(sci_b_tofll5 && sci_b_tofrr5){
     sci_tof5        =   (sci->tof_bll5 * sci_tofll5 + sci->tof_a5 + sci->tof_brr5 * sci_tofrr5)/2.0 ;  // tof_a5  is essentially unnecessary (even confusing) = 0
     sci_tof5_calib   =  -1.0*sci_tof5 + id->id_tofoff5;
     
    
   }
   /* check for polygon in raw detof spectrum of SC41 */
   
   sci_b_detof = Check_PolyCond_X_Y(sci_tof2, sci_e[5], cSCI_detof, 4); // cSCI_detof->Test(sci_tof2, sci_e[5]);
 ///==================================================================================///
                                /// Start of ID analysis
///==================================================================================///
  float speed_light = 0.299792458; //m/ns
  float temp_tm_to_MeV = 299.792458;
  float temp_mu = 931.4940954; //MeV
  
  // Extraction of position to be used for momentum analysis
  float temp_s8x = mhtdc_sc81lr_x;
  float temp_s4x = -999.;
  if(b_tpc_xy[4] && b_tpc_xy[5]){
    temp_s4x = tpc_x_s4 ;
  }
  float temp_s2x = -999.; //fill in the next if part
  if(1== id->mhtdc_s2pos_option){//SC21X from multihit tdc is used for S2X
    temp_s2x = mhtdc_sc21lr_x;
   // cout<<"mhtdc_sc21lr_x " <<mhtdc_sc21lr_x << endl;
  }
  if(2== id->mhtdc_s2pos_option){//TPCX is used
    if(b_tpc_xy[2] && b_tpc_xy[3]){//tpc2324
      temp_s2x = tpc_x_s2_foc_23_24;
    }else if(b_tpc_xy[1] && b_tpc_xy[3]){//tpc2224
      temp_s2x = tpc_x_s2_foc_22_24;
    }else if (b_tpc_xy[0] && b_tpc_xy[1]){//tpc2122
      temp_s2x = tpc_x_s2_foc_21_22;
    }
  }
 
  
  ////=======================================
  ////   S2S4 MultihitTDC ID analysis
  
  id_mhtdc_beta_s2s4=0;
  id_mhtdc_gamma_s2s4=0;
  float mean_brho_s2s4 =0;
  id_mhtdc_delta_s2s4=0;
  id_mhtdc_aoq_s2s4=0;
  id_mhtdc_v_cor_music41=0;
  id_mhtdc_z_music41=0;
  id_mhtdc_v_cor_music42=0;
  id_mhtdc_z_music42=0;
  
  
   // Calculation of velocity beta and gamma
  id_mhtdc_beta_s2s4   =  ( id->mhtdc_length_s2s4 / mhtdc_tof4121) / speed_light;
  id_mhtdc_gamma_s2s4  = 1./sqrt(1. - id_mhtdc_beta_s2s4*id_mhtdc_beta_s2s4);
  
  // calculation of delta(momentum_deviation) and AoQ
   mean_brho_s2s4 = 0.5*( frs->bfield[2] + frs->bfield[3] );
   
    if(temp_s4x==-999 || temp_s2x==-999){
        id_mhtdc_aoq_s2s4=0;
        id_mhtdc_aoq_corr_s2s4=0;
    }
        else if( -200<temp_s4x && temp_s4x<200. && -200.< temp_s2x && temp_s2x<200. ){
            id_mhtdc_delta_s2s4 = ( temp_s4x - (temp_s2x * frs->magnification[1] ))/(-1.0 * frs->dispersion[1] *1000.0 ) ; //1000 is dispertsion from meter to mm. -1.0 is sign definition.
   // cout<<"id_mhtdc_beta_s2s4 " << id_mhtdc_beta_s2s4 <<endl;
        if(0.0 < id_mhtdc_beta_s2s4 && id_mhtdc_beta_s2s4 < 1.0){
            id_mhtdc_aoq_s2s4 = mean_brho_s2s4 *( 1. + id_mhtdc_delta_s2s4   ) * temp_tm_to_MeV / (temp_mu * id_mhtdc_beta_s2s4 * id_mhtdc_gamma_s2s4);
            
            id_mhtdc_aoq_corr_s2s4 = id_mhtdc_aoq_s2s4 - id->a2AoQCorr * id_a4;
 
    }
    
  }
  else id_mhtdc_aoq_s2s4=0;
    // cout<<"id_mhtdc_aoq_s2s4 " << id_mhtdc_aoq_s2s4<<" mean_brho_s2s4 " <<mean_brho_s2s4 << " id_mhtdc_delta_s2s4 " << id_mhtdc_delta_s2s4 << " temp_tm_to_MeV " <<temp_tm_to_MeV<< " id_mhtdc_beta_s2s4 " <<id_mhtdc_beta_s2s4 << " id_mhtdc_gamma_s2s4  " <<id_mhtdc_gamma_s2s4 << endl;

  // calculation of dE and Z
  // from MUSIC41
  float temp_music41_de = de[0]>0.0;
  if( (temp_music41_de>0.0)  && (id_mhtdc_beta_s2s4>0.0) && (id_mhtdc_beta_s2s4<1.0)){
    Double_t power = 1., sum = 0.;
    for (int i=0;i<4;i++){
      sum += power * id->mhtdc_vel_a_music41[i];
      power *= id_mhtdc_beta_s2s4;
    }
    id_mhtdc_v_cor_music41 = sum;
    if (id_mhtdc_v_cor_music41 > 0.0){
      id_mhtdc_z_music41 = frs->primary_z * sqrt(de[0]/id_mhtdc_v_cor_music41 ) + id->mhtdc_offset_z_music41;
    }
  }

  float temp_music42_de = de[1]>0.0;
  if( (temp_music42_de>0.0)  && (id_mhtdc_beta_s2s4>0.0) && (id_mhtdc_beta_s2s4<1.0)){
    Double_t power = 1., sum = 0.;
    for (int i=0;i<4;i++){
      sum += power * id->mhtdc_vel_a_music42[i];
      power *= id_mhtdc_beta_s2s4;
    }
    id_mhtdc_v_cor_music42 = sum;
    if (id_mhtdc_v_cor_music42 > 0.0){
      id_mhtdc_z_music42 = frs->primary_z * sqrt(de[1]/id_mhtdc_v_cor_music42 ) + id->mhtdc_offset_z_music42;
    }
  }
 
 
  id_trigger=trigger;
  
   float mhtdc_gamma1square = 1.0 + TMath::Power(((299.792458/931.494)*(id_brho[0]/id_mhtdc_aoq_s2s4)),2);
    id_mhtdc_gamma_ta_s2 = TMath::Sqrt(mhtdc_gamma1square);
    id_mhtdc_dEdegoQ     =  (id_mhtdc_gamma_ta_s2  - id_mhtdc_gamma_s2s4)*id_mhtdc_aoq_s2s4;
    id_mhtdc_dEdeg       =  id_mhtdc_dEdegoQ * id_mhtdc_z_music41;
  // focal plane information
  // S2 priority: tpc2324 -> tpc2224 -> tpc2122 -> sc22 -> sc21
  if(b_tpc_xy[2] && b_tpc_xy[3]){//tpc2324
    id_x2 = tpc_x_s2_foc_23_24;
    id_y2 = tpc_y_s2_foc_23_24;
    id_a2 = tpc_angle_x_s2_foc_23_24;
    id_b2 = tpc_angle_y_s2_foc_23_24;
  }else if (b_tpc_xy[1] && b_tpc_xy[3]){//tpc2224
    id_x2 = tpc_x_s2_foc_22_24;
    id_y2 = tpc_y_s2_foc_22_24;
    id_a2 = tpc_angle_x_s2_foc_22_24;
    id_b2 = tpc_angle_y_s2_foc_22_24;
  }else if (b_tpc_xy[0] && b_tpc_xy[1]){//tpc2122
    id_x2 = tpc_x_s2_foc_21_22;
    id_y2 = tpc_y_s2_foc_21_22;
    id_a2 = tpc_angle_x_s2_foc_21_22;
    id_b2 = tpc_angle_y_s2_foc_21_22;
  }else if (sci_b_x[2]){//sc21
    id_x2 = sci_x[2];
    id_y2 = 0.0;
    id_a2 = 0.0;
    id_b2 = 0.0;
  }

  // S4 only 1 possibility =  TPC4142
  if(b_tpc_xy[4] && b_tpc_xy[5]){
    id_x4 = tpc_x_s4;
    id_a4 = tpc_angle_x_s4;
    id_y4 = tpc_y_s4;
    id_b4 = tpc_angle_y_s4;
  }

  // S8 only 1 possibility =  SC81x
  if( sci_b_x[10]){
    id_x8 = sci_x[10];
    id_a8 = 0.0;
    id_y8 = 0.0;
    id_b8 = 0.0;
  }

  /*  check that the positions are OK   */
  id_b_x2 = Check_WinCond(id_x2, cID_x2);// cID_x2->Test(id_x2);
  id_b_x4 = Check_WinCond(id_x4, cID_x4);// cID_x4->Test(id_x4);
  
 // id_b_x8 = Check_WinCond(id_x8, cID_x8);// cID_x4->Test(id_x4);
  

  // remove temporarily
  // hID_E_Xs4->Fill(id_x4,de[0]);// added by 2016Jun.16
  // hID_E_Xs2->Fill(id_x2,de[0]);// added by 2016Jun.16

  /*----------------------------------------------------------*/
  /* Determination of beta                                    */
  /* ID.TofOff(i)                   Flight time offset [ps]   */
  /* ID.Path(i)                     Flight path/c [ps]        */
  /* TOF(i)        BIN FLOAT(24),   Flight time  [ps]         */
  /*----------------------------------------------------------*/

if(id->tof_s4_select   == 1){
  ///SC21-SC41
  if (sci_b_tofll2 && sci_b_tofrr2){
    //// id_beta = id->id_path2 /(id->id_tofoff2 - sci_tof2);
    id_beta = id->id_path2 /  sci_tof2_calib ;// calculate non-inverted "real" tof already in sci analysis.
//  
  }
}

if(id->tof_s4_select   == 2){
  ///SC21-SC41
  if (sci_b_tofll3 && sci_b_tofrr3){
    //// id_beta = id->id_path2 /(id->id_tofoff2 - sci_tof2);
    id_beta = (id->id_path3 /  sci_tof3_calib) -0.1;// calculate non-inverted "real" tof already in sci analysis.
 
  }
}


if(id->tof_s4_select   == 3){

  ///SC22-SC41
    
  if (sci_b_tofll5 && sci_b_tofrr5){
    //// id_beta = id->id_path2 /(id->id_tofoff2 - sci_tof2);
    id_beta = (id->id_path5 /  sci_tof5_calib)/2 ;// calculate non-inverted "real" tof already in sci analysis.
//cout<<"id_beta " << id_beta<<endl;
   
  }
}

  // /*------------------------------------------------------*/
  // /* Determination of Brho                                */
  // /* Dispersion and magnification are still the same      */
  // /* variable for S41-S21 and S42-S41, adjust in setup.C  */
  // /*------------------------------------------------------*/

  // first half of FRS, TA-S2
  if (id_b_x2){
      id_rho[0]  = frs->rho0[0] * (1. - id_x2/1000./frs->dispersion[0]);
      id_brho[0] = (fabs(frs->bfield[0]) + fabs(frs->bfield[1]))/ 2. * id_rho[0];
      
//      if(bDrawHist){  hID_BRho[0]->Fill(id_brho[0]); }
  }
  // second half S2-S4
  
  if (id_b_x2 && id_b_x4){
      id_rho[1] = frs->rho0[1] * (1. - (id_x4 - frs->magnification[1] * id_x2) / 1000. / frs->dispersion[1]) ;
      id_brho[1] = (fabs(frs->bfield[2]) + fabs(frs->bfield[3]))/ 2. * id_rho[1];
     
      //cout<<"111 id_brho[1]  " << id_brho[1] << "frs->bfield[0]" << frs->bfield[2] << " frs->bfield[3] " << frs->bfield[3] << endl;
     // if(bDrawHist){  hID_BRho[1]->Fill(id_brho[1]); }
  }
  
  //cout<<"222 id_brho[1]  " << id_brho[1] << "frs->bfield[0]" << frs->bfield[2] << " frs->bfield[3] " << frs->bfield[3] << endl;
 // cout<<"id_rho[1] " << id_rho[1] << " frs->rho0[1] " << frs->rho0[1] << " id_x4 " << id_x4 << " frs->magnification[1] " << frs->magnification[1] << " id_x2 " << id_x2 << " frs->dispersion[1] " << frs->dispersion[1]<< endl;

 /*--------------------------------------------------------------*/
  /* Determination of A/Q                                         */
  /*--------------------------------------------------------------*/
  /* Beta(i)       BIN FLOAT(24),   Beta = v/c                    */
  /* Gamma(i)      BIN FLOAT(24),   Gamma= sqrt(1/1-beta**2)      */
  /*--------------------------------------------------------------*/
  Float_t f = 931.4940 / 299.792458 ;    /* factor needed for aoq calculation.. the u/(c*10^-6) factor  */

  /* for S2-S4 */
  if (sci_b_tofll2 && sci_b_tofrr2 && id_b_x2 && id_b_x4){
      if ((id_beta>0.0) && (id_beta<1.0)){
  	  id_gamma = 1./sqrt(1. - id_beta * id_beta);
  	  id_AoQ   = id_brho[1]/id_beta/id_gamma/f ;
  	 
      id_AoQ_corr = id_AoQ - id->a2AoQCorr * id_a4;  //correction for id_a2, JK 16.9.11
  	
  	  
  	  
  	  // if (!b_tpc_xy[4] || !b_tpc_xy[5]) // no sense to do "if TPC4142 dont work, correct with S4 angle ... "
  	  //   id_AoQ_corr = id_AoQ - id->a4AoQCorr * id_a4;
//   	  if(bDrawHist)
//   	    {
//   	      hID_AoQ->Fill(id_AoQ);
//   	      hID_AoQcorr->Fill(id_AoQ_corr);
// 	      hID_DeltaBrho_AoQ->Fill(id_AoQ,id_brho[0]-id_brho[1]);
//   	    }
  	  id_b_AoQ = kTRUE;
      }
  }

  /* for S2-S8  */
//   if (sci_b_tofll4 && sci_b_tofrr4 && id_b_x2 && id_b_x8){
//     if ((id_beta_s2s8>0.0) && (id_beta_s2s8<1.0)){
//       id_gamma_s2s8 = 1./sqrt(1. - id_beta_s2s8 * id_beta_s2s8);
//       id_AoQ_s2s8      = id_brho[2]/id_beta_s2s8/id_gamma_s2s8/ f ;
//       if(bDrawHist)
// 	{
// 	  hID_AoQ_s2s8->Fill(id_AoQ_s2s8);
// 	}
//       id_b_AoQ_s2s8 = kTRUE;
//     }
//   }


  /*------------------------------------------------*/
  /* Determination of Z                             */
  /*------------------------------------------------*/
  /****  S4  (MUSIC 1)   */
  //  if((de_cor[0]>0.0) && (id_beta>0.0) && (id_beta<1.0)) {

  if((de[0]>0.0) && (id_beta>0.0) && (id_beta<1.0))
    {
      Double_t power = 1., sum = 0.;
      for (int i=0;i<4;i++){
  	  sum += power * id->vel_a[i];
  	  power *= id_beta;
      }
      id_v_cor = sum;
      id_znocorr = frs->primary_z * sqrt(de[0]) + id->offset_z;
      
      if (id_v_cor > 0.0){
  	id_z = frs->primary_z * sqrt(de[0]/id_v_cor) + id->offset_z;
      }
      if ((id_z>0.0) && (id_z<100.0)){
	id_b_z = kTRUE;
      }
    }


  /****  S4  (MUSIC 2)   */
  
  if((de[1]>0.0) && (id_beta>0.0) && (id_beta<1.0)){
    Double_t power = 1., sum = 0.;
    for (int i=0;i<4;i++)
      {
	sum += power * id->vel_a2[i];
	power *= id_beta;
      }
    id_v_cor2 = sum;
      id_z2nocorr = frs->primary_z * sqrt(de[1]) + id->offset_z2;
    if (id_v_cor2 > 0.0){
      id_z2 = frs->primary_z * sqrt(de[1]/id_v_cor2) + id->offset_z2;
     
    }
    if ((id_z2>0.0) && (id_z2<100.0)){
      id_b_z2 = kTRUE;
    }
  }


  /****  S4  (MUSIC)   */
  if((de[2]>0.0) && (id_beta>0.0) && (id_beta<1.0)){
    Double_t power = 1., sum = 0.;
    for (int i=0;i<4;i++){
      sum += power * id->vel_a3[i];
      power *= id_beta;
    }
    id_v_cor3 = sum;
    if (id_v_cor3 > 0.0){
      id_z3 = frs->primary_z * sqrt(de[2]/id_v_cor3) + id->offset_z3;
    }

    if ((id_z3>0.0) && (id_z3<100.0)){
      id_b_z3 = kTRUE;
    }
  }


  //Z from sc81
//   if((sci_e[10]>0.0) && (id_beta_s2s8>0.0) && (id_beta_s2s8<1.0)){
//     Double_t power = 1., sum = 0.;
//     for (int i=0;i<4;i++){
//       sum += power * id->vel_a_sc81[i];
//       power *= id_beta_s2s8;
//     }
//     id_v_cor_sc81 = sum;
//     if (id_v_cor_sc81 > 0.0){
//       id_z_sc81 = frs->primary_z * sqrt(sci_e[10]/id_v_cor_sc81) + id->offset_z_sc81;
//     }
// 
//     if ((id_z_sc81>0.0) && (id_z_sc81<100.0)){
//       id_b_z_sc81 = kTRUE;
//     }
//   }

  //Z from tpc21222324
  double temp_s2tpc_de=1.0; int temp_s2tpc_de_count=0;
  for(int ii=0; ii<4; ii++){
    if(b_tpc_de[ii]){
      temp_s2tpc_de *=tpc_de[ii]; temp_s2tpc_de_count++;
    }
  }
  if(temp_s2tpc_de_count == 2){
    id_de_s2tpc = pow(temp_s2tpc_de, 1./temp_s2tpc_de_count);
    id_b_de_s2tpc = kTRUE;
  }

//   if(id_b_de_s2tpc && (id_beta_s2s8>0.0) && (id_beta_s2s8<1.0)){
//     Double_t power = 1., sum = 0.;
//     for (int i=0;i<4;i++){
//       sum += power * id->vel_a_s2tpc[i];
//       power *= id_beta_s2s8;
//     }
//     id_v_cor_s2tpc = sum;
//     if (id_v_cor_s2tpc > 0.0){
//       id_z_s2tpc = frs->primary_z * sqrt( id_de_s2tpc/id_v_cor_s2tpc ) + id->offset_z_s2tpc;
//       //printf("id_de_s2tpc = %f, v_cor = %f, z_s2tpc = %f\n",id_de_s2tpc, id_v_cor_s2tpc, id_z_s2tpc);
//     }
// 
//     if ((id_z_s2tpc>0.0) && (id_z_s2tpc<100.0)){
//       id_b_z_s2tpc = kTRUE;
//     }
//   }

  
  /// charge state for high Z
  /// S468 2020 April
  if(id_b_AoQ && id_b_x2 && id_b_z){
    float gamma1square = 1.0 + TMath::Power(((299.792458/931.494)*(id_brho[0]/id_AoQ)),2);
    id_gamma_ta_s2 = TMath::Sqrt(gamma1square);
    id_dEdegoQ     =  (id_gamma_ta_s2  - id_gamma)*id_AoQ;
    id_dEdeg       =  id_dEdegoQ * id_z;
   
//     if(bDrawHist){
//       hdEdegoQ_Z        ->Fill(id_z, id_dEdegoQ);
//       hdEdeg_Z          ->Fill(id_z, id_dEdeg);
//       if(cID_dEdeg_Z1->Test(id_z,id_dEdeg)){ // dEdeg_vs_Z check
//     hID_Z1_AoQ_cdEdegZ->Fill(id_AoQ, id_z);;
//     hID_Z1_AoQcorr_cdEdegZ->Fill(id_AoQ_corr, id_z);
//    
//     if(id_b_z && id_b_z2 && fabs(id_z2-id_z)<0.4){ // in addition Z same
//       hID_Z1_AoQ_zsame_cdEdegZ->Fill(id_AoQ, id_z);
//       hID_Z1_AoQcorr_zsame_cdEdegZ->Fill(id_AoQ_corr, id_z);
//             }
//       }
//     }
  }
  
  
   /*------------------------------------------------*/
   /* Identification Plots                           */
   /*------------------------------------------------*/

   /****  for S2-S4  ****/

   for (int i=0;i<5;i++){
     id_b_x4AoQ[i] = Check_PolyCond_Multi_X_Y(id_AoQ, id_x4, cID_x4AoQ, 5, 0);
     id_b_x2AoQ[i] = Check_PolyCond_Multi_X_Y(id_AoQ, id_x2, cID_x2AoQ, 5, 0);
     id_b_z_AoQ[i] =  Check_PolyCond_Multi_X_Y(id_AoQ, id_z, cID_Z_AoQ, 5, i);
   }
}

    
void FRS_Detector_System::get_Event_data(Raw_Event* RAW){
        
    RAW->set_DATA_MUSIC(de, de_cor, music_e1, music_e2, music_t1, music_t2); 
    RAW->set_DATA_SCI(sci_l, sci_r, sci_e, sci_tx, sci_x); 
    RAW->set_DATA_SCI_dT(dt_21l_21r, dt_41l_41r, dt_21l_41l, dt_21r_41r, dt_42l_42r, dt_43l_43r,dt_42l_21l, dt_42r_21r, dt_81l_81r, dt_21l_81l, dt_21r_81r);
    RAW->set_DATA_SCI_ToF(sci_tofll2, sci_tofll3, sci_tof2_calib, sci_tofrr2, sci_tofrr3, sci_tof3_calib, sci_tofll5, sci_tofrr5, sci_tof5_calib);
    RAW->set_DATA_TPC(tpc_lt,tpc_rt, tpc_x, tpc_y, x0, x1);
    RAW->set_DATA_ID_2_4(id_x2, id_y2, id_a2, id_b2, id_x4, id_y4, id_a4, id_b4);
    RAW->set_DATA_ID_Beta_Rho(id_brho, id_rho, id_beta, id_beta3, id_gamma);
    RAW->set_DATA_ID_Z_AoQ(id_AoQ, id_AoQ_corr, id_z, id_z2, id_z3,id_znocorr, id_z2nocorr, id_dEdeg,id_dEdegoQ);
    RAW->set_DATA_ID_Timestamp(timestamp, ts, ts2);
    RAW->set_DATA_FRS_SCALERS(time_in_ms,spill_count,ibin_for_s,ibin_for_100ms,ibin_for_spill,ibin_clean_for_s,ibin_clean_for_100ms, ibin_clean_for_spill,increase_scaler_temp);
    RAW->set_DATA_VFTX(TRaw_vftx_21l,TRaw_vftx_21r,TRaw_vftx_22l,TRaw_vftx_22r,TRaw_vftx_41l,TRaw_vftx_41r,TRaw_vftx_42l,TRaw_vftx_42r,TRaw_vftx);

     RAW->set_DATA_ID_MHTDC(id_mhtdc_aoq_s2s4, id_mhtdc_aoq_corr_s2s4, id_mhtdc_z_music41, id_mhtdc_z_music42, id_mhtdc_dEdeg, id_mhtdc_dEdegoQ, id_mhtdc_beta_s2s4);
}

//---------------------------------------------------------------

int* FRS_Detector_System::get_pdata(){return pdata;}

//---------------------------------------------------------------

Float_t FRS_Detector_System::rand3(){
    
  return random3.Uniform(-0.5,0.5);
    
}

//---------------------------------------------------------------

Int_t FRS_Detector_System::getbits(Int_t value, int nword, int start, int length){
    
  UInt_t buf = (UInt_t) value;
  buf = buf >> ((nword-1)*16 + start - 1);
  buf = buf & ((1 << length) - 1);
  return buf;
  
}

//---------------------------------------------------------------

Int_t FRS_Detector_System::get2bits(Int_t value, int nword, int start, int length){
    
  UInt_t buf = (UInt_t) value;
  buf = buf >> ((nword-1)*32 + start - 1);
  buf = buf & ((1 << length) - 1);
  return buf;
  
}

//---------------------------------------------------------------
Bool_t FRS_Detector_System::Check_WinCond(Float_t P, Float_t* V){
    
    if(P >= V[0] && P <= V[1]) return true;
    else return false;

}

//---------------------------------------------------------------
Bool_t FRS_Detector_System::Check_WinCond_Multi(Float_t P, Float_t** V, int cond_num){
    
    if(P >= V[cond_num][0] && P <= V[cond_num][1]) return true;
    else return false;

}

//---------------------------------------------------------------
Bool_t FRS_Detector_System::Check_WinCond_Multi_Multi(Float_t P, Float_t*** V, int cond_num, int cond_num_2){
    
    if(P >= V[cond_num][cond_num_2][0] && P <= V[cond_num][cond_num_2][1]) return true;
    else return false;

}

//---------------------------------------------------------------

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
Bool_t FRS_Detector_System::Check_PolyCond(Float_t* P, Float_t** V, int n ){
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
       if (((V[i][1] <= P[1]) && (V[i+1][1] > P[1]))     // an upward crossing
        || ((V[i][1] > P[1]) && (V[i+1][1] <=  P[1]))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(P[1]  - V[i][1]) / (V[i+1][1] - V[i][1]);
            if (P[0] <  V[i][0] + vt * (V[i+1][0] - V[i][0])) // P.x < intersect
                 ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    
    if((cn&1) == 0) return false;
    if((cn&1) == 1) return true;
    else return false;

    //return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}

//---------------------------------------------------------------

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
Bool_t FRS_Detector_System::Check_PolyCond_Multi(Float_t* P, Float_t*** V, int n, int cond_num )
{
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
       if (((V[cond_num][i][1] <= P[1]) && (V[cond_num][i+1][1] > P[1]))     // an upward crossing
        || ((V[cond_num][i][1] > P[1]) && (V[cond_num][i+1][1] <=  P[1]))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(P[1]  - V[cond_num][i][1]) / (V[cond_num][i+1][1] - V[cond_num][i][1]);
            if (P[0] <  V[cond_num][i][0] + vt * (V[cond_num][i+1][0] - V[cond_num][i][0])) // P.x < intersect
                 ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    
    if((cn&1) == 0) return false;
    if((cn&1) == 1) return true;
    else return false;

    //return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}
//---------------------------------------------------------------

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
Bool_t FRS_Detector_System::Check_PolyCond_X_Y(Float_t X, Float_t Y, Float_t** V, int n ){
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
       if (((V[i][1] <= Y) && (V[i+1][1] > Y))     // an upward crossing
        || ((V[i][1] > Y) && (V[i+1][1] <=  Y))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(Y  - V[i][1]) / (V[i+1][1] - V[i][1]);
            if (X <  V[i][0] + vt * (V[i+1][0] - V[i][0])) // P.x < intersect
                 ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    
    if((cn&1) == 0) return false;
    if((cn&1) == 1) return true;
    else return false;

    //return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}

//---------------------------------------------------------------

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
Bool_t FRS_Detector_System::Check_PolyCond_Multi_X_Y(Float_t X, Float_t Y, Float_t*** V, int n, int cond_num )
{
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
       if (((V[cond_num][i][1] <= Y) && (V[cond_num][i+1][1] > Y))     // an upward crossing
        || ((V[cond_num][i][1] > Y) && (V[cond_num][i+1][1] <=  Y))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(Y  - V[cond_num][i][1]) / (V[cond_num][i+1][1] - V[cond_num][i][1]);
            if (X <  V[cond_num][i][0] + vt * (V[cond_num][i+1][0] - V[cond_num][i][0])) // P.x < intersect
                 ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    
    if((cn&1) == 0) return false;
    if((cn&1) == 1) return true;
    else return false;

    //return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}

//---------------------------------------------------------------

void FRS_Detector_System::Setup_Conditions(){
        
    string line;
    int line_number = 0;
    
    const char* format = "%f %f %f %f %f %f %f %f %f %f %f %f %f %f";

    ifstream cond_a("Configuration_Files/FRS/FRS_Window_Conditions/lim_csum.txt");
    
    lim_csum = new Float_t**[4];
    
    for(int i = 0; i < 4; ++i){
    
    lim_csum[i] = new Float_t*[7];
    
        for(int j = 0; j < 7; ++j){
        
        
        lim_csum[i][j] = new Float_t[2];
        
        }
    }
    
    
    while(/*cond_a.good()*/getline(cond_a,line,'\n')){
    
    //getline(cond_a,line,'\n');
    if(line[0] == '#') continue;
    

        sscanf(line.c_str(),format,&lim_csum[line_number][0][0],&lim_csum[line_number][0][1]
                    ,&lim_csum[line_number][1][0],&lim_csum[line_number][1][1]
                    ,&lim_csum[line_number][2][0],&lim_csum[line_number][2][1]
                    ,&lim_csum[line_number][3][0],&lim_csum[line_number][3][1]
                    ,&lim_csum[line_number][4][0],&lim_csum[line_number][4][1]
                    ,&lim_csum[line_number][5][0],&lim_csum[line_number][5][1]
                    ,&lim_csum[line_number][6][0],&lim_csum[line_number][6][1]);
                    
    line_number++;
    
    
    }
    
   /*lim_csum = {{{300,1800},{350,1800},{300,1800},{300,1700},{300,2000},{300,2000},{300,2000}},
                 {{300,1800},{300,1800},{300,1800},{300,1700},{300,2000},{300,2000},{300,2000}},
                 {{300,1800},{300,1840},{300,1800},{300,1700},{300,2000},{300,2000},{300,2000}},
                 {{300,1800},{300,1810},{300,1800},{300,1700},{300,2000},{300,2000},{300,2000}}};
    */           
                 

    
    lim_xsum = new Float_t*[13];
    lim_ysum = new Float_t*[13];
    
    for(int i = 0; i < 13; ++i){
    
    
    lim_xsum[i] = new Float_t[2];
    lim_ysum[i] = new Float_t[2];
    
    //lim_xsum[i][0] = 1;
    //lim_xsum[i][1] = 8000;
    //lim_ysum[i][0] = 2;
    //lim_ysum[i][1] = 8000;
    
    
    }
    
    line_number = 0;
    
    format = "%f %f";

    ifstream cond_b("Configuration_Files/FRS/FRS_Window_Conditions/lim_xsum.txt");
    
    while(/*cond_b.good()*/getline(cond_b,line,'\n')){
    
    //getline(cond_b,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&lim_xsum[line_number][0],&lim_xsum[line_number][1]);
                    
    line_number++;
    }
    
    line_number = 0;
    
    format = "%f %f";

     ifstream cond_c("Configuration_Files/FRS/FRS_Window_Conditions/lim_ysum.txt");
    
    while(/*cond_c.good()*/getline(cond_c,line,'\n')){
    
    //getline(cond_c,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&lim_ysum[line_number][0],&lim_ysum[line_number][1]);
                    
    line_number++;
    }
    
    /*lim_xsum = {  {1,8000},  // MW11
                {1,8000},  // MW21
                {1,8000},  // MW22
                {1,8000},  // MW31
                {1,8000},  // MW41
                {1,8000},  // MW42
                {1,8000},  // MW51
                {1,8000},  // MW61
                {1,8000},  // MW71
                {1,8000},  // MW81
                {1,8000},  // MW82
                {1,8000},  // MWB1
                {1,8000}};  // MWB2

    lim_ysum = { {2,8000},  // MW11
                {2,8000},  // MW21
                {2,8000},  // MW22
                {2,8000},  // MW31
                {2,8000},  // MW41
                {2,8000},  // MW42
                {2,8000},  // MW51
                {2,8000},  // MW61
                {2,8000},  // MW71
                {2,8000},  // MW81
                {2,8000},  // MW82
                {2,8000},  // MWB1
                {2,8000}};  // MWB2*/
                
    
    
    /*for(int i=0;i<7;i++){
    
    //    cTPC_CSUM[i][0]=MakeWindowCond("TPC/CSUM1",condname,lim_csum1[i][0],
    //                 lim_csum1[i][1],"CSUM1");
    char condname[40];
    
    sprintf(condname,"TPC/CSUM1%s",tpc_name_ext1[i]);
    cTPC_CSUM[i][0]=MakeWinCond(condname,lim_csum1[i][0],lim_csum1[i][1],"CSUM1");
    
    sprintf(condname,"CSUM2%s",tpc_name_ext1[i]);
    cTPC_CSUM[i][1]=MakeWindowCond("TPC/CSUM2",condname,lim_csum2[i][0],lim_csum2[i][1],"CSUM2");
    
    sprintf(condname,"CSUM3%s",tpc_name_ext1[i]);
    cTPC_CSUM[i][2]=MakeWindowCond("TPC/CSUM3",condname,lim_csum3[i][0],lim_csum3[i][1],"CSUM3");
    
    sprintf(condname,"CSUM4%s",tpc_name_ext1[i]);
    cTPC_CSUM[i][3]=MakeWindowCond("TPC/CSUM4",condname,lim_csum4[i][0],lim_csum4[i][1],"CSUM4");
    }
    
      for(int i=0;i<13;i++){  // up to MW31
       //up to MWB2 (09.07.2018)

      char condname[40];
      sprintf(condname,"XSUM%s",mw_name_ext[i]);
      cMW_XSUM[i] = MakeWindowCond("MW/XSUM", condname, lim_xsum[i][0], lim_xsum[i][1], MW_XSUM);

      sprintf(condname,"YSUM%s",mw_name_ext[i]);
      cMW_YSUM[i] = MakeWindowCond("MW/YSUM", condname, lim_ysum[i][0], lim_ysum[i][1], MW_YSUM);
    }*/
    
    /*** MUSIC Conditions ***/
    
    cMusic1_E = new Float_t*[8];
    cMusic1_T = new Float_t*[8];
    cMusic2_E = new Float_t*[8];
    cMusic2_T = new Float_t*[8];
    cMusic3_T = new Float_t*[4];
    cMusic3_E = new Float_t*[4];

    cMusic3_dec = new Float_t[2];
    
    for(int i = 0; i < 8; ++i){
    
    
    cMusic1_E[i] = new Float_t[2];
    cMusic1_T[i] = new Float_t[2];
    cMusic2_E[i] = new Float_t[2];
    cMusic2_T[i] = new Float_t[2];

    if(i < 4){
        cMusic3_E[i] = new Float_t[2];
        cMusic3_T[i] = new Float_t[2];
    }
    
    }
    
    
    line_number = 0;
    
    format = "%f %f %f %f";

     ifstream cond_d("Configuration_Files/FRS/FRS_Window_Conditions/MUSIC1.txt");
    
    while(/*cond_d.good()*/getline(cond_d,line,'\n')){
    
    //getline(cond_d,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cMusic1_E[line_number][0],&cMusic1_E[line_number][1],&cMusic1_T[line_number][0],&cMusic1_T[line_number][1]);
                    
    line_number++;
    }

    line_number = 0;
    
    format = "%f %f %f %f";

     ifstream cond_e("Configuration_Files/FRS/FRS_Window_Conditions/MUSIC2.txt");
    
    while(/*cond_e.good()*/getline(cond_e,line,'\n')){
    
    //getline(cond_e,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cMusic2_E[line_number][0],&cMusic2_E[line_number][1],&cMusic2_T[line_number][0],&cMusic2_T[line_number][1]);
                    
    line_number++;
    }

    line_number = 0;
    
    format = "%f %f %f %f";

     ifstream cond_f("Configuration_Files/FRS/FRS_Window_Conditions/MUSIC3.txt");
    
    while(/*cond_f.good()*/getline(cond_f,line,'\n')){
    
    //getline(cond_f,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cMusic3_E[line_number][0],&cMusic3_E[line_number][1],&cMusic3_T[line_number][0],&cMusic3_T[line_number][1]);
                    
    line_number++;
    }

    
    line_number = 0;
    
    format = "%f %f";

     ifstream cond_g("Configuration_Files/FRS/FRS_Window_Conditions/MUSIC_dEc3.txt");
    
    while(/*cond_g.good()*/getline(cond_g,line,'\n')){
    
    //getline(cond_g,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cMusic3_dec[0],&cMusic3_dec[1]);
    }

    /***SCINTILATOR Condtions***/
    
    cSCI_L = new Float_t[2];
    cSCI_R = new Float_t[2];
    cSCI_E = new Float_t[2];
    cSCI_Tx = new Float_t[2];
    cSCI_X = new Float_t[2];
    
    
    line_number = 0;
    
    format = "%f %f";

     ifstream cond_h("Configuration_Files/FRS/FRS_Window_Conditions/SCI_Cons.txt");
    
    while(/*cond_h.good()*/getline(cond_h,line,'\n')){
    
    //getline(cond_h,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cSCI_L[0],&cSCI_L[1]);
        
    getline(cond_h,line,'\n');
        sscanf(line.c_str(),format,&cSCI_R[0],&cSCI_R[1]);
        
    getline(cond_h,line,'\n');
        sscanf(line.c_str(),format,&cSCI_E[0],&cSCI_E[1]);
        
    getline(cond_h,line,'\n');
        sscanf(line.c_str(),format,&cSCI_Tx[0],&cSCI_Tx[1]);
        
    getline(cond_h,line,'\n');
        sscanf(line.c_str(),format,&cSCI_X[0],&cSCI_X[1]);
        
    }
    
    cSCI_LL2 = new Float_t[2];
    cSCI_RR2 = new Float_t[2];
    cSCI_LL3 = new Float_t[2];
    cSCI_RR3 = new Float_t[2];
    cSCI_LL4 = new Float_t[2];
    cSCI_RR4 = new Float_t[2];
    
    cSCI_LL5 = new Float_t[2];
    cSCI_RR5 = new Float_t[2];
    format = "%f %f";

    ifstream cond_i("Configuration_Files/FRS/FRS_Window_Conditions/SCI_LLRR.txt");

    while(/*cond_i.good()*/getline(cond_i,line,'\n')){
    
    
    //getline(cond_i,line,'\n');
    if(line[0] == '#') continue;    
        sscanf(line.c_str(),format,&cSCI_LL2[0],&cSCI_LL2[1]);
    getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_RR2[0],&cSCI_RR2[1]);
        
    getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_LL3[0],&cSCI_LL3[1]);
        
    getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_RR3[0],&cSCI_RR3[1]);
        
        
    getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_LL4[0],&cSCI_LL4[1]);
        
    getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_RR4[0],&cSCI_RR4[1]);
        
        getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_LL5[0],&cSCI_LL5[1]);
        
        getline(cond_i,line,'\n');
        sscanf(line.c_str(),format,&cSCI_RR5[0],&cSCI_RR5[1]);
        
    }
    
    cSCI_detof = new Float_t*[5];
    
    for(int i = 0; i < 5; ++i){
    
    cSCI_detof[i] = new Float_t[2];
    
    }
    
    line_number = 0;
    
    format = "%f %f";

     ifstream cond_j("Configuration_Files/FRS/FRS_Window_Conditions/SCI_dEToF.txt");
    
    while(/*cond_j.good()*/getline(cond_j,line,'\n')){
    
    //getline(cond_j,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cSCI_detof[line_number][0],&cSCI_detof[line_number][1]);
                    
    line_number++;
    }

    /***ID Condtions***/
    
    cID_x2 = new Float_t[2];
    cID_x4 = new Float_t[2];
   // cID_x8 = new Float_t[2];
    cID_Z_Z = new Float_t[2];
    
    format = "%f %f";

     ifstream cond_k("Configuration_Files/FRS/FRS_Window_Conditions/ID_x2.txt");


    while(/*cond_k.good()*/getline(cond_k,line,'\n')){
    
    //getline(cond_k,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cID_x2[0],&cID_x2[1]);
        
    }
    
     ifstream cond_l("Configuration_Files/FRS/FRS_Window_Conditions/ID_x4.txt");

    while(/*cond_l.good()*/getline(cond_l,line,'\n')){
    
    //getline(cond_l,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cID_x4[0],&cID_x4[1]);
        
    }
    
//      ifstream cond_l("Configuration_Files/FRS_Window_Conditions/ID_x8.txt");
// 
//     while(/*cond_l.good()*/getline(cond_l,line,'\n')){
//     
//     //getline(cond_l,line,'\n');
//     if(line[0] == '#') continue;
//         sscanf(line.c_str(),format,&cID_x8[0],&cID_x8[1]);
//         
//     }
    
     ifstream cond_m("Configuration_Files/FRS/FRS_Window_Conditions/ID_Z_Z.txt");

    while(/*cond_m.good()*/getline(cond_m,line,'\n')){
    
    //getline(cond_m,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cID_Z_Z[0],&cID_Z_Z[1]);
        
    }
    
    
    cID_x4AoQ_Z = new Float_t*[5];
    
    for(int i = 0; i < 5; ++i){
    
    cID_x4AoQ_Z[i] = new Float_t[2];
    
    }
    
    
    line_number = 0;

     ifstream cond_n("Configuration_Files/FRS/FRS_Window_Conditions/ID_Z_Z.txt");

    while(/*cond_n.good()*/getline(cond_n,line,'\n')){
    
    //getline(cond_n,line,'\n');
    if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&cID_x4AoQ_Z[line_number][0],&cID_x4AoQ_Z[line_number][1]);
        
    line_number++;

    }
    
    cID_x2AoQ = new Float_t**[6];
    cID_x4AoQ = new Float_t**[6];
    cID_Z_AoQ = new Float_t**[5];
    
    for(int i = 0; i < 6; ++i){
    
    cID_x2AoQ[i] = new Float_t*[6];
    cID_x4AoQ[i] = new Float_t*[6];
    if (i < 5) cID_Z_AoQ[i] = new Float_t*[6];
    
        for(int j = 0; j < 6; ++j){
        
        cID_x2AoQ[i][j] = new Float_t[2];
        cID_x4AoQ[i][j] = new Float_t[2];
        if (i < 5) cID_Z_AoQ[i][j] = new Float_t[2];
        
        }
    
    }
    
    line_number = 0;
    int selection_number = 0;

     ifstream cond_o("Configuration_Files/FRS/FRS_Window_Conditions/ID_x2AoQ.txt");

    while(/*cond_o.good()*/getline(cond_o,line,'\n')){
    
    //getline(cond_o,line,'\n');
    if(line[0] == '#') continue;
    if(line[0] == '&'){
        selection_number++;
        line_number = 0;
        continue;
    }
    
        sscanf(line.c_str(),format,&cID_x2AoQ[selection_number][line_number][0],&cID_x2AoQ[selection_number][line_number][1]);
        
    line_number++;

    }
    
    line_number = 0;
    int selection_number_2 = 0;

     ifstream cond_z("Configuration_Files/FRS/FRS_Window_Conditions/ID_x4AoQ.txt");

    while(/*cond_b.good()*/getline(cond_z,line,'\n')){
    
    //getline(cond_o,line,'\n');
    if(line[0] == '#') continue;
    if(line[0] == '&'){
        selection_number_2++;
        line_number = 0;
        continue;
    }
    
        sscanf(line.c_str(),format,&cID_x4AoQ[selection_number_2][line_number][0],&cID_x4AoQ[selection_number_2][line_number][1]);
        
    line_number++;

    }
    
    
    line_number = 0;
    selection_number = 0;

     ifstream cond_p("Configuration_Files/FRS/FRS_Window_Conditions/ID_Z_AoQ.txt");

    while(/*cond_p.good()*/getline(cond_p,line,'\n')){
    
    //getline(cond_p,line,'\n');
    if(line[0] == '#') continue;
    if(line[0] == '&'){
        selection_number++;
        line_number = 0;
        continue;
    }
        sscanf(line.c_str(),format,&cID_Z_AoQ[selection_number][line_number][0],&cID_Z_AoQ[selection_number][line_number][1]);
        
    line_number++;

    }
    
    cID_dEToF = new Float_t*[5];

    for(int i = 0; i < 5; ++i){
    
    cID_dEToF[i] = new Float_t[2];
    
    }
    line_number = 0;
    selection_number = 0;

     ifstream cond_q("Configuration_Files/FRS/FRS_Window_Conditions/ID_dEToF.txt");

    while(/*cond_q.good()*/getline(cond_q,line,'\n')){
    
    //getline(cond_p,line,'\n');
    if(line[0] == '#') continue;

        sscanf(line.c_str(),format,&cID_dEToF[line_number][0],&cID_dEToF[line_number][1]);
        
    line_number++;

    }
    
    
    /*for(int i=0;i<8;i++){
    
       sprintf(name,"Music1_E(%f)",i);
       cMusic1_E[i] = MakeWindowCond("MUSIC/MUSIC(1)/E",name, 10, 4086, "hMUSIC1_E");

        sprintf(name,"Music1_T(%d)",i);
       cMusic1_T[i] = MakeWindowCond("MUSIC/MUSIC(1)/T",name, 10, 4086,"hMUSIC1_T");

       sprintf(name,"Music2_E(%d)",i);
       cMusic2_E[i] = MakeWindowCond("MUSIC/MUSIC(2)/E",name, 10, 4086, "hMUSIC2_E");

       sprintf(name,"Music2_T(%d)",i);
       cMusic2_T[i] = MakeWindowCond("MUSIC/MUSIC(2)/T",name, 10, 4086, "hMUSIC2_T");
     }

   for(int i=0;i<4;i++){
       
       sprintf(name,"Music3_E(%d)",i);
       cMusic3_E[i] = MakeWindowCond("MUSIC/MUSIC(3)/E",name, 10, 4086, "hMUSIC3_E");

       sprintf(name,"Music3_T(%d)",i);
       cMusic3_T[i] = MakeWindowCond("MUSIC/MUSIC(3)/T",name, 10, 4086, "hMUSIC3_T");
    }
   
   
    cMusic3_dec = MakeWindowCond("MUSIC/MUSIC 3","Music3_dec", 10., 4000., "hMUSIC3_dECOR");*/

}


#ifdef CALIBRATION_VFTX
void FRS_Detector_System::VFTX_Calibration(int module, int channel) {
    cout<<"!!!!VFTX CALIBRATION!!!"<<endl;
  int l_VFTX_SN[VFTX_N] = VFTX_SN;
  double integral_ft;
  //double integral_tot_ft;
  int    FirstBin_ft;
  double integral [1000];
  double calib[1000];
  for (int bin=1; bin<=1000; bin++) {
    calib[bin-1] = 0. ;
    integral [bin-1] = 0. ;
  }
  std::ofstream f_int;
  std::ofstream f_cal;
  f_cal.open(Form("Configuration_Files/FRS/VFTX_Calib/VFTX_%05d_Bin2Ps_ch%02d.dat",l_VFTX_SN[module],channel),std::ios::out);
  f_int.open(Form("Configuration_Files/FRS/VFTX_Int/VFTX_%05d_INT_ch%02d.txt",l_VFTX_SN[module],channel),std::ios::out);
  // calculation of the integral
  //integral_tot_ft = (double)h1_vftx_ft[module][channel]->Integral();
  integral_ft     = (double)h1_vftx_ft[module][channel]->Integral(5,1000);
  cout<<"integral_ft "<<integral_ft<<" module " << module << " channel " << channel << endl;
  // calculation of the first bin
  FirstBin_ft = 1 ;
  for (int bin=2; bin<=1000; bin++)     {
    if (h1_vftx_ft[module][channel]->GetBinContent(bin)>0) {
      FirstBin_ft = bin;
      break;
    }
  }//end of for(bin) 
  // calculation of the calibration parameters
  if(integral_ft>0){
    for (int bin=FirstBin_ft; bin<=1000; bin++) {
      integral[bin-1]  = integral[bin-2]+h1_vftx_ft[module][channel]->GetBinContent(bin);
      calib[bin-1]     = (integral[bin-1]*5000.)/integral_ft ;
    }
  }//end of if spectrum not empty
  //SAVE
  for(int bin=1; bin<=1000;bin++)       {
    f_cal <<std::setw(4) <<bin-1 <<" " <<std::setw(15) <<std::fixed <<std::setprecision(9) <<calib[bin-1]    <<std::endl;
    f_int <<std::setw(4) <<bin-1 <<" " <<std::setw(15) <<std::fixed <<std::setprecision(9) <<integral[bin-1] <<std::endl;
  }
  f_cal.close();
  f_int.close();
  return;
}
#endif

// --- ----------------------------------- --- //
// --- RECONSTRUCT Traw FROM THE VFTX DATA --- //
// --- ----------------------------------- --- //
void FRS_Detector_System::m_VFTX_Bin2Ps(){
  int l_VFTX_SN[VFTX_N]  = VFTX_SN;

  for(int i=0; i<VFTX_N; i++)
    for(int j=0; j<VFTX_MAX_CHN; j++)
      for(int k=0; k<1000; k++)
	VFTX_Bin2Ps[i][j][k]=0.;

  int b; double ft_ps; 
  for(int mod=0; mod<VFTX_N; mod++){
    for(int ch=0; ch<VFTX_MAX_CHN; ch++){
      std::ifstream in;
      in.open(Form("Configuration_Files/FRS/VFTX_Calib/VFTX_%05d_Bin2Ps_ch%02d.dat",l_VFTX_SN[mod],ch));
      if(!in.is_open()){
	for(int bin=0; bin<1000; bin++)
	  VFTX_Bin2Ps[mod][ch][bin] = 0.; // no data in ps if we don't have the files
	printf("WARNING : VFTX %05d ch %02d file not found, you will not have precise data\n",l_VFTX_SN[mod],ch);
      }
      else {
	while(!in.eof()) {
	  in >>b >>ft_ps;
	  VFTX_Bin2Ps[mod][ch][b] = ft_ps;
	  if(b>1000) printf(" !!! WARNING !!!! file Configuration_Files/FRS/VFTX_Calib/VFTX%02d_Bin2Ps_ch%2d.dat, overflow b=%i \n",mod,ch,b);
	}// end of for (ch over VFTX_CHN)
	in.close();
      }
    }//end of for (ch)
  }// end of for (mod over VFTX_N
  return;
}

Double_t FRS_Detector_System::VFTX_GetTraw_ps(int module, int channel, int cc, int ft, float rand) {
  Double_t gain;
  Double_t calib = (Double_t)VFTX_Bin2Ps[module][channel][ft];
 
  if (calib==0) calib = ft;
  if(rand<0) {
    Double_t calib_prev = (Double_t)VFTX_Bin2Ps[module][channel][ft-1];
    gain = calib - calib_prev; 
  }
  else {
    Double_t calib_next = (Double_t)VFTX_Bin2Ps[module][channel][ft+1];
    gain = calib_next - calib; 
  }
  double ft_ps = calib + gain*(double)rand ;
  double cc_ps = (double)(cc);
 // cout<<"ft_ps " <<ft_ps << " cc_ps " <<" cc_ps "<< cc_ps <<endl;
   
  return (5000.*cc_ps - ft_ps);  
}
