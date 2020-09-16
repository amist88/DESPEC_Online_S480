#include "CorrelParameter.h"
#include <iostream>

using namespace std;

CorrelParameter::CorrelParameter()
: TGo4Parameter("AnaOnlineCorPar")
{
    GSetup_corr_FRS_Ge =false;
    GSetup_corr_FRS_Ge_long =false;
    GFRS_Ge_LongIso_incprmt = false;
    GSetup_corr_FRS_fat = false;
    GSetup_corr_FRS_Aida =false;
    GSetup_corr_Beta_Gamma=false;
    GSetup_corr_Beta_Gamma_Gamma=false;
    GSetup_corr_Beta_Gamma_Fatima=false;
   
    GbPlast_Egate_low = 0;
    GbPlast_Egate_high = 0;
    
    GFRS_AIDA_TLow = 0;
    GFRS_AIDA_THigh = 0;
    GAIDA_bPlas_TLow = 0;
    GAIDA_bPlas_THigh = 0;
    
    GLongGate = 0;
    GFRS_Ge_LongIso_THigh = 0;
    GFRS_Ge_LongIso_TLow = 0;
    GFRS_Ge_LongIso_HBin = 0;
    GFRS_Ge_LongIso_HLow = 0;
    GFRS_Ge_LongIso_HHigh = 0;
    GFRS_Ge_LongIso_TScale = 0;
    
    GAidaImpDecT_Low = 0;
    GAidaImpDecT_High = 0;
    GAidaImpDecT_HBin = 0;
    GAidaFB_dT = 0;
    GAidaFB_dE = 0;
    GAIDA_DecEFront = 0;
    GAIDA_DecEBack = 0; 
    
    GAida_Ge_WRdT_Low=0;
    GAida_Ge_WRdT_High=0;
    
    GAida_Fat_WRdT_High=0;
    GAida_Fat_WRdT_Low=0;
    
    GFRS_Fat_TLow = 0;
    GFRS_Fat_THigh = 0;
    GFRS_Ge_TLow = 0;
    GFRS_Ge_THigh = 0;
    
    
   
//     GbPlast_GAL_TLow = 0;
//     GbPlast_GAL_THigh = 0;
//     
    
//     GFRS_AIDA_DSSD1 = false;
//     GFRS_AIDA_DSSD2 = false;
//     GFRS_AIDA_DSSD3 = false;
    
}
//-------------------------------------------------------------------------------------------//
// constructor
// reads correlation parameters from the file Correlations_config.dat
CorrelParameter::CorrelParameter(const Text_t* name)
: TGo4Parameter(name)
{
  ifstream    file;
  
      file.open("Configuration_Files/DESPEC_General_Setup/Correlations_config.dat");
  if (file.fail()) {
        cout << "ERROR: CorrelParameter - Could not open file: Configuration_Files/DESPEC_General_Setup/Correlations_config.dat ! Correlations disabled!!\n"; 
         GSetup_corr_FRS_Ge =false;
         GSetup_corr_FRS_Ge_long=false;
         GFRS_Ge_LongIso_incprmt=false;
         GSetup_corr_FRS_fat=false;
         GSetup_corr_FRS_Aida =false;
         GSetup_corr_Beta_Gamma=false;
         GSetup_corr_Beta_Gamma_Gamma=false;
         GSetup_corr_Beta_Gamma_Fatima=false;
    
         GFRS_AIDA_TLow = -1000;
         GFRS_AIDA_THigh = 1000;
         GAIDA_bPlas_TLow = -1000;
         GAIDA_bPlas_THigh = 1000;
         
         GFRS_Ge_TLow = -1000;
         GFRS_Ge_THigh = 1000;
         GLongGate=0;
         GFRS_Ge_LongIso_TLow = 0;
         GFRS_Ge_LongIso_THigh = 20000;
         GFRS_Ge_LongIso_HBin = 700;
         GFRS_Ge_LongIso_HLow =0;
         GFRS_Ge_LongIso_HHigh = 700; 
         GFRS_Ge_LongIso_TScale = 1000;
         
         GAidaImpDecT_Low = 0;
         GAidaImpDecT_High = 5;
         GAidaImpDecT_HBin = 50;
         GAidaFB_dT = 3000;
         GAidaFB_dE = 200;
         GAIDA_DecEFront = 1500;
         GAIDA_DecEBack = 1500; 
         
         GAida_Ge_WRdT_Low=0;
         GAida_Ge_WRdT_High=20000;
         
         GAida_Fat_WRdT_High=0;
         GAida_Fat_WRdT_Low=20000;
        
         
         GFRS_Fat_TLow = -1000;
         GFRS_Fat_THigh = 1000;
         
//          GbPlast_GAL_TLow = -1000;
//          GbPlast_GAL_THigh = 1000;
        
         
//          GFRS_AIDA_DSSD1 = true;
//          GFRS_AIDA_DSSD2 = true;
//          GFRS_AIDA_DSSD3 = true;
//          
       
         GbPlast_Egate_low = 0;
         GbPlast_Egate_high = 5000;
  }
  
else {
cout << "CorrelParameter - reading from Configuration_Files/DESPEC_General_Setup/Correlations_config.dat";
        ///Detector system correlation pairs
       if(IsData(file)) file >> GSetup_corr_FRS_Aida;
       if(IsData(file)) file >> GSetup_corr_FRS_Ge;
       if(IsData(file)) file >> GSetup_corr_FRS_Ge_long >>GFRS_Ge_LongIso_incprmt;
       if(IsData(file)) file >> GSetup_corr_FRS_fat;
       if(IsData(file)) file >> GSetup_corr_Beta_Gamma >> GSetup_corr_Beta_Gamma_Gamma >> GSetup_corr_Beta_Gamma_Fatima;

    
       ///AIDA FRS WR Time Gate
       if(IsData(file)) file >> GFRS_AIDA_TLow >>GFRS_AIDA_THigh;
       ///FRS Galileo WR Time Gate
       if(IsData(file)) file >> GFRS_Ge_TLow >>GFRS_Ge_THigh;
        ///FRS Galileo PID Gate Long isomer
       if(IsData(file)) file >> GLongGate;
      
       ///FRS Galileo WR Time Gate Long isomer
       if(IsData(file)) file >> GFRS_Ge_LongIso_TLow >> GFRS_Ge_LongIso_THigh>> GFRS_Ge_LongIso_HBin>> GFRS_Ge_LongIso_HLow>> GFRS_Ge_LongIso_HHigh>> GFRS_Ge_LongIso_TScale;
      
         ///FRS Fatima WR Time Gate
       if(IsData(file)) file >> GFRS_Fat_TLow >>GFRS_Fat_THigh;
       
       ///AIDA Implant Decay time window
       if(IsData(file)) file >>GAidaImpDecT_Low >>GAidaImpDecT_High>>GAidaImpDecT_HBin;
       
       ///AIDA Front-back conditions
       if(IsData(file)) file >>GAidaFB_dT >>GAidaFB_dE;
       
       ///AIDA Max decay energy conditions
       if(IsData(file)) file >>GAIDA_DecEFront >>GAIDA_DecEBack;
       
       /// Plastic energy gates 
       if(IsData(file)) file >> GbPlast_Egate_low>> GbPlast_Egate_high;

       ///AIDA(Decays)-Ge WR dT 
       if(IsData(file)) file >>GAida_Ge_WRdT_Low >>GAida_Ge_WRdT_High;

       ///AIDA(Decays)-Fatima WR dT 
       if(IsData(file)) file >>GAida_Fat_WRdT_Low >>GAida_Fat_WRdT_High;
       
       ///AIDA bPlas WR Time Gate
       if(IsData(file)) file >> GAIDA_bPlas_TLow >>GAIDA_bPlas_THigh;
       
//         ///bPlas Galileo WR Time Gate
//        if(IsData(file)) file >> GbPlast_GAL_TLow >>GbPlast_GAL_THigh;
       
//        ///Use AIDA  DSSD Gate
//        if(IsData(file)) file >> GFRS_AIDA_DSSD1 >>GFRS_AIDA_DSSD2 >>GFRS_AIDA_DSSD3;
//        ///Fatima energy gates

       
   
       
       
       if (file.fail()) cout << " !!!There is a problem in the correlation file !!!";
    
  }
  file.close();
}
//--------------------------------------------------------------------------//
CorrelParameter::~CorrelParameter()
{}
//---------------------------------------------------------------------------//
Int_t CorrelParameter::PrintParameter(Text_t *buf, Int_t)
{
  cout << "\n Online Analysis Correlation Parameters: " << GetName() << endl;
     
  cout << "//////Detector Correlation pairs\n";
  
  cout<<"////FRS-Aida: " <<  GSetup_corr_FRS_Aida << endl;
  
  cout << "////FRS-Germanium Short isomers: " << GSetup_corr_FRS_Ge << ", FRS-Germanium Long isomers: " << GSetup_corr_FRS_Ge_long <<" (prompt gammas for LongCorr ("<<GFRS_Ge_LongIso_incprmt <<"))"<<  endl;
  
  cout << "////FRS-Fatima Short isomers: " <<GSetup_corr_FRS_fat <<endl;
  
  cout<<"////Beta-delayed Gammas: " <<  GSetup_corr_Beta_Gamma  <<"  Ge Gamma-Gamma: " << GSetup_corr_Beta_Gamma_Gamma << " BDG with Fatima: " << GSetup_corr_Beta_Gamma_Fatima << endl;

  cout << "//////AIDA FRS WR Time Gate\n";
  cout << "////Lower: " << GFRS_AIDA_TLow << "  \t Upper: = " <<  GFRS_AIDA_THigh <<  endl;
  
  cout << "//////FRS Galileo Time Gate\n";
  cout << "////Lower: " << GFRS_Ge_TLow << "  \t Upper:  " <<  GFRS_Ge_THigh <<  endl;
  
  cout << "//////FRS Galileo Time Gate Long isomer\n";
  cout << "////PID Gate Num: " << GLongGate << "  Lower: " << GFRS_Ge_LongIso_TLow << "  Upper: " <<  GFRS_Ge_LongIso_THigh << "\t Long Isomer Histo: (Bin, Low, High) ("<< GFRS_Ge_LongIso_HBin <<" " << GFRS_Ge_LongIso_HLow <<" " << GFRS_Ge_LongIso_HHigh << ") Long Isomer Scale Time "<< GFRS_Ge_LongIso_TScale  <<endl;
  
  cout << "//////FRS Fatima WR Time Gate\n";
  cout << "////Lower: " << GFRS_Fat_TLow << "  \t Upper: = " <<  GFRS_Fat_THigh <<  endl;
   
  cout << "//////Beta-delayed gamma conditions\n";
  cout << "//////AIDA Implant-Decay Time window\n";
  cout <<"////Lower: " << GAidaImpDecT_Low << "  \t Upper: " <<  GAidaImpDecT_High <<" \t dT Histo Binning " << GAidaImpDecT_HBin<<endl;
  
  cout << "//////AIDA(decays) Ge WR Time Gate \n";
  cout <<"////Lower: " << GAida_Ge_WRdT_Low << "  \t Upper: " <<  GAida_Ge_WRdT_High<<endl;
  
  cout << "//////AIDA(decays) Fatima WR Time Gate \n";
  cout <<"////Lower: " << GAida_Fat_WRdT_Low << "  \t Upper: " <<  GAida_Fat_WRdT_High<<endl;
  
  
  cout << "//////AIDA Front-Back Time Conditions \n";
  cout <<"////Aida Front-Back dT: " << GAidaFB_dT << "  \t Aida Front-Back dE: " <<  GAidaFB_dE<<endl;
  
  cout << "//////AIDA Max Decay Energy Conditions \n";
  cout <<"////Aida ('good') Decay Energy Max Front: " << GAIDA_DecEFront << "  \t Aida ('good') Decay Energy Max Back " <<  GAIDA_DecEBack<<endl;

    //   
  cout << "//////AIDA bPlas WR Time Gate\n";
  cout << "////Lower: " << GAIDA_bPlas_TLow << "  \t Upper: " <<  GAIDA_bPlas_THigh <<  endl;
  
  cout << "////bPlastic Energy Gate limits\n";
  cout << "////Lower: " << GbPlast_Egate_low << "  \t Upper: " << GbPlast_Egate_high <<  endl;
  
 
  
  
//   cout << "//////bPlas Galileo WR Time Gate\n";
//   cout << "////Lower: " << GbPlast_GAL_TLow << "  \t Upper: " <<  GbPlast_GAL_THigh <<  endl;

  
//   cout << "//////Use AIDA DSSD 2D position Gate\n";
//   cout << "////DSSD1: " << GFRS_AIDA_DSSD1 << "  \t DSSD2: " <<  GFRS_AIDA_DSSD2 << "  \t DSSD3: "<<GFRS_AIDA_DSSD3 <<endl;


  return 1;
}
//---------------------------------------------------------------------------//
Bool_t CorrelParameter::UpdateFrom(TGo4Parameter *pp)
{
  if(pp->InheritsFrom("CorrelParameter"))
  {
    CorrelParameter *from = (CorrelParameter *) pp;
      GSetup_corr_FRS_Aida = from->GSetup_corr_FRS_Aida;
      GSetup_corr_FRS_Ge = from->GSetup_corr_FRS_Ge;
      GSetup_corr_FRS_Ge_long = from->GSetup_corr_FRS_Ge_long;
      GSetup_corr_FRS_fat = from->GSetup_corr_FRS_fat;
      GFRS_Ge_LongIso_incprmt = from->GFRS_Ge_LongIso_incprmt;

      GSetup_corr_Beta_Gamma = from->GSetup_corr_Beta_Gamma;
      GSetup_corr_Beta_Gamma_Gamma = from->GSetup_corr_Beta_Gamma_Gamma;
      GSetup_corr_Beta_Gamma_Fatima = from->GSetup_corr_Beta_Gamma_Fatima;
      
      GFRS_AIDA_TLow = from->GFRS_AIDA_TLow;
      GFRS_AIDA_THigh = from->GFRS_AIDA_THigh;
     
      
      GFRS_Ge_TLow = from->GFRS_Ge_TLow; 
      GFRS_Ge_THigh = from->GFRS_Ge_THigh; 
      GLongGate = from-> GLongGate;
      GFRS_Ge_LongIso_TLow = from-> GFRS_Ge_LongIso_TLow;
      GFRS_Ge_LongIso_THigh = from-> GFRS_Ge_LongIso_THigh;
      GFRS_Ge_LongIso_HBin = from-> GFRS_Ge_LongIso_HBin;
      GFRS_Ge_LongIso_HLow = from->GFRS_Ge_LongIso_HLow;
      GFRS_Ge_LongIso_HHigh = from->GFRS_Ge_LongIso_HHigh;
      GFRS_Ge_LongIso_TScale = from-> GFRS_Ge_LongIso_TScale;
      
      GAidaImpDecT_Low  = from-> GAidaImpDecT_Low;
      GAidaImpDecT_High = from-> GAidaImpDecT_High;
      GAidaImpDecT_HBin  = from->GAidaImpDecT_HBin;
      
      GAidaFB_dT = from->GAidaFB_dT;
      GAidaFB_dE = from->GAidaFB_dE;
      
      GAIDA_DecEFront = from->GAIDA_DecEFront;
      GAIDA_DecEBack = from->GAIDA_DecEBack;
  
//       GbPlast_GAL_TLow = from->GbPlast_GAL_TLow; 
//       GbPlast_GAL_THigh = from->GbPlast_GAL_THigh; 
      
      GFRS_Fat_TLow = from->GFRS_Fat_TLow; 
      GFRS_Fat_THigh = from->GFRS_Fat_THigh; 
     
      
      GAIDA_bPlas_TLow = from->GAIDA_bPlas_TLow;
      GAIDA_bPlas_THigh = from->GAIDA_bPlas_THigh;
      
//       GFRS_AIDA_DSSD1 = from->GFRS_AIDA_DSSD1;
//       GFRS_AIDA_DSSD2 = from->GFRS_AIDA_DSSD2;
//       GFRS_AIDA_DSSD3 = from->GFRS_AIDA_DSSD3;
//       
    
      GbPlast_Egate_low  = from->GbPlast_Egate_low;
      GbPlast_Egate_high = from->GbPlast_Egate_high;

       cout << "CorrelParameter - Parameter : " << GetName() << " UPDATED\n";
  }
  else {
      cout << "ERROR: CorrelParameter - Wrong parameter object: " << pp->ClassName() << endl;}
        return kTRUE;
  }
  //---------------------------------------------------------------------------//
  int CorrelParameter::IsData(ifstream &f) {
  char dum;
  char dumstr[300];
  int retval = 0;

  /* 'operator >>' does not read End-of-line, therefore check if read 
      character is not EOL (10) */
  do {
    dum = f.get();
    if (dum == '#')    // comment line => read whole line and throw it away
      f.getline(dumstr,300);
  }
  while ((dum == '#') || ((int)dum == 10) || (dum == ' ')); 

  f.unget();   // go one character back
  retval = 1;
  return retval;
}
ClassImp(CorrelParameter)
