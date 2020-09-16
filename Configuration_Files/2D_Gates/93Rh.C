{
//========= Macro generated from object: cID_Z_AoQ0/
//========= by ROOT version6.16/00
   
   TCutG *cutg = new TCutG("cID_Z_AoQ0",8);
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   cutg->SetFillColor(ci);
   cutg->SetFillStyle(3444);

   ci = TColor::GetColor("#ff0000");
   cutg->SetLineColor(ci);
   cutg->SetPoint(0,2.14893,48.2871);
   cutg->SetPoint(1,2.15693,47.7168);
   cutg->SetPoint(2,2.16793,47.5113);
   cutg->SetPoint(3,2.17893,47.8132);
   cutg->SetPoint(4,2.17713,48.1739);
   cutg->SetPoint(5,2.17073,48.5639);
   cutg->SetPoint(6,2.15693,48.6226);
   cutg->SetPoint(7,2.14893,48.2871);
   cutg->Draw("");
}
