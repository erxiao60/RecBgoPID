#include "DmpAlgBgoPID.h"
#include "DmpDataBuffer.h"
#include "DmpBgoBase.h"
#include "TMath.h"
#include "TVector3.h"
#include "TFile.h"
#include "TF1.h"
//-------------------------------------------------------------------
DmpAlgBgoPID::DmpAlgBgoPID()
 :DmpVAlg("BgoPID"),
  fBgoPID(0),
  fBgoHits(0)
{ 
}

//-------------------------------------------------------------------
DmpAlgBgoPID::~DmpAlgBgoPID(){
} 

//-------------------------------------------------------------------
void DmpAlgBgoPID::Reset(){
memset(x_lc,0,sizeof(x_lc));
memset(Ex_lc,0,sizeof(Ex_lc));
memset(rms,0,sizeof(rms));
memset(LC_E,0,sizeof(LC_E));
memset(LC_Size,0,sizeof(LC_Size));
Energy=0.;

}
//-------------------------------------------------------------------
bool   DmpAlgBgoPID::Initialize(){
  
  // read input data
  fBgoHits= new DmpEvtBgoHits();
  gDataBuffer->LinkRootFile("Event/Cal/Hits",fBgoHits);
  fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/Cal/Hits"));
  if(!fBgoHits){
    gDataBuffer->LinkRootFile("Event/MCTruth/BgoFDigit",fBgoHits);
    fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/MCTruth/BgoFDigit"));
  } 

  fBgoPID = new DmpEvtBgoPID();
  gDataBuffer->RegisterObject("Event/Rec/PID",fBgoPID,"DmpEvtBgoPID");
  
  //Define histograms
TEnergy=new TH1D("TEnergy","TEnergy;Energy(MeV);Counts",300000,0,300000);

EnergyLgDis=new TH2D("Energy_Distribution","EnergyDepositVsLayers;Layers;Energy(MeV)",16,-1,15,300000,0,300000);
nHitsVsE=new TH2D("nHits","nHits;Energy(MeV);nHits(>3*PedestalSigma)",300000,0,300000,200,0,200);
Lg_Center=new TH2D("LongituralWeigtingCenter","LgWCenter;Energy(MeV);Lg_Center(mm)",300000,0,300000,100,0,100);//mm
  for (int i=0;i<14;i++){   
  char ff[80];
  char rr[80];
  char ee[80];
  sprintf(ff,"FValue_Layer%02d",i);
  sprintf(rr,"RMS_Versus_EnergyFraction_L%02d",i);
  sprintf(ee,"EnergyTransverseDistrabution_Layer%02d",i);
  FValue[i]=new TH1D(ff,ff,100,0,500);
  RMSVsEFraction[i]=new TH2D(rr,rr,200,0,2000,100,0,0.8);
  EnergyTrDis[i]=new TH2D(ee,ee,24,-1,23,300000,0,300000);
  }
  
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoPID::ProcessThisEvent(){
  fBgoPID->Reset();
  Reset();

  //Loop Hits
  short nHits=fBgoHits->fGlobalBarID.size();
  short gid=0,l=0,b=0;
  fBgoPID->fSize=nHits;
  //firt,Loop for x layer centers
  for(short n=0;n<nHits;++n){
  gid=fBgoHits->fGlobalBarID[n];
  l=DmpBgoBase::GetLayerID(gid);
  b=DmpBgoBase::GetBarID(gid);
  EnergyTrDis[l]->Fill(b,fBgoHits->fEnergy[n]);
  if(l%2==0)
  Ex_lc[l]+= fBgoHits->fEnergy[n]*fBgoHits->fPosition[n].y();  
  else
  Ex_lc[l]+= fBgoHits->fEnergy[n]*fBgoHits->fPosition[n].x();

  LC_E[l]+= fBgoHits->fEnergy[n];
  LC_Size[l]++;
  Energy+= fBgoHits->fEnergy[n];
  }
  for(short nl=0;nl<14;nl++){ 
    x_lc[nl]=Ex_lc[nl]/LC_E[nl];
  }
  //then,Loop for fSqr_RMS
  for(short n=0;n<nHits;++n){ 
  gid=fBgoHits->fGlobalBarID[n];
  l=DmpBgoBase::GetLayerID(gid);
  b=DmpBgoBase::GetBarID(gid);
  if(l%2==0)
  rms[l]+=fBgoHits->fEnergy[n]*TMath::Power((fBgoHits->fPosition[n].y()-x_lc[l]),2);
  else
  rms[l]+=fBgoHits->fEnergy[n]*TMath::Power((fBgoHits->fPosition[n].x()-x_lc[l]),2);
  }
  //Fill event class
  for(short nl=0;nl<14;nl++){ 
    fBgoPID->fLayer.push_back(nl);
    fBgoPID->fSqr_RMS.push_back(rms[nl]/LC_E[nl]);
    fBgoPID->fFValue.push_back(rms[nl]/Energy);
    fBgoPID->fLC_E.push_back(LC_E[nl]);//MeV
    fBgoPID->fLC_Size.push_back(LC_Size[nl]);
    fBgoPID->fLgCenter+=LC_E[nl]*nl/Energy;
    fBgoPID->fTrWidth+=LC_E[nl]*LC_Size[nl]/Energy;
  } 
  fBgoPID->fEnergy=Energy;

  //Fill Histograms
  //TH1D *TEnergy=new TH1D("TEnergy","TEnergy;Energy(MeV);Counts",200,0,500000);
  //TH2D *EnergyLgDis=new TH2D("Energy_Distribution","EnergyDepositVsLayers;Layers;Energy(MeV)",16-1,15,200,0,20000);
  //TH2D *nHitsVsE=new TH2D("nHits","nHits;Energy(MeV);nHits(>3*PedestalSigma)",200,0,20000,200,0,200);
  //TH2D *Lg_Center=new TH2D("LongituralWeigtingCenter","LgWCenter;Energy(MeV);Lg_Center(mm)",200,0,20000,800,0,800);//mm
  //TH1D *FValue[14];
  //TH2D *RMSVsEFraction[14];
  //TH2D *EnergyTrDis[14];

  TEnergy->Fill(Energy);
  nHitsVsE->Fill(Energy,fBgoPID->fSize);
  Lg_Center->Fill(Energy,fBgoPID->fLgCenter);

for(int i=0;i<14;i++){
EnergyLgDis->Fill(i,LC_E[i]);
FValue[i]->Fill(fBgoPID->fFValue[i]);
RMSVsEFraction[i]->Fill(fBgoPID->fSqr_RMS[i],LC_E[i]/Energy);
}
// EnergyTrDis[i]->Fill();
// }
  return true;
} 

//-------------------------------------------------------------------
bool DmpAlgBgoPID::Finalize(){
  
  std::string histFileName ="./PID/Histograms/"+gRootIOSvc->GetInputStem()+"_PID_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  double max=TEnergy->GetMean();
  int nBin=1;
  std::cout<<"mean: "<<max<<std::endl;
  max=((int)((max*5)/1000))*1000;
  if (max<1000){ max=1000;nBin=10;}
  else if(max>=1000&&max<10000){nBin=10;}
  else if(max>=10000&&max<100000){nBin=100;}
  else{nBin=300;}
  TEnergy->RebinX(nBin);
  TEnergy->SetAxisRange(0,max,"X");
  TEnergy->Write();
  delete TEnergy;
  
  
  EnergyLgDis->RebinY(nBin);
  EnergyLgDis->SetAxisRange(0,0.5*max,"Y");
  EnergyLgDis->Write();
  delete EnergyLgDis;

  nHitsVsE->RebinX(nBin);
  nHitsVsE->SetAxisRange(0,max,"X");
  nHitsVsE->Write();
  delete nHitsVsE;

  Lg_Center->RebinX(nBin);
  Lg_Center->SetAxisRange(0,max,"X");
  Lg_Center->Write();
  delete Lg_Center;
   
  for(int i=0;i<14;i++){
    FValue[i]->Write();
    delete FValue[i];
  }
  for(int i=0;i<14;i++){
    RMSVsEFraction[i]->Write();
    delete RMSVsEFraction[i];
  }
  for(int i=0;i<14;i++){
    EnergyTrDis[i]->RebinX(nBin);
    EnergyTrDis[i]->SetAxisRange(0,0.5*max,"Y");
    EnergyTrDis[i]->Write();
    delete EnergyTrDis[i];
  }

  histFile->Close(); 
  return true;
}

