#include "DmpAlgBgoPID.h"
#include "DmpDataBuffer.h"
#include "DmpBgoBase.h"
#include "TMath.h"
#include "TVector3.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
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
memset(LayerXC,0,sizeof(LayerXC));
memset(LayerEX,0,sizeof(LayerEX));
memset(rms,0,sizeof(rms));
memset(BarE,0,sizeof(BarE));
memset(LayerE,0,sizeof(LayerE));
memset(LayerCS,0,sizeof(LayerCS));
LgCenter=0.;
LgCperGeV=0.;
nHitsperGeV=0.;
Energy=0.;
maximum=0.;

}
//-------------------------------------------------------------------
bool   DmpAlgBgoPID::Initialize(){
  
  // read input data
  fBgoHits= new DmpEvtBgoHits();
  fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/Cal/Hits"));
  if(!fBgoHits){
    gDataBuffer->LinkRootFile("Event/Cal/Hits",fBgoHits);
    fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/Cal/Hits"));
  }
  if(!fBgoHits){
    gDataBuffer->LinkRootFile("Event/MCTruth/BgoFDigit",fBgoHits);
    fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/MCTruth/BgoFDigit"));
  }
  fBgoPID = new DmpEvtBgoPID();
  gDataBuffer->RegisterObject("Event/Rec/PID",fBgoPID,"DmpEvtBgoPID");
  
  //Define histograms
TEnergy=new TH1D("TEnergy","TEnergy;Energy(MeV);Counts",300000,0,300000);
TEnergyCor=new TH1D("EnergyCor","EnergyCor;Energy(MeV);Counts",300000,0,300000);

EnergyLg=new TH2D("Energy_Distribution","EnergyDepositVsLayers;Layers;Energy(MeV)",16,-1,15,300000,0,300000);
LgECDeE=new TH1D("LgCenter/E","LC/E;LC/E;Counts",100,0,10);
nHitsVsE=new TH2D("nHits","nHits;Energy(MeV);nHits(>3*PedestalSigma)",300000,0,300000,200,0,200);
nHitsDeE=new TH1D("nHitsPerGeV","nHitsPeGeV;nHits;counts",200,0,50);


LgEC=new TH2D("LongituralWeigtingCenter","LgWCenter;Energy(MeV);LgEC(mm)",300000,0,300000,100,0,100);//mm
  for (int i=0;i<14;i++){   
  char ff[80];
  char rr[80];
  char ee[80];
  char ll[80];
  sprintf(ff,"FValue_Layer%02d",i);
  sprintf(rr,"RMS_Versus_EnergyFraction_L%02d",i);
  sprintf(ee,"EnergyTransverseDistrabution_Layer%02d",i);
  sprintf(ll,"Energy_Layer%02d",i);
  FValue[i]=new TH1D(ff,ff,100,0,500);
  RMSVsEFraction[i]=new TH2D(rr,rr,200,0,2000,100,0,0.8);
  EnergyTrDis[i]=new TH2D(ee,ee,24,-1,23,300000,0,300000);
  LayerEnergy[i]=new TH1D(ll,ll,300000,0,300000);
  }
  EnergyLgCut=new TProfile("LE_DisProfile","LE_DisProfile;Layers;Energy(MeV)",16,-1,15);
  
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
  BarE[l][b]=fBgoHits->fEnergy[n];
  EnergyTrDis[l]->Fill(b,BarE[l][b]);

  if(fBgoHits->fEnergy[n]>300000){std::cout<<"Error: "<<" L:"<<l<<" B:"<<b<<" Energy:"<<fBgoHits->fEnergy[n];continue;}
  if(l%2==0)
 // LayerEX[l]+= fBgoHits->fEnergy[n]*fBgoHits->fPosition[n].y();  
  LayerEX[l]+= fBgoHits->fEnergy[n]*b; 
  else
 // LayerEX[l]+= fBgoHits->fEnergy[n]*fBgoHits->fPosition[n].x();
  LayerEX[l]+= fBgoHits->fEnergy[n]*b;

  LayerE[l]+= fBgoHits->fEnergy[n];
  LayerCS[l]++;
  Energy+= fBgoHits->fEnergy[n];
  }
  for(short nl=0;nl<14;nl++){ 
    LayerXC[nl]=LayerEX[nl]/LayerE[nl];
  }
  //then,Loop for fSqr_RMS
  for(short n=0;n<nHits;++n){ 
  gid=fBgoHits->fGlobalBarID[n];
  l=DmpBgoBase::GetLayerID(gid);
  b=DmpBgoBase::GetBarID(gid);
  if(l%2==0)
  rms[l]+=fBgoHits->fEnergy[n]*TMath::Power((fBgoHits->fPosition[n].y()-LayerXC[l]),2);
  else
  rms[l]+=fBgoHits->fEnergy[n]*TMath::Power((fBgoHits->fPosition[n].x()-LayerXC[l]),2);
  }
  //Fill event class
  for(short nl=0;nl<14;nl++){ 
    fBgoPID->fLayer.push_back(nl);
    fBgoPID->fSqr_RMS.push_back(rms[nl]/LayerE[nl]);
    fBgoPID->fFValue.push_back(rms[nl]/Energy);
    fBgoPID->fLayerE.push_back(LayerE[nl]);//MeV
    fBgoPID->fLayerCS.push_back(LayerCS[nl]);
    fBgoPID->fLgCenter+=LayerE[nl]*nl/Energy;
    fBgoPID->fTrWidth+=LayerE[nl]*LayerCS[nl]/Energy;
} 
  
  fBgoPID->fEnergy=Energy;
  if(Energy>maximum){maximum=Energy;}

  //Fill Histograms
  //TH1D *TEnergy=new TH1D("TEnergy","TEnergy;Energy(MeV);Counts",200,0,500000);
  //TH2D *EnergyLg=new TH2D("Energy_Distribution","EnergyDepositVsLayers;Layers;Energy(MeV)",16-1,15,200,0,20000);
  //TH2D *nHitsVsE=new TH2D("nHits","nHits;Energy(MeV);nHits(>3*PedestalSigma)",200,0,20000,200,0,200);
  //TH2D *LgEC=new TH2D("LongituralWeigtingCenter","LgWCenter;Energy(MeV);LgEC(mm)",200,0,20000,800,0,800);//mm
  //TH1D *FValue[14];
  //TH2D *RMSVsEFraction[14];
  //TH2D *EnergyTrDis[14];

  TEnergy->Fill(Energy);
  TEnergyCor->Fill(Energy+EnergyCor());
  nHitsVsE->Fill(Energy,fBgoPID->fSize);
  nHitsDeE->Fill(fBgoPID->fSize/Energy*1000);
  LgEC->Fill(Energy,fBgoPID->fLgCenter);
  LgECDeE->Fill(fBgoPID->fLgCenter/Energy*1000);

for(int i=0;i<14;i++){
EnergyLg->Fill(i,LayerE[i]);
EnergyLgCut->Fill(i,LayerE[i]);
FValue[i]->Fill(fBgoPID->fFValue[i]);
RMSVsEFraction[i]->Fill(fBgoPID->fSqr_RMS[i],LayerE[i]/Energy);
LayerEnergy[i]->Fill(LayerE[i]);

}
// EnergyTrDis[i]->Fill();
// }
//  if(Energy>200000){CorFitting();}
  return true;
} 
//-------------------------------------------------------------------
double DmpAlgBgoPID::EnergyCor(){
  double CorL_E[14];
  memset(CorL_E,0,sizeof(CorL_E));
  double CorE=0.;
  double Tratio=(2.5/1.12)/(0.15/26.6+0.1/22.9);//radiation length ratio of two crystal/their gap
  double Lratio=(2.5/1.12)/(0.30/26.6+0.1/22.9);//radiation length ratio of two crystal/their gap
  for(int layer=0;layer<14;layer++){
    //transverse correction
    for(int bar=0;bar<21;bar++){
      double GapL_E=TMath::Sqrt(BarE[layer][bar]*BarE[layer][bar+1])/Tratio;
   //   double GapL_E=(BarE[layer][bar]+BarE[layer][bar+1])/2;
   //   double GapL_E=0.;
   //   if(BarE[layer][bar]>BarE[layer][bar+1]){GapL_E=BarE[layer][bar];}
   //   else{GapL_E=BarE[layer][bar+1];}
      CorE+=GapL_E;
      CorL_E[layer]+=GapL_E;
    }
  }
    //longitudinal correction
  for(int layer=0;layer<13;layer++){
   double GapE=TMath::Sqrt((LayerE[layer]+CorL_E[layer])*(LayerE[layer]+CorL_E[layer]))/Lratio;
   // double GapE=((LayerE[layer]+CorL_E[layer])+(LayerE[layer]+CorL_E[layer]))/2;
   // double GapE=LayerE[layer]+CorL_E[layer];
    CorE+=GapE;
  }
  return CorE;
}
//-------------------------------------------------------------------
bool DmpAlgBgoPID::Finalize(){
  
  std::string histFileName ="./PID/Histograms/"+gRootIOSvc->GetInputStem()+"_PID_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
 //gStyle->SetOptFit(1111);
  TF1 *mygaus=new TF1("mygaus","gaus",80,300000);
  double max=maximum;
  int nBin=1;
  std::cout<<"Max: "<<max<<std::endl;
  max=((int)((max*1.5)/1000+0.5))*1000;
  if (max<1000){ max=1000;nBin=2;}
  else if(max>=1000&&max<5000){nBin=5;}
  else if(max>=5000&&max<10000){nBin=10;}
  else if(max>=10000&&max<50000){nBin=50;}
  else if(max>=50000&&max<100000){nBin=100;}
  else{nBin=200;}
  TEnergy->RebinX(nBin);
  TEnergy->SetAxisRange(0,max,"X");
  double mean=TEnergy->GetMean();
  double par[3];
  
  mygaus->SetParameter(1,mean);
  mygaus->SetRange(0.5*mean,1.5*mean);
  TEnergy->Fit(mygaus,"QR0");
  mygaus->SetRange(par[1]-2.5*par[2],par[1]+3*par[2]);
  TEnergy->Fit(mygaus,"QR");
  TEnergy->Write();
  delete TEnergy;
  
  TEnergyCor->RebinX(nBin);
  TEnergyCor->SetAxisRange(0,max,"X");
  TEnergyCor->Write();
  delete TEnergyCor;
  
  //EnergyLg->RebinY(nBin);
  EnergyLg->SetAxisRange(0,0.5*max,"Y");
  EnergyLg->Write();
  delete EnergyLg;

  EnergyLgCut->Write();
  delete EnergyLgCut;


  nHitsVsE->RebinX(nBin);
  nHitsVsE->SetAxisRange(0,max,"X");
  nHitsVsE->Write();
  delete nHitsVsE;
  
  nHitsDeE->Write();
  delete nHitsDeE;



  LgEC->RebinX(nBin);
  LgEC->SetAxisRange(0,max,"X");
  LgEC->Write();
  delete LgEC;

  LgECDeE->Write();
  delete LgECDeE;
   
  for(int i=0;i<14;i++){
    FValue[i]->Write();
    delete FValue[i];
  }
  for(int i=0;i<14;i++){
    RMSVsEFraction[i]->Write();
    delete RMSVsEFraction[i];
  }
  for(int i=0;i<14;i++){
    LayerEnergy[i]->RebinX(nBin);
    LayerEnergy[i]->SetAxisRange(0,0.5*max,"X");
    LayerEnergy[i]->Write();
    delete LayerEnergy[i];
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

