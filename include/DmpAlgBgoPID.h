#ifndef DmpAlgBgoPID_H
#define DmpAlgBgoPID_H

#include "DmpVAlg.h"
#include "DmpEvtBgoPID.h"
#include "DmpEvtBgoHits.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
class DmpEvtBgoPID;
class DmpEvtBgoHits;

 class DmpAlgBgoPID : public DmpVAlg{
/*
 *  DmpAlgBgoPID
 *
 */
public:
  DmpAlgBgoPID();
  ~DmpAlgBgoPID();

  //void Set(const std::string &type,const std::string &value);
  // if you need to set some options for your algorithm at run time. Overload Set()
  bool Initialize();
  bool ProcessThisEvent();    // only for algorithm
  bool Finalize();
  void Reset();
  double EnergyCor();

private:
  DmpEvtBgoPID *fBgoPID;
  DmpEvtBgoHits *fBgoHits;
  double LayerXC[14];            //layer transverse weighting center
  double LayerEX[14];            //layer sum of Energy*(transverse position)
  double rms[14];                //sum Energy*Power((x_li-xlc),2)
  double LgCenter;               //longitudinal energy center
  double LgCperGeV;              //LgCenter/Energy; mm/GeV
  double nHitsperGeV;

//  double RMS;//sum rms
  double Energy;
  double maximum;
  double BarE[14][22];
  double LayerE[14];//Cluster energy
  short  LayerCS[14];
  TH1D *TEnergy;

  TH2D *EnergyLg;
  TH1D *LgECDeE;
  TH2D *nHitsVsE;
  TH1D *nHitsDeE;

  TH2D *LgEC;//mm                //longitudinal energy weighting center
  TH1D *FValue[14];
  TH2D *RMSVsEFraction[14];
  TH2D *EnergyTrDis[14];

  TH1D *LayerEnergy[14];         //fill events that with total energy cut
  TProfile *EnergyLgCut; 

  //Corrected Energy
  TH1D *TEnergyCor;
  TH1D *fEnergy_LongitudinalCor;
};

#endif
