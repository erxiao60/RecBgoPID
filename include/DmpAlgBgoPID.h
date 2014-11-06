#ifndef DmpAlgBgoPID_H
#define DmpAlgBgoPID_H

#include "DmpVAlg.h"
#include "DmpEvtBgoPID.h"
#include "DmpEvtBgoHits.h"
#include "TH1D.h"
#include "TH2D.h"

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

private:
  DmpEvtBgoPID *fBgoPID;
  DmpEvtBgoHits *fBgoHits;
  double x_lc[14];
  double Ex_lc[14];
  double rms[14];//sum Energy*Power((x_li-xlc),2)

//  double RMS;//sum rms
  double Energy;
  double LC_E[14];//Cluster energy
  short  LC_Size[14];
  TH1D *TEnergy;

  TH2D *EnergyLgDis;
  TH2D *nHitsVsE;
  TH2D *Lg_Center;//mm
  TH1D *FValue[14];
  TH2D *RMSVsEFraction[14];
  TH2D *EnergyTrDis[14];
};

#endif
