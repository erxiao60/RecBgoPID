/*   $Id: DmpEvtBgoPID.h, 2014-10-19 09:43:16+08:00 DAMPE $
 *--------------------------------------------------------
 *  Author(s):
 *
 *--------------------------------------------------------
*/

#ifndef DmpEvtBgoPID_H
#define DmpEvtBgoPID_H

#include "TObject.h"

class DmpEvtBgoPID : public TObject{
/*
 *  DmpEvtBgoPID
 *
 */
public:
  DmpEvtBgoPID();
  ~DmpEvtBgoPID();

  void Reset();

public:
  std::vector<short> fLayer;
  std::vector<double> fFValue;
  std::vector<double> fSqr_RMS;
  std::vector<double> fLC_E;//Layer Cluster Energy
  std::vector<short> fLC_Size;//Layer Cluster Size

  short  fSize;
  double fEnergy;
  double fLgCenter;//longitudinal energy weigting center
  double fTrWidth;//fTransverseWidth; energy weigting
private:
  /*
   *  Define your data members at here
   *
   */
  
  ClassDef(DmpEvtBgoPID,1)

};

#endif
