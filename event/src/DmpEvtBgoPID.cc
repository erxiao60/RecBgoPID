/*   $Id: DmpEvtBgoPID.cc, 2014-10-19 09:43:16+08:00 DAMPE $
 *--------------------------------------------------------
 *  Author(s):
 *
 *--------------------------------------------------------
*/

#include "DmpEvtBgoPID.h"

ClassImp(DmpEvtBgoPID)

DmpEvtBgoPID::DmpEvtBgoPID()
{
}

DmpEvtBgoPID::~DmpEvtBgoPID()
{
}

void DmpEvtBgoPID::Reset()
{
fLayer.clear();
fFValue.clear();
fSqr_RMS.clear();
fLC_E.clear();
fLC_Size.clear();

fSize=0;
fEnergy=0.;
fLgCenter=0.;
fTrWidth=0.;
}
