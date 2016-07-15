using namespace std;
#ifndef _ANA_BFEEDDOWN_SAVEMASSHISTO_H_
#define _ANA_BFEEDDOWN_SAVEMASSHISTO_H_

#include "../include/uti.h"
#include "../include/parameters.h"
#include "../include/setBranches.h"

int findptbin(Int_t j)
{
  Int_t ipt = -1;
  for(int i=0;i<nPtBins;i++)
    {
      if(Dpt[j]>=ptBins[i] && Dpt[j]<ptBins[i+1])
        {
          ipt = i;
          break;
        }
    }
  return ipt;
}

bool passtriggersel()
{
  return true;
}

bool passeventsel(Float_t centmin, Float_t centmax)
{
  Bool_t returnval = false;
  if(pcollisionEventSelection &&
     pprimaryVertexFilter &&
     phfCoincFilter3 &&
     pclusterCompatibilityFilter &&
     TMath::Abs(PVz) < 15.0)
    {
      if(hiBin >= centmin*2 && hiBin < centmax*2)
        {
          returnval = true;
        }
    }
  return returnval;
}

bool passcutsel(Int_t iptbin, Int_t j)
{
  Bool_t returnval = false;
  if(Dtrk1highPurity[j] && Dtrk1highPurity[j] &&
     TMath::Abs(Dtrk1Eta[j]) < 1.2 && TMath::Abs(Dtrk2Eta[j]) < 1.2 &&
     Dtrk1Pt[j] > 0.7 && Dtrk2Pt[j] > 0.7 &&
     TMath::Abs(Dtrk1PtErr[j]/Dtrk1Pt[j]) < 0.1 && TMath::Abs(Dtrk2PtErr[j]/Dtrk2Pt[j]) < 0.1 &&
     (Dtrk1PixelHit[j]+Dtrk1StripHit[j]) >= 10.5 && (Dtrk2PixelHit[j]+Dtrk2StripHit[j]) >= 10 &&
     (Dtrk1Chi2ndf[j]/(Dtrk1nStripLayer[j]+Dtrk1nPixelLayer[j])) < 0.15 && (Dtrk2Chi2ndf[j]/(Dtrk2nStripLayer[j]+Dtrk2nPixelLayer[j])) < 0.15 )
    {
      if(TMath::Abs(Dy[j]) < 1.)
        {
          if(Dalpha[j] < 0.12 && 
             (DsvpvDistance[j]/DsvpvDisErr[j]) > ffls3dcut[iptbin] && 
             Dchi2cl[j] > vprobcut[iptbin] )
            {
              returnval = true;
            }
        }
    }
  return returnval;
}

bool passsigreg(Int_t j)
{
  if(TMath::Abs(Dmass[j]-1.8649)<0.025) return true;
  else return false;
}

bool passsidbnd(Int_t j)
{
  if(TMath::Abs(Dmass[j]-1.8649)>0.075 && TMath::Abs(Dmass[j]-1.8649)<0.1) return true;
  else return false;
}

bool passprompt(Int_t j)
{
  if(DgenBAncestorpt[j]<=0) return true;
  else return false;
}

int findcentbin(Float_t centmin, Float_t centmax)
{
  for(int i=0;i<nCentBins;i++)
    {
      if(centmin>=centBins[i] && centmax<=centBins[i+1]) return i;
    }
  cout<<"Error: Cannot find valid centrality bin"<<endl;
  return -1;
}

void caldeltaphi(Int_t j, Float_t* result)
{
  result[0] = -1;
  result[1] = -1;
  Float_t v2delphi,v3delphi;
  if(Dy[j]>0)
    {
      v2delphi = TMath::Abs(Dphi[j]-hiEvtPlanes[0]);
      v3delphi = TMath::Abs(Dphi[j]-hiEvtPlanes[17]);
    }
  else
    {
      v2delphi = TMath::Abs(Dphi[j]-hiEvtPlanes[1]);
      v3delphi = TMath::Abs(Dphi[j]-hiEvtPlanes[18]);
    }

  if(v2delphi<(PI/2.)) result[0] = v2delphi;
  else if(v2delphi>=(PI/2.) && v2delphi<PI) result[0] = PI-v2delphi;
  else if(v2delphi>=PI && v2delphi<(3.*PI/2.)) result[0] = v2delphi-PI;
  else
    {
      cout<<"Error: Invalid v2 delta phi - delta phi = "<<(v2delphi/PI)<<" Pi"<<endl;
    }

  if(v3delphi<(PI/3.)) result[1] = v3delphi;
  else if(v3delphi>=(PI/3.) && v3delphi<(2.*PI/3.)) result[1] = (2./3.)*PI-v3delphi;
  else if(v3delphi>=(2.*PI/3.) && v3delphi<PI) result[1] = v3delphi-(2./3.)*PI;
  else if(v3delphi>=PI && v3delphi<(4.*PI/3.)) result[1] = (4./3.)*PI-v3delphi;
  else
    {
      cout<<"Error: Invalid v3 delta phi - delta phi = "<<(v3delphi/PI)<<" Pi"<<endl;
    }
}

#endif
