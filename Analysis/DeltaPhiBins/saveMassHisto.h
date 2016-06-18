#ifndef _ANA_DELTAPHIBINS_SAVEMASSHISTO_H_
#define _ANA_DELTAPHIBINS_SAVEMASSHISTO_H_

#include "uti.h"
#include "setBranches.h"
const float PI = 3.14159265359;

Double_t minhisto=1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

const int nPtBins = 8;
Float_t ptBins[nPtBins+1] = {2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 25.0, 40.0};
Float_t ffls3dcut[nPtBins] = {5.86,  5.86,  4.86,  4.54,  4.42,  4.06,  3.50,  3.00};
Float_t vprobcut[nPtBins] =  {0.224, 0.224, 0.170, 0.125, 0.091, 0.069, 0.054, 0.050};
const int nPhiBins = 5;

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

void findhistno(Int_t j, Bool_t isData, Int_t* histno)
{
  histno[0] = -1;
  histno[1] = -1;
  for(int i=0;i<nPtBins;i++)
    {
      if(Dpt[j]>ptBins[i] && Dpt[j]<ptBins[i+1])
        {
          histno[0] = i;
          break;
        }
    }
  if(isData)
    {
      Float_t dphistep = (PI/2.)/nPhiBins;
      Float_t deltaphi = -1;
      if(Dy[j]>0.)
        {
          if(TMath::Abs(Dphi[j]-hiEvtPlanes[0])<(PI/2.)) deltaphi = TMath::Abs(Dphi[j]-hiEvtPlanes[0]);
          else if(TMath::Abs(Dphi[j]-hiEvtPlanes[0])>=(PI/2.) && TMath::Abs(Dphi[j]-hiEvtPlanes[0])<PI) deltaphi = PI-TMath::Abs(Dphi[j]-hiEvtPlanes[0]);
          else if(TMath::Abs(Dphi[j]-hiEvtPlanes[0])>=PI && TMath::Abs(Dphi[j]-hiEvtPlanes[0])<(3.*PI/2.)) deltaphi = TMath::Abs(Dphi[j]-hiEvtPlanes[0])-PI;
          else
            {              
              cout<<"Error: Invalid delta phi - [positive] delta phi = "<<(TMath::Abs(Dphi[j]-hiEvtPlanes[0])/PI)<<" Pi"<<endl;
            }
        }
      else
        {
          if(TMath::Abs(Dphi[j]-hiEvtPlanes[1])<(PI/2.)) deltaphi = TMath::Abs(Dphi[j]-hiEvtPlanes[1]);
          else if(TMath::Abs(Dphi[j]-hiEvtPlanes[1])>=(PI/2.) && TMath::Abs(Dphi[j]-hiEvtPlanes[1])<PI) deltaphi = PI-TMath::Abs(Dphi[j]-hiEvtPlanes[1]);
          else if(TMath::Abs(Dphi[j]-hiEvtPlanes[1])>=PI && TMath::Abs(Dphi[j]-hiEvtPlanes[1])<(3.*PI/2.)) deltaphi = TMath::Abs(Dphi[j]-hiEvtPlanes[1])-PI;
          else
            {
              cout<<"Error: Invalid delta phi - [negative] delta phi = "<<(TMath::Abs(Dphi[j]-hiEvtPlanes[1])/PI)<<" Pi"<<endl;
            }
        }
      for(int i=0;i<nPhiBins;i++)
        {
          if(deltaphi>(dphistep*i) && deltaphi<(dphistep*(i+1)))
            {
              histno[1] = i;
              break;
            }
        }
    }
}

#endif
