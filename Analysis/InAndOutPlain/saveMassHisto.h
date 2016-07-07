using namespace std;
#ifndef _ANA_INANDOUTPLAIN_SAVEMASSHISTO_H_
#define _ANA_INANDOUTPLAIN_SAVEMASSHISTO_H_

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

const int nCentBins = 6;
Int_t centBins[nCentBins+1] = {0, 10, 30, 50, 70, 90, 100};
Float_t EPm_resolution_v2_etagap[nCentBins] = {0.685732, 0.859684, 0.805492, 0.566930, 0.211378, 0.0307577};
Float_t EPp_resolution_v2_etagap[nCentBins] = {0.685895, 0.859866, 0.805762, 0.567147, 0.210694, 0.0329058};
Float_t EPm_resolution_v3_etagap[nCentBins] = {0.526727, 0.493562, 0.371182, 0.182194, 0.0346885, 99999.};
Float_t EPp_resolution_v3_etagap[nCentBins] = {0.528257, 0.495857, 0.373729, 0.183724, 0.0316327, 99999.};
Float_t SP_EPm_resolution_v2_etagap[nCentBins] = {2.20732, 4.50646, 5.36544, 3.85047, 1.16114, 0.121428};
Float_t SP_EPp_resolution_v2_etagap[nCentBins] = {2.22527, 4.54705, 5.41751, 3.87233, 1.15763, 0.123545};
Float_t SP_EPm_resolution_v3_etagap[nCentBins] = {1.1864, 1.44074, 1.42084, 0.878359, 0.164648, 99999.};
Float_t SP_EPp_resolution_v3_etagap[nCentBins] = {1.19751, 1.45886, 1.4413, 0.889368, 0.16584, 99999.};

TString tfend[4] = {"v2_inpl","v2_outpl","v3_inpl","v3_outpl"};
TString tpr[2] = {"Inclusive","Prompt"};

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

void DrawCmsTlatex(TString collision)
{
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} #bf{#it{Preliminary}}");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();

  TLatex* texCol = new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV",collision.Data()));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();
}

#endif
