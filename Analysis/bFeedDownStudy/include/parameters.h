using namespace std;
#ifndef _ANA_BFEEDDOWN_PARAMETERS_H_
#define _ANA_BFEEDDOWN_PARAMETERS_H_

#include "uti.h"
#include "setBranches.h"
const float PI = 3.14159265359;

const int nPtBins = 9;
Double_t ptBins[nPtBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 25.0, 40.0};
Double_t ffls3dcut[nPtBins] = {6.00,  5.86,  5.86,  4.86,  4.54,  4.42,  4.06,  3.50,  3.00};
Double_t vprobcut[nPtBins] =  {0.250, 0.224, 0.224, 0.170, 0.125, 0.091, 0.069, 0.054, 0.050};
const int nPhiBins = 5;

const int nCentBins = 6;
Int_t centBins[nCentBins+1] = {0, 10, 30, 50, 70, 90, 100};
Double_t EPm_resolution_v2_etagap[nCentBins] = {0.685732, 0.859684, 0.805492, 0.566930, 0.211378, 0.0307577};
Double_t EPp_resolution_v2_etagap[nCentBins] = {0.685895, 0.859866, 0.805762, 0.567147, 0.210694, 0.0329058};

Double_t minhisto=1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;
const int nMassBins = 60;
Double_t fstMassBin = minhisto;
Double_t widMassBin = binwidthmass;
Double_t massBins[nMassBins+1];

const int nDcaBins = 12;
Double_t fstDcaBin = 0.001;
Double_t widDcaBin = 1.27;
//Double_t dcaBins[nDcaBins+1] = {0, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.012, 0.016, 0.024, 0.032, 0.045, 0.070};
Double_t dcaBins[nDcaBins+1] = {0, 0.001, 0.00227, 0.0038829, 0.00593128, 0.008, 0.0118366, 0.0160324, 0.0213612, 0.0281287, 0.0367235, 0.0476388, 0.07};
const int nD0Bins = 20;
Double_t fstD0Bin = 3.5;
Double_t widD0Bin = 5;
Double_t d0Bins[nD0Bins+1];

TString tfname[3][2] = {{"v2_inpl","v2_outpl"},{"v3_inpl","v3_outpl"},{"inclusive",""}};

const int nFonllBins = 400;
Double_t fstFonllBins = 0;
Double_t lstFonllBins = 100;
Double_t widFonllBins = (lstFonllBins-fstFonllBins)/nFonllBins;

void initBins(bool resetDca=false)
{
  if(resetDca)
    {
      dcaBins[0] = 0;
      for(int i=1;i<nDcaBins+1;i++) dcaBins[i] = dcaBins[i-1]+fstDcaBin*pow(widDcaBin,i-1);
    }
  for(int i=0;i<nMassBins+1;i++) massBins[i] = fstMassBin+i*widMassBin;
  for(int i=0;i<nD0Bins+1;i++) d0Bins[i] = fstD0Bin+i*widD0Bin;
}

void DrawCmsTlatex(TString collision, Float_t tsize=0.04)
{
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} #bf{#it{Preliminary}}");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(tsize);
  texCms->SetTextFont(42);
  texCms->Draw();

  TLatex* texCol = new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV",collision.Data()));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(tsize);
  texCol->SetTextFont(42);
  texCol->Draw();
}

#endif
