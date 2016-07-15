using namespace std;
#include "../include/uti.h"
#include "../include/parameters.h"

void bFeedDownFONLL(TString outputFonllP, TString outputFonllNP, TString outputEffP, TString outputEffNP, TString outputFraction, TString tfend, Float_t centmin, Float_t centmax)
{
  TFile* inputFileFonllP = new TFile(outputFonllP);
  TFile* inputFileFonllNP = new TFile(outputFonllNP);
  TFile* inputFileEffP = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputEffP.Data(),centmin,centmax,tfend.Data()));
  TFile* inputFileEffNP = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputEffNP.Data(),centmin,centmax,tfend.Data()));

  TH1D* hEffP = (TH1D*)inputFileEffP->Get("hEff");                                          hEffP->SetName("hEffP");
  TH1D* hEffNP = (TH1D*)inputFileEffNP->Get("hEff");                                        hEffNP->SetName("hEffNP");
  TGraphAsymmErrors* gFonllP = (TGraphAsymmErrors*)inputFileFonllP->Get("gaeSigmaDzero");   gFonllP->SetName("gFonllP");
  TGraphAsymmErrors* gFonllNP = (TGraphAsymmErrors*)inputFileFonllNP->Get("gaeSigmaDzero"); gFonllNP->SetName("gFonllNP");
  TGraphAsymmErrors* grFraction = (TGraphAsymmErrors*)gFonllP->Clone("grPromptFraction");
  for(int i=0;i<nPtBins;i++)
    {
      Double_t effP = hEffP->GetBinContent(i+1);
      Double_t effNP = hEffNP->GetBinContent(i+1);
      Double_t effPErr = hEffP->GetBinError(i+1)/effP;
      Double_t effNPErr = hEffNP->GetBinError(i+1)/effNP;
      Double_t fonllP,fonllNP,pt;
      gFonllP->GetPoint(i,pt,fonllP);
      gFonllNP->GetPoint(i,pt,fonllNP);
      Double_t fonllPErrl = gFonllP->GetErrorYlow(i)/fonllP;
      Double_t fonllPErrh = gFonllP->GetErrorYhigh(i)/fonllP;
      Double_t fonllNPErrl = gFonllNP->GetErrorYlow(i)/fonllNP;
      Double_t fonllNPErrh = gFonllNP->GetErrorYhigh(i)/fonllNP;
      Double_t nfP = fonllP*effP;
      Double_t nfNP = fonllNP*effNP;
      Double_t nfPErrl = TMath::Sqrt(fonllPErrl*fonllPErrl+effPErr*effPErr);
      Double_t nfPErrh = TMath::Sqrt(fonllPErrh*fonllPErrh+effPErr*effPErr);
      Double_t nfNPErrl = TMath::Sqrt(fonllNPErrl*fonllNPErrl+effNPErr*effNPErr);
      Double_t nfNPErrh = TMath::Sqrt(fonllNPErrh*fonllNPErrh+effNPErr*effNPErr);
      Double_t pFraction = nfP/(nfP+nfNP);
      Double_t pFractionErrl = TMath::Sqrt(nfPErrl*nfPErrl+(nfPErrl*nfPErrl*nfNP*nfNP+nfNPErrh*nfNPErrh*nfP*nfP)/((nfP+nfNP)*(nfP+nfNP)));
      Double_t pFractionErrh = TMath::Sqrt(nfPErrh*nfPErrh+(nfPErrh*nfPErrh*nfNP*nfNP+nfNPErrl*nfNPErrl*nfP*nfP)/((nfP+nfNP)*(nfP+nfNP)));
      grFraction->SetPoint(i,pt,pFraction);
      grFraction->SetPointEYlow(i,0);
      grFraction->SetPointEYhigh(i,0);
      /*
      grFraction->SetPointEYlow(i,pFractionErrl*pFraction);
      grFraction->SetPointEYhigh(i,pFractionErrh*pFraction);
      */
    }
  TFile* fout = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputFraction.Data(),centmin,centmax,tfend.Data()), "recreate");
  grFraction->Write();
  fout->Close();

}

int main(int argc, char *argv[])
{
  if(argc==9)
    {
      bFeedDownFONLL(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],atof(argv[7]),atof(argv[8]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
