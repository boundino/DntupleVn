using namespace std;
#include "saveMassHisto.h"

TString infname;
Float_t centMin,centMax;
int calVn(TString inputfile="outfiles/YieldHisto", TString outputfile="outfiles/VnPtHisto",
          Float_t centmin=0., Float_t centmax=100.,
          Int_t isPrompt=0, TString tfprFraction="../bFeedDown/outfilesResult/bFeedDownResult")
{
  infname = inputfile;
  centMin = centmin;
  centMax = centmax;

  void calvn(Float_t Nin, Float_t NinErr, Float_t Nout, Float_t NoutErr, Float_t Rn, Float_t* results);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TH1D* hV2 = new TH1D("hV2",";D^{0} p_{T} (GeV/c);v_{2}",nPtBins,ptBins);
  TH1D* hV3 = new TH1D("hV3",";D^{0} p_{T} (GeV/c);v_{3}",nPtBins,ptBins);

  Int_t centbin = findcentbin(centmin,centmax);
  if(centbin<0) return 1;
  Float_t Rn2 = (EPm_resolution_v2_etagap[centbin]+EPp_resolution_v2_etagap[centbin])/2.;
  Float_t Rn3 = (EPm_resolution_v3_etagap[centbin]+EPp_resolution_v3_etagap[centbin])/2.;

  TFile* infile = new TFile(Form("%s_cent_%.0f_%.0f.root",infname.Data(),centMin,centMax));
  TH1D** hYield = new TH1D*[4];
  Double_t temptBins;
  Double_t prFraction[4][nPtBins];
  for(int i=0;i<4;i++) 
    {
      hYield[i] = (TH1D*)infile->Get(Form("hYield_%s",tfend[i].Data()));
      if(isPrompt)
        {
          TFile* infprFraction = new TFile(Form("%s_%s.root",tfprFraction.Data(),tfend[i].Data()));
          TGraphErrors* grFraction = (TGraphErrors*)infprFraction->Get("grPromptFraction");
          for(int j=0;j<nPtBins;j++) grFraction->GetPoint(j,temptBins,prFraction[i][j]);
          delete grFraction;
          delete infprFraction;
        }
      else
        {
          for(int j=0;j<nPtBins;j++) prFraction[i][j] = 1.;
        }
    }
  //
  for(int i=0;i<nPtBins;i++)
    {
      Float_t Nyield[4];
      Float_t NyieldErr[4];
      for(int j=0;j<4;j++)
        {
          Nyield[j] = hYield[j]->GetBinContent(i+1);
          NyieldErr[j] = hYield[j]->GetBinError(i+1);
        }
      Float_t* V2results = new Float_t[2];
      Float_t* V3results = new Float_t[2];
      calvn(Nyield[0]*prFraction[0][i],NyieldErr[0]*prFraction[0][i],Nyield[1]*prFraction[1][i],NyieldErr[1]*prFraction[1][i],Rn2,V2results);
      calvn(Nyield[2]*prFraction[2][i],NyieldErr[2]*prFraction[2][i],Nyield[3]*prFraction[3][i],NyieldErr[3]*prFraction[3][i],Rn3,V3results);
      hV2->SetBinContent(i+1,V2results[0]);
      hV2->SetBinError(i+1,V2results[1]);
      hV3->SetBinContent(i+1,V3results[0]);
      hV3->SetBinError(i+1,V3results[1]);
    }

  TFile* outf = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputfile.Data(),centmin,centmax,tpr[isPrompt].Data()),"recreate");
  outf->cd();
  hV2->Write();
  hV3->Write();
  outf->Close();

  return 0;
}

void calvn(Float_t Nin, Float_t NinErr, Float_t Nout, Float_t NoutErr, Float_t Rn, Float_t* results)
{
  Float_t vn,vnErr;
  vn = (PI/4.)*(1./Rn)*((Nin-Nout)/(Nin+Nout));
  vnErr = vn*(NinErr+NoutErr)*(1./TMath::Abs(Nin-Nout)+1./TMath::Abs(Nin+Nout));
  results[0] = vn;
  results[1] = vnErr;
}


int main(int argc, char *argv[])
{
  if(argc==7)
    {
      calVn(argv[1], argv[2], atof(argv[3]), atof(argv[4]), atoi(argv[5]), argv[6]);
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
