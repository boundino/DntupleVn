using namespace std;
#include "../include/uti.h"
#include "../include/parameters.h"

void plotPnNP(TString outputfileP,TString outputfileNP,TString tfend,Float_t centmin,Float_t centmax)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TFile* infP = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputfileP.Data(),centmin,centmax,tfend.Data()));
  TFile* infNP = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputfileNP.Data(),centmin,centmax,tfend.Data()));

  TH1D* hEffP = (TH1D*)infP->Get("hEff");          hEffP->SetName("hEffP");
  TH1D* hEffNP = (TH1D*)infNP->Get("hEff");        hEffNP->SetName("hEffNP");

  hEffP->SetTitle(";p_{T} (GeV/c);#alpha x #epsilon_{reco} x #epsilon_{sel}");
  hEffP->SetMinimum(0.);
  hEffP->SetMaximum(1.5);
  hEffP->SetLineColor(2);
  hEffP->SetMarkerStyle(20);
  hEffP->SetMarkerSize(1.);
  hEffP->SetMarkerColor(2);
  hEffP->GetXaxis()->SetTitleOffset(0.95);//0.9
  hEffP->GetYaxis()->SetTitleOffset(1.15);//1.
  hEffP->GetXaxis()->SetTitleSize(0.055);//0.045
  hEffP->GetYaxis()->SetTitleSize(0.055);//0.045
  hEffP->GetXaxis()->SetTitleFont(42);
  hEffP->GetYaxis()->SetTitleFont(42);
  hEffNP->SetMarkerStyle(20);
  hEffNP->SetMarkerSize(1.);
  hEffNP->SetLineColor(4);
  hEffNP->SetMarkerColor(4);

  TString texper = "%";
  TLatex* texCent = new TLatex(0.57,0.84, Form("Centrality %.0f - %.0f%s",centmin,centmax,texper.Data()));
  texCent->SetNDC();
  texCent->SetTextAlign(12);
  texCent->SetTextSize(0.045);
  texCent->SetTextFont(42);

  TCanvas* cEff = new TCanvas("cEff","",600,600);
  hEffP->Draw();
  hEffNP->Draw("same");
  DrawCmsTlatex("PbPb");
  texCent->Draw();
  TLegend* legEff = new TLegend(0.55,0.70,0.90,0.81);
  legEff->SetBorderSize(0);
  legEff->SetFillStyle(0);
  legEff->AddEntry(hEffP,"Prompt","pl");
  legEff->AddEntry(hEffNP,"Non-prompt","pl");
  legEff->Draw("same");
  cEff->SaveAs(Form("plots/cEff_cent_%.0f_%.0f_%s.pdf",centmin,centmax,tfend.Data()));
}

int main(int argc, char* argv[])
{
  if(argc==6)
    {
      plotPnNP(argv[1],argv[2],argv[3],atof(argv[4]),atof(argv[5]));
      return 1;
    }
  else
    {
      std::cout<<"Invalid parameter number"<<std::endl;
      return 0;
    }
}
