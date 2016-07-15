using namespace std;
#include "../include/uti.h"
#include "../include/parameters.h"
Int_t pcolor[4] = {kBlue,kBlue-2,kGray+3,kGray+1};
TString tfend[4] = {"MC","MC_DCA","data","data_DCA"};

void plotFractions(TString inputfrMC, TString inputfrDa, Float_t centmin, Float_t centmax)
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

  TFile** inputfileMC = new TFile*[2];
  inputfileMC[0] = new TFile(Form("%s_cent_%.0f_%.0f_inclusive.root",inputfrMC.Data(),centmin,centmax));
  inputfileMC[1] = new TFile(Form("%s_cent_%.0f_%.0f_DCA.root",inputfrMC.Data(),centmin,centmax));
  TFile* inputfileDa = new TFile(Form("%s_cent_%.0f_%.0f_inclusive.root",inputfrDa.Data(),centmin,centmax));
  cout<<Form("%s_cent_%.0f_%.0f_inclusive.root",inputfrDa.Data(),centmin,centmax)<<endl;
  TGraphAsymmErrors** grFraction = new TGraphAsymmErrors*[2];
  for(int i=0;i<2;i++)
    {
      grFraction[i] = (TGraphAsymmErrors*)inputfileMC[i]->Get("grPromptFraction");
      grFraction[i]->SetName(Form("grPromptFraction_%s",tfend[i].Data()));
      grFraction[i]->SetMarkerSize(1.1);
      grFraction[i]->SetMarkerStyle(20);
      grFraction[i]->SetLineColor(pcolor[i]);
      grFraction[i]->SetMarkerColor(pcolor[i]);
    }
  
  TGraphErrors** gFraction = new TGraphErrors*[2];
  gFraction[0] = (TGraphErrors*)inputfileDa->Get("grPromptFraction");
  gFraction[1] = (TGraphErrors*)inputfileDa->Get("grPromptFraction_DCA");
  for(int i=0;i<2;i++)
    {
      gFraction[i]->SetName(Form("grPromptFraction_%s",tfend[i+2].Data()));
      gFraction[i]->SetMarkerSize(1.1);
      gFraction[i]->SetMarkerStyle(20);
      gFraction[i]->SetLineColor(pcolor[i+2]);
      gFraction[i]->SetMarkerColor(pcolor[i+2]);
    }
  
  TH2F* hempty = new TH2F("hempty","",20,0.,42.,10.,0.,1.);
  hempty->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
  hempty->GetYaxis()->SetTitle("Prompt fraction");
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.3);
  hempty->GetYaxis()->SetTitleOffset(1.8);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);
  TString tper = "%";
  TLatex* texCent = new TLatex(0.32,0.58, Form("Cent. %.0f-%.0f%s",centmin,centmax,tper.Data()));
  texCent->SetNDC();
  texCent->SetTextAlign(12);
  texCent->SetTextSize(0.04);
  texCent->SetTextFont(42);

  TCanvas* cMC = new TCanvas("cMC","",600,600);
  hempty->Draw();
  grFraction[0]->Draw("psame");
  grFraction[1]->Draw("psame");
  texCent->Draw();
  DrawCmsTlatex("PbPb");
  TLegend* legMC = new TLegend(0.30, 0.25, 0.80, 0.55);
  legMC->SetBorderSize(0);
  legMC->SetTextSize(0.04);
  legMC->SetTextFont(42);
  legMC->SetFillStyle(0);
  legMC->AddEntry((TObject*)0, "FONLL + MC eff", NULL);
  legMC->AddEntry(grFraction[0], "Inclusive", "pl");
  legMC->AddEntry(grFraction[1], "DCA < 0.008 cm", "pl");
  legMC->Draw();
  cMC->SaveAs(Form("plotsResult/cPromptFraction_FONLL_cent_%.0f_%.0f.pdf",centmin,centmax));

}

int main(int argc, char *argv[])
{
  if(argc==5)
    {
      plotFractions(argv[1],argv[2],atof(argv[3]), atof(argv[4]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
