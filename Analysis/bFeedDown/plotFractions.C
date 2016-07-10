using namespace std;

#include "uti.h"
#include "saveMassHisto.h"
Int_t pcolor[5] = {kBlack,kRed,kBlue,kRed,kBlue};
TString tfend[5] = {"inclusive","v2_inpl","v2_outpl","v3_inpl","v3_outpl"};
void plotFractions(Float_t centmin, Float_t centmax)
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

  TFile* inputfiles[5];
  for(int i=0;i<5;i++) inputfiles[i] = new TFile(Form("outfilesResult/bFeedDownResult_cent_%.0f_%.0f_%s.root",centmin,centmax,tfend[i].Data()));
  TGraphErrors* grFraction[5];
  for(int i=0;i<5;i++)
    {
      grFraction[i] = (TGraphErrors*)inputfiles[i]->Get("grPromptFraction");
      grFraction[i]->SetName(Form("grFraction_%d",i));
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
  for(int i=0;i<5;i++)
    {
      grFraction[i]->SetMarkerSize(1.1);
      grFraction[i]->SetMarkerStyle(21);
      grFraction[i]->SetLineColor(pcolor[i]);
      grFraction[i]->SetMarkerColor(pcolor[i]);
    }

  TString tper = "%";
  TLatex* texCent = new TLatex(0.51,0.57, Form("Cent. %.0f-%.0f%s",centmin,centmax,tper.Data()));
  texCent->SetNDC();
  texCent->SetTextAlign(12);
  texCent->SetTextSize(0.04);
  texCent->SetTextFont(42);

  TCanvas* cV2 = new TCanvas("cV2","",600,600);
  hempty->Draw();
  DrawCmsTlatex("PbPb");
  for(int i=0;i<5;i++)
    {
      if(i==3||i==4) continue;
      grFraction[i]->Draw("psame");
    }
  TLegend* legV2 = new TLegend(0.5, 0.40, 0.8, 0.55);
  legV2->SetBorderSize(0);
  legV2->SetTextSize(0.04);
  legV2->SetTextFont(42);
  legV2->SetFillStyle(0);
  legV2->AddEntry(grFraction[0], "Inclusive", "pl");
  legV2->AddEntry(grFraction[1], "v_{2} in plain", "pl");
  legV2->AddEntry(grFraction[2], "v_{2} out of plain", "pl");
  legV2->Draw();
  texCent->Draw();
  cV2->SaveAs(Form("plotsResult/cV2_cent_%.0f_%.0f.pdf",centmin,centmax));

  TCanvas* cV3 = new TCanvas("cV3","",600,600);
  hempty->Draw();
  DrawCmsTlatex("PbPb");
  for(int i=0;i<5;i++)
    {
      if(i==1||i==2) continue;
      grFraction[i]->Draw("psame");
    }
  TLegend* legV3 = new TLegend(0.5, 0.4, 0.8, 0.55);
  legV3->SetBorderSize(0);
  legV3->SetTextSize(0.04);
  legV3->SetTextFont(42);
  legV3->SetFillStyle(0);
  legV3->AddEntry(grFraction[0], "Inclusive", "pl");
  legV3->AddEntry(grFraction[1], "v_{3} in plain", "pl");
  legV3->AddEntry(grFraction[2], "v_{3} out of plain", "pl");
  legV3->Draw();
  texCent->Draw();
  cV3->SaveAs(Form("plotsResult/cV3_cent_%.0f_%.0f.pdf",centmin,centmax));
}

int main(int argc, char *argv[])
{
  if(argc==3)
    {
      plotFractions(atof(argv[1]), atof(argv[2]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
