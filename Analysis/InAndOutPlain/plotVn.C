using namespace std;
#include "saveMassHisto.h"

TString infname;
Float_t centMin,centMax;
int plotVn(TString inputfile="outfiles/VnPtHisto",
          Float_t centmin=30., Float_t centmax=50.)
{
  infname = inputfile;
  centMin = centmin;
  centMax = centmax;

  void plotComparison(TH1D* hInclusive, TH1D* hPrompt);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TFile* infile_Inclusive = new TFile(Form("%s_cent_%.0f_%.0f_Inclusive.root",infname.Data(),centMin,centMax));
  TH1D* hV2_Inclusive = (TH1D*)infile_Inclusive->Get("hV2"); hV2_Inclusive->SetName("hV2_Inclusive");
  TH1D* hV3_Inclusive = (TH1D*)infile_Inclusive->Get("hV3"); hV3_Inclusive->SetName("hV3_Inclusive");
  TFile* infile_Prompt = new TFile(Form("%s_cent_%.0f_%.0f_Prompt.root",infname.Data(),centMin,centMax));
  TH1D* hV2_Prompt = (TH1D*)infile_Prompt->Get("hV2"); hV2_Prompt->SetName("hV2_Prompt");
  TH1D* hV3_Prompt = (TH1D*)infile_Prompt->Get("hV3"); hV3_Prompt->SetName("hV3_Prompt");

  TCanvas* cV2 = new TCanvas("cV2","",600,600);
  hV2_Inclusive->SetXTitle("D^{0} p_{T} / GeV/c");
  hV2_Inclusive->SetYTitle("v_{2}");
  plotComparison(hV2_Inclusive,hV2_Prompt);
  cV2->SaveAs(Form("plotsResult/V2_InAndOutPlain_cent_%.0f_%.0f.pdf",centMin,centMax));
  cV2->SaveAs(Form("plotsResult/V2_InAndOutPlain_cent_%.0f_%.0f.png",centMin,centMax));

  TCanvas* cV3 = new TCanvas("cV3","",600,600);
  hV3_Inclusive->SetXTitle("D^{0} p_{T} / GeV/c");
  hV3_Inclusive->SetYTitle("v_{3}");
  plotComparison(hV3_Inclusive,hV3_Prompt);
  cV3->SaveAs(Form("plotsResult/V3_InAndOutPlain_cent_%.0f_%.0f.pdf",centMin,centMax));
  cV3->SaveAs(Form("plotsResult/V3_InAndOutPlain_cent_%.0f_%.0f.png",centMin,centMax));

  return 0;
}

void plotComparison(TH1D* hInclusive, TH1D* hPrompt)
{
  TLine* l = new TLine(2., 0., 40., 0.);
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);

  hInclusive->GetXaxis()->CenterTitle();
  hInclusive->GetYaxis()->CenterTitle();
  hInclusive->SetAxisRange(-0.2,0.4,"Y");
  hInclusive->GetXaxis()->SetTitleOffset(1.1);
  hInclusive->GetYaxis()->SetTitleOffset(1.2);
  hInclusive->GetXaxis()->SetLabelOffset(0.007);
  hInclusive->GetYaxis()->SetLabelOffset(0.007);
  hInclusive->GetXaxis()->SetTitleSize(0.050);
  hInclusive->GetYaxis()->SetTitleSize(0.050);
  hInclusive->GetXaxis()->SetTitleFont(42);
  hInclusive->GetYaxis()->SetTitleFont(42);
  hInclusive->GetXaxis()->SetLabelFont(42);
  hInclusive->GetYaxis()->SetLabelFont(42);
  hInclusive->GetXaxis()->SetLabelSize(0.04);
  hInclusive->GetYaxis()->SetLabelSize(0.04);

  hInclusive->SetLineColor(kBlack);
  hInclusive->SetMarkerColor(kBlack);
  hInclusive->SetMarkerSize(1.0);
  hInclusive->SetMarkerStyle(20);
  hInclusive->SetStats(0);
  hInclusive->Draw("e");
  hPrompt->SetLineColor(kBlue);
  hPrompt->SetMarkerColor(kBlue);
  hPrompt->SetMarkerSize(1.0);
  hPrompt->SetMarkerStyle(20);
  hPrompt->SetStats(0);
  hPrompt->Draw("esame");

  l->Draw();
  DrawCmsTlatex("PbPb");
  TLegend* legV2 = new TLegend(0.5, 0.2, 0.8, 0.3);
  legV2->SetBorderSize(0);
  legV2->SetTextSize(0.04);
  legV2->SetTextFont(42);
  legV2->SetFillStyle(0);
  legV2->AddEntry(hInclusive, "Inclusive D^{0}", "pl");
  legV2->AddEntry(hPrompt, "Prompt D^{0}", "pl");
  legV2->Draw();

  TLatex* tex;

  tex = new TLatex(0.50,0.83,"|y| < 1.0");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  TString texper="%";
  tex = new TLatex(0.50,0.78,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.50,0.73,"In and out of plain method");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

}

int main(int argc, char *argv[])
{
  if(argc==4)
    {
      plotVn(argv[1], atof(argv[2]), atof(argv[3]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
