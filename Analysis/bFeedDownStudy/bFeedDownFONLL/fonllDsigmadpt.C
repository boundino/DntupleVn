using namespace std;
#include "../include/uti.h"
#include "../include/parameters.h"

void fonllDsigmadpt(TString inputFONLL="fonlls/FONLL_pp_promptDzero_5TeV_y1.dat", TString outputFONLL="outfiles/FONLL_pp_promptDzero_5TeV_y1.root")
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

  TString tlabel = inputFONLL;
  tlabel.ReplaceAll("fonlls/","");
  tlabel.ReplaceAll(".dat","");
  TString tp = inputFONLL;
  tp.ReplaceAll("fonlls/FONLL_pp_","");
  tp.ReplaceAll("Dzero_5TeV_y1.dat","");

  cout<<endl;
  cout<<"  -- Processing FONLL: "<<tlabel<<endl;

  Float_t DzeroFF=0.557;           //FF of D->D0
  ifstream getdata(Form("%s",inputFONLL.Data()));
  if(!getdata.is_open())
    {
      cout<<"Opening the file fails"<<endl;
    }

  Float_t central[nFonllBins];
  Float_t min_all[nFonllBins],max_all[nFonllBins],min_sc[nFonllBins],max_sc[nFonllBins],min_mass[nFonllBins],max_mass[nFonllBins];
  Float_t tem;
  for(int i=0;i<nFonllBins;i++)
    {
      getdata>>tem;
      getdata>>central[i];
      getdata>>min_all[i];
      getdata>>max_all[i];
      getdata>>min_sc[i];
      getdata>>max_sc[i];
      getdata>>min_mass[i];
      getdata>>max_mass[i];
    }

  TH1F* hpt = new TH1F("hpt","",nFonllBins,fstFonllBins,lstFonllBins);
  TH1F* hminall = new TH1F("hminall","",nFonllBins,fstFonllBins,lstFonllBins);
  TH1F* hmaxall = new TH1F("hmaxall","",nFonllBins,fstFonllBins,lstFonllBins);
  for(int i=0;i<nFonllBins;i++)
    {
      hpt->SetBinContent(i+1,central[i]);
      hminall->SetBinContent(i+1,min_all[i]);
      hmaxall->SetBinContent(i+1,max_all[i]);
    }
  TH1F* hpt_rebin = (TH1F*)hpt->Rebin(nPtBins,"hpt_rebin",ptBins);
  TH1F* hminall_rebin = (TH1F*)hminall->Rebin(nPtBins,"hminall_rebin",ptBins);
  TH1F* hmaxall_rebin = (TH1F*)hmaxall->Rebin(nPtBins,"hmaxall_rebin",ptBins);
  hpt_rebin->Scale(widFonllBins);
  hminall_rebin->Scale(widFonllBins);
  hmaxall_rebin->Scale(widFonllBins);
  divideBinWidth(hpt_rebin,false);
  divideBinWidth(hminall_rebin,false);
  divideBinWidth(hmaxall_rebin,false);

  Float_t asigma[nPtBins],aerrorl[nPtBins],aerrorh[nPtBins]; 
  Float_t apt[nPtBins],aptl[nPtBins];
  for(int i=0;i<nPtBins;i++)
    {
      apt[i] = (ptBins[i+1]+ptBins[i])/2.;
      aptl[i] = (ptBins[i+1]-ptBins[i])/2.;
      asigma[i] = hpt_rebin->GetBinContent(i+1);
      aerrorl[i] = hpt_rebin->GetBinContent(i+1)-hminall_rebin->GetBinContent(i+1);
      aerrorh[i] = hmaxall_rebin->GetBinContent(i+1)-hpt_rebin->GetBinContent(i+1);
    }

  TGraphAsymmErrors* gaeSigma = new TGraphAsymmErrors(nPtBins, apt, asigma, aptl, aptl, aerrorl, aerrorh);
  gaeSigma->SetFillColor(2);
  gaeSigma->SetFillStyle(3001);

  TGraphAsymmErrors* gaeSigmaDzero=(TGraphAsymmErrors*)gaeSigma->Clone();
  gaeSigmaDzero->SetName("gaeSigmaDzero");
  gaeSigmaDzero->SetFillColor(2);
  gaeSigmaDzero->SetFillStyle(3001);
  
  for(int i=0;i<gaeSigmaDzero->GetN();i++)
    {
      gaeSigmaDzero->GetY()[i] *= DzeroFF;
      gaeSigmaDzero->SetPointEYhigh(i,gaeSigmaDzero->GetErrorYhigh(i)*DzeroFF);
      gaeSigmaDzero->SetPointEYlow(i,gaeSigmaDzero->GetErrorYlow(i)*DzeroFF); 
    }
     
  TCanvas* cFonll = new TCanvas("cFonll","",600,500);
  cFonll->SetLogy();
  TH2F* hempty = new TH2F("hempty","",10,0,50.,10,1.e+2,1.e+11);
  hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hempty->GetYaxis()->SetTitle("d#sigma(D)/dp_{T} (pb#upointGeV^{-1}c)");
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.);
  hempty->GetYaxis()->SetTitleOffset(1.2);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);  
  hempty->Draw();
  hminall->SetLineColor(4);
  hmaxall->SetLineColor(4);
  hpt->SetLineColor(4);
  hminall->Draw("same");
  hmaxall->Draw("same");
  hpt->Draw("same");
  gaeSigma->SetLineWidth(3);
  gaeSigma->Draw("psame");
  gaeSigmaDzero->SetLineWidth(3);
  gaeSigmaDzero->SetLineColor(2);
  gaeSigmaDzero->Draw("psame");

  TLatex* tlatex = new TLatex(0.37,0.85,"pp 5TeV FONLL calculation");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.042);
  tlatex->Draw();
  
  TLatex* tlatextotunc = new TLatex(0.37,0.80,"Total syst uncertainties shown");
  tlatextotunc->SetNDC();
  tlatextotunc->SetTextColor(1);
  tlatextotunc->SetTextFont(42);
  tlatextotunc->SetTextSize(0.042);
  tlatextotunc->Draw();
  
  TLatex* tlatexD0 = new TLatex(0.37,0.75,"D^{0},|y|<1, BR unc not shown");
  tlatexD0->SetNDC();
  tlatexD0->SetTextColor(1);
  tlatexD0->SetTextFont(42);
  tlatexD0->SetTextSize(0.042);
  tlatexD0->Draw();
  
  TLegend* leg = new TLegend(0.35,0.60,0.90,0.70);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.042);
  TLegendEntry* ent_gaeSigma = leg->AddEntry(gaeSigma,"Frag.Frac.=1.0 (pure FONLL)","PL");
  TLegendEntry* ent_gaeSigmaDzero = leg->AddEntry(gaeSigmaDzero,"Multiplied by Frag.Frac.=0.577","PL");
  leg->Draw();

  TString tprompt = "Prompt D meson";
  if(tp=="nonprompt") tprompt = "B to D meson";
  TLatex* tPrompt = new TLatex(0.37,0.55,tprompt);
  tPrompt->SetNDC();
  tPrompt->SetTextColor(1);
  tPrompt->SetTextFont(42);
  tPrompt->SetTextSize(0.042);
  tPrompt->Draw();

  gaeSigma->SetName("gaeSigma");
  gaeSigmaDzero->SetName("gaeSigmaDzero");
  cFonll->SaveAs(Form("plots/c%s.pdf",tlabel.Data()));
  
  TFile* foutput = new TFile(outputFONLL.Data(),"recreate");
  foutput->cd();
  gaeSigma->Write();
  gaeSigmaDzero->Write();
  hpt->Write();
  hminall->Write();
  hmaxall->Write();
  foutput->Close();
}

int main(int argc, char *argv[])
{
  if(argc==3)
    {
      fonllDsigmadpt(argv[1], argv[2]);
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
