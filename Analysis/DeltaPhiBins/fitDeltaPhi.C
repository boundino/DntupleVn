using namespace std;
#include "saveMassHisto.h"

TString infname;
Float_t centMin,centMax;
int fitDeltaPhi(TString inputfile="outfiles/DeltaPhiHisto", TString outputfile="outfiles/V2PtHisto",
                Float_t centmin=0., Float_t centmax=100.)
{
  infname = inputfile;
  centMin = centmin;
  centMax = centmax;

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TF1* fit(Float_t ptmin, Float_t ptmax);

  TH1D* hV2 = new TH1D("hV2",";D^{0} p_{T} (GeV/c);v_{2}",nPtBins,ptBins);
  Int_t centbin = findcentbin(centmin,centmax);
  if(centbin<0) return 1;
  Float_t Rn = (EPm_resolution_v2_etagap[centbin]+EPp_resolution_v2_etagap[centbin])/2.;
  for(int i=0;i<nPtBins;i++)
    {
      TF1* f = fit(ptBins[i],ptBins[i+1]);
      Double_t v2 = f->GetParameter(1) / Rn;
      Double_t v2Err = f->GetParError(1) /Rn;
      hV2->SetBinContent(i+1,v2);
      hV2->SetBinError(i+1,v2Err);
    }
  TCanvas* cV2 = new TCanvas("cV2","",600,600);
  hV2->SetXTitle("D^{0} p_{T} / GeV/c");
  hV2->SetYTitle("v_{2}");
  hV2->GetXaxis()->CenterTitle();
  hV2->GetYaxis()->CenterTitle();
  hV2->SetAxisRange(-0.2,0.4,"Y");
  hV2->GetXaxis()->SetTitleOffset(1.1);
  hV2->GetYaxis()->SetTitleOffset(1.2);
  hV2->GetXaxis()->SetLabelOffset(0.007);
  hV2->GetYaxis()->SetLabelOffset(0.007);
  hV2->GetXaxis()->SetTitleSize(0.050);
  hV2->GetYaxis()->SetTitleSize(0.050);
  hV2->GetXaxis()->SetTitleFont(42);
  hV2->GetYaxis()->SetTitleFont(42);
  hV2->GetXaxis()->SetLabelFont(42);
  hV2->GetYaxis()->SetLabelFont(42);
  hV2->GetXaxis()->SetLabelSize(0.04);
  hV2->GetYaxis()->SetLabelSize(0.04);
  hV2->SetLineColor(kBlack);
  hV2->SetMarkerColor(kBlack);
  hV2->SetMarkerSize(1.0);
  hV2->SetMarkerStyle(20);
  hV2->SetStats(0);
  hV2->Draw("e");

  TLine* l = new TLine(2., 0., 40., 0.);
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);
  l->Draw();

  DrawCmsTlatex("PbPb");
  TLatex* tex;

  tex = new TLatex(0.58,0.83,"|y| < 1.0");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  TString texper="%";
  tex = new TLatex(0.58,0.78,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.58,0.73,"#Delta#Phi method");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  cV2->SaveAs(Form("plots/V2_deltaphibins_cent_%.0f_%.0f.pdf",centMin,centMax));
  cV2->SaveAs(Form("plots/V2_deltaphibins_cent_%.0f_%.0f.png",centMin,centMax));
  TFile* outf = new TFile(Form("%s_cent_%.0f_%.0f.root",outputfile.Data(),centmin,centmax),"recreate");
  outf->cd();
  hV2->Write();
  outf->Close();

  return 0;
}

TF1* fit(Float_t ptmin, Float_t ptmax)
{
  TCanvas* c = new TCanvas(Form("c_%.0f_%.0f",ptmin,ptmax),"",600,600);
  TFile* infile = new TFile(Form("%s_cent_%.0f_%.0f_pt_%.0f_%.0f.root",infname.Data(),centMin,centMax,ptmin,ptmax));
  TH1D* hPhi = (TH1D*)infile->Get("hPhi");

  TF1* f = new TF1(Form("f_%.0f_%.0f",ptmin,ptmax),"[0]*(1+2*[1]*TMath::Cos(2*x))",0,PI/2.);
  f->SetParLimits(0,0,1.e+10);
  f->SetParLimits(1,0,1.);
  f->SetParameter(0.,1.e+3);
  f->SetParameter(1.,0.1);
  f->SetLineColor(kRed);

  hPhi->Fit(Form("f_%.0f_%.0f",ptmin,ptmax),"I m","",0.,PI/2.);
  hPhi->Fit(Form("f_%.0f_%.0f",ptmin,ptmax),"I m","",0.,PI/2.);

  hPhi->SetXTitle("#Delta#Phi");
  hPhi->SetYTitle("d^{2}N / (dp_{T}d#Delta#Phi)");
  hPhi->GetXaxis()->CenterTitle();
  hPhi->GetYaxis()->CenterTitle();
  Float_t vhMin = hPhi->GetMinimum();
  Float_t vhMax = hPhi->GetMaximum();
  Float_t vhIntev = vhMax-vhMin;
  hPhi->SetAxisRange(vhMin-vhIntev*0.25,vhMax+vhIntev*0.2,"Y");
  hPhi->GetXaxis()->SetTitleOffset(1.3);
  hPhi->GetYaxis()->SetTitleOffset(1.8);
  hPhi->GetXaxis()->SetLabelOffset(0.007);
  hPhi->GetYaxis()->SetLabelOffset(0.007);
  hPhi->GetXaxis()->SetTitleSize(0.045);
  hPhi->GetYaxis()->SetTitleSize(0.045);
  hPhi->GetXaxis()->SetTitleFont(42);
  hPhi->GetYaxis()->SetTitleFont(42);
  hPhi->GetXaxis()->SetLabelFont(42);
  hPhi->GetYaxis()->SetLabelFont(42);
  hPhi->GetXaxis()->SetLabelSize(0.04);
  hPhi->GetYaxis()->SetLabelSize(0.04);
  hPhi->SetLineColor(kBlack);
  hPhi->SetMarkerColor(kBlack);
  hPhi->SetMarkerSize(0.8);
  hPhi->SetMarkerStyle(20);
  hPhi->SetStats(0);
  hPhi->Draw("e");
  f->Draw("same");

  Float_t v2ob = f->GetParameter(1);
  Float_t v2obErr = f->GetParError(1);

  DrawCmsTlatex("PbPb");
  TLatex* tex;

  tex = new TLatex(0.58,0.83,"|y| < 1.0");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.58,0.78,Form("%.1f < p_{T} < %.1f GeV/c",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  TString texper="%";
  tex = new TLatex(0.58,0.72,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.58,0.67,Form("v_{2}^{obs} = %.3f #pm %.3f",v2ob,v2obErr));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  c->SaveAs(Form("plots/DeltaPhi_cent_%.0f_%.0f_pt_%.0f_%.0f.pdf",centMin,centMax,ptmin,ptmax));

  return f;

}

int main(int argc, char *argv[])
{
  if(argc==5)
    {
      fitDeltaPhi(argv[1], argv[2], atof(argv[3]), atof(argv[4]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
