using namespace std;
#ifndef _ANA_BFEEDDOWN_BFEEDDOWNFRACTION_H_
#define _ANA_BFEEDDOWN_BFEEDDOWNFRACTION_H_

#include "uti.h"
//#include "saveMassHisto.h"
TH1D* hD0DcaMCPSignal;
TH1D* hD0DcaMCNPSignal;

Double_t funMix(Double_t* x_, Double_t* para)
{
  float x = x_[0];
  float A = para[0];
  float R = para[1];
  float promptYield = 0;
  float nonPromptYield = 0;

  promptYield = hD0DcaMCPSignal->GetBinContent(hD0DcaMCPSignal->GetXaxis()->FindBin(x));
  nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));

  return A*(R*promptYield+(1-R)*nonPromptYield);
}

Double_t funNonPrompt(Double_t* x_, Double_t* para)
{
  float x = x_[0];
  float A = para[0];
  float R = para[1];
  float nonPromptYield = 0;
  nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
  return A*(1-R)*nonPromptYield;
}

void RandomSmear(TH1D* h0, TH1D* h)
{
  for(int i=1; i<h0->GetXaxis()->GetNbins(); i++)
    {
      h->SetBinContent(i, gRandom->Gaus(h0->GetBinContent(i), h0->GetBinError(i)));
    }
}

void normalize(TH1D* h)
{
  h->Sumw2();
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      Float_t val=h->GetBinContent(i);
      Float_t valErr=h->GetBinError(i);
      if (val==0) valErr=1;
      h->SetBinContent(i,val/h->GetBinWidth(i));
      h->SetBinError(i,valErr/h->GetBinWidth(i));
    }
  h->Scale(1./h->Integral(0,100,"width"));
}

void setColorTitleLabel(TH1* h, int color)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(20);
  h->SetStats(0);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
}

TF1* fitMass(TH1D* hData, TH1D* hMCSignal, TH1D* hMCSwapped)
{
  Double_t setparam0=100.;
  Double_t setparam1=1.865;
  Double_t setparam2=0.03;
  Double_t setparam10=0.005;
  Double_t setparam8=0.1;
  Double_t setparam9=0.1;
  Double_t fixparam1=1.865;
  Double_t minhisto=1.7;
  Double_t maxhisto=2.0;

  TF1* f = new TF1("fMass","[0]*([7]*([9]*Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11])))+(1-[7])*Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x", 1.7, 2.0);
  f->SetParLimits(4,-1000,1000);
  f->SetParLimits(10,0.005,0.05);
  f->SetParLimits(2,0.01,0.1);
  f->SetParLimits(8,0.02,0.2);
  f->SetParLimits(7,0,1);
  f->SetParLimits(9,0,1);

  f->SetParameter(0,setparam0);
  f->SetParameter(1,setparam1);
  f->SetParameter(2,setparam2);
  f->SetParameter(10,setparam10);
  f->SetParameter(9,setparam9);

  f->FixParameter(8,setparam8);
  f->FixParameter(7,1);
  f->FixParameter(1,fixparam1);
  f->FixParameter(3,0);
  f->FixParameter(4,0);
  f->FixParameter(5,0);
  f->FixParameter(6,0);
  f->FixParameter(11,0);

  hMCSignal->Fit("fMass","q","",minhisto,maxhisto);
  hMCSignal->Fit("fMass","q","",minhisto,maxhisto);
  f->ReleaseParameter(1);
  hMCSignal->Fit("fMass","L q","",minhisto,maxhisto);
  hMCSignal->Fit("fMass","L q","",minhisto,maxhisto);
  hMCSignal->Fit("fMass","L m","",minhisto,maxhisto);

  f->FixParameter(1,f->GetParameter(1));
  f->FixParameter(2,f->GetParameter(2));
  f->FixParameter(10,f->GetParameter(10));
  f->FixParameter(9,f->GetParameter(9));
  f->FixParameter(7,0);
  f->ReleaseParameter(8);
  f->SetParameter(8,setparam8);

  hMCSwapped->Fit("fMass","L q","",minhisto,maxhisto);
  hMCSwapped->Fit("fMass","L q","",minhisto,maxhisto);
  hMCSwapped->Fit("fMass","L q","",minhisto,maxhisto);
  hMCSwapped->Fit("fMass","L m","",minhisto,maxhisto);

  f->FixParameter(7,hMCSignal->Integral(0,1000)/(hMCSwapped->Integral(0,1000)+hMCSignal->Integral(0,1000)));
  f->FixParameter(8,f->GetParameter(8));
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);
  f->ReleaseParameter(6);

  f->SetLineColor(kRed);

  hData->Fit("fMass","q","",minhisto,maxhisto);
  hData->Fit("fMass","q","",minhisto,maxhisto);
  f->ReleaseParameter(1);
  f->SetParLimits(1,1.86,1.87);
  f->ReleaseParameter(11);
  f->SetParLimits(11,-0.2,0.2);
  hData->Fit("fMass","L q","",minhisto,maxhisto);
  hData->Fit("fMass","L q","",minhisto,maxhisto);
  hData->Fit("fMass","L q","",minhisto,maxhisto);
  hData->Fit("fMass","L m","",minhisto,maxhisto);

  TF1* background = new TF1("fBackground","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
  background->SetParameter(0,f->GetParameter(3));
  background->SetParameter(1,f->GetParameter(4));
  background->SetParameter(2,f->GetParameter(5));
  background->SetParameter(3,f->GetParameter(6));
  background->SetLineColor(4);
  background->SetRange(minhisto,maxhisto);
  background->SetLineStyle(2);

  TF1* mass = new TF1("fSignal","[0]*([3]*([4]*Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))))");
  mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(7),f->GetParameter(9),f->GetParameter(10),f->GetParameter(11));
  mass->SetParError(0,f->GetParError(0));
  mass->SetParError(1,f->GetParError(1));
  mass->SetParError(2,f->GetParError(2));
  mass->SetParError(3,f->GetParError(7));
  mass->SetParError(4,f->GetParError(9));
  mass->SetParError(5,f->GetParError(10));
  mass->SetFillColor(kOrange-3);
  mass->SetFillStyle(3002);
  mass->SetLineColor(kOrange-3);
  mass->SetLineWidth(3);
  mass->SetLineStyle(2);

  TF1* massSwap = new TF1("fBackground","[0]*(1-[2])*Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))");
  massSwap->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(7),f->GetParameter(8),f->GetParameter(11));
  massSwap->SetParError(0,f->GetParError(0));
  massSwap->SetParError(1,f->GetParError(1));
  massSwap->SetParError(2,f->GetParError(7));
  massSwap->SetParError(3,f->GetParError(8));
  massSwap->SetFillColor(kGreen+4);
  massSwap->SetFillStyle(3005);
  massSwap->SetLineColor(kGreen+4);
  massSwap->SetLineWidth(3);
  massSwap->SetLineStyle(1);

  hData->SetXTitle("m_{#piK} (GeV/c^{2})");
  hData->SetYTitle("Entries / (5 MeV/c^{2})");
  hData->SetAxisRange(0,hData->GetBinContent(hData->GetMaximumBin())*1.4*1.2,"Y");
  hData->SetMarkerSize(0.3);
  hData->Draw("e");

  cout<<"hData->GetMaximum(): "<<hData->GetMaximum()<<endl;

  background->Draw("same");
  mass->SetRange(minhisto,maxhisto);
  mass->Draw("same");
  massSwap->SetRange(minhisto,maxhisto);
  massSwap->Draw("same");
  f->Draw("same");

  Double_t yield = mass->Integral(minhisto,maxhisto)/hData->GetBinWidth(1);
  Double_t yieldErr = mass->Integral(minhisto,maxhisto)/hData->GetBinWidth(1)*mass->GetParError(0)/mass->GetParameter(0);

  std::cout<<"integral function yield: "<<yield<<"    fit yield: "<<f->GetParameter(0)*f->GetParameter(7)/hData->GetBinWidth(1)<<" +- "<<f->GetParError(0)*f->GetParameter(7)/hData->GetBinWidth(1)<<std::endl;

  TLegend* leg = new TLegend(0.65,0.5,0.82,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(hData,"Data","pl");
  leg->AddEntry(f,"Fit","l");
  leg->AddEntry(mass,"D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
  leg->AddEntry(massSwap,"K-#pi swapped","f");
  leg->AddEntry(background,"Combinatorial","l");
  leg->Draw("same");

  hData->GetFunction("fMass")->Delete();
  TH1D* hDataNoFitFun = (TH1D*) hData->Clone("hDataNoFitFun");
  hDataNoFitFun->Draw("esame");

  return f;
}

#endif
