using namespace std;
#include "uti.h"
#include "parameters.h"

Double_t minhisto=1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

TString weight = "pthatweight";
TString weightgen;
TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
Float_t hiBinMin,hiBinMax,centMin,centMax;

#define PI 3.1415926

void fitD(TString inputdata="", TString inputmc="",
          TString trgselection="",  TString cut="", TString cutmcgen="",
          int doweight=0, TString outputfile="outfiles/mytest",
          Float_t centmin=0., Float_t centmax=100.)
{
  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  seldata = Form("%s&&%s&&hiBin>%f&&hiBin<%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
  selmceff = Form("%s&&hiBin>%f&&hiBin<%f",cut.Data(),hiBinMin,hiBinMax);
  selmcgen = Form("%s&&hiBin>%f&&hiBin<%f",cutmcgen.Data(),hiBinMin,hiBinMax);
  selmc = Form("%s",cut.Data());

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TF1* fit(TTree* nt, TTree* ntMC, Float_t ptmin, Float_t ptmax, Float_t dphimin, Float_t dphimax, Double_t count);

  if(doweight==0)
    {
      weightgen="1";
      weight="1";
    }
  
  TFile* inf = new TFile(inputdata.Data());
  TFile* infMC = new TFile(inputmc.Data());

  TTree* nt = (TTree*)inf->Get("ntDkpi");
  nt->AddFriend("ntHlt");
  nt->AddFriend("ntHi");
  TTree* ntMC = (TTree*)infMC->Get("ntDkpi");
  ntMC->AddFriend("ntHlt");
  ntMC->AddFriend("ntHi");

  Int_t nPhiBins = 5;
  Float_t dphimi=0,dphima=PI/2.;
  for(int i=0;i<nBins;i++)
    {
      TH1D* hPhi = new TH1D("hPhi","",nPhiBins,0,PI/2.);
      for(int j=0;j<nPhiBins;j++)
        {
          dphimi = (j*1./nPhiBins)*(PI/2.);
          dphima = ((j+1)*1./nPhiBins)*(PI/2.);
          TF1* f = fit(nt, ntMC, ptBins[i], ptBins[i+1], dphimi, dphima, j);
          double yield = f->Integral(minhisto,maxhisto)/binwidthmass;
          double yieldErr = f->Integral(minhisto,maxhisto)/binwidthmass*f->GetParError(0)/f->GetParameter(0);
          hPhi->SetBinContent(j+1,yield/((dphima-dphimi)*(ptBins[i+1]-ptBins[i])));
          hPhi->SetBinError(j+1,yieldErr/((dphima-dphimi)*(ptBins[i+1]-ptBins[i])));
        }
      TFile* outf = new TFile(Form("%s_cent_%.0f_%.0f_pt_%.0f_%.0f.root",outputfile.Data(),centMin,centMax,ptBins[i],ptBins[i+1]),"recreate");
      outf->cd();
      hPhi->Write();
      outf->Close();
      delete h;
      delete outf;
    }
  cout<<endl;
}

TF1* fit(TTree* nt, TTree* ntMC, Float_t ptmin, Float_t ptmax, Float_t dphimin, Float_t dphimax, Double_t count)
{
  TCanvas* c = new TCanvas(Form("c%d",count),"",600,600);
  TH1D* h = new TH1D(Form("h-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  TH1D* hMCSignal = new TH1D(Form("hMCSignal-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  TH1D* hMCSwapped = new TH1D(Form("hMCSwapped-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  
  TF1* f = new TF1(Form("f%d",count),"[0]*([7]*([9]*Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11])))+(1-[7])*Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x", 1.7, 2.0);

  TString dphipo = "abs(Dphi-hiEvtPlanes[0])";
  TString dphine = "abs(Dphi-hiEvtPlanes[1])";
  TString cutdphipo = Form("(%s>0&&%s<1.570796&&%s>%f&&%s<%f)||(%s>1.570796&&%s<3.141593&&(3.141593-%s)>%f&&(3.141593-%s)<%f)||(%s>3.141593&&%s<4.712389&&(%s-3.141593)>%f&&(%s-3.141593)<%f)",dphipo.Data(),dphipo.Data(),dphipo.Data(),dphimin,dphipo.Data(),dphimax,dphipo.Data(),dphipo.Data(),dphipo.Data(),dphimin,dphipo.Data(),dphimax,dphipo.Data(),dphipo.Data(),dphipo.Data(),dphimin,dphipo.Data(),dphimax);
  TString cutdphine = Form("(%s>0&&%s<1.570796&&%s>%f&&%s<%f)||(%s>1.570796&&%s<3.141593&&(3.141593-%s)>%f&&(3.141593-%s)<%f)||(%s>3.141593&&%s<4.712389&&(%s-3.141593)>%f&&(%s-3.141593)<%f)",dphine.Data(),dphine.Data(),dphine.Data(),dphimin,dphine.Data(),dphimax,dphine.Data(),dphine.Data(),dphine.Data(),dphimin,dphine.Data(),dphimax,dphine.Data(),dphine.Data(),dphine.Data(),dphimin,dphine.Data(),dphimax);
  
  nt->Project(Form("h-%d",count),"Dmass",Form("(%s&&Dpt>%f&&Dpt<%f&&((Dy>0&&(%s))||(Dy<0&&(%s))))",seldata.Data(),ptmin,ptmax,cutdphipo.Data(),cutdphine.Data()));////////////////////////////////////////////////////////////delta phi branch name
  ntMC->Project(Form("hMCSignal-%d",count),"Dmass",Form("%s*(%s&&Dpt>%f&&Dpt<%f&&(Dgen==23333))","1",selmc.Data(),ptmin,ptmax));   
  ntMC->Project(Form("hMCSwapped-%d",count),"Dmass",Form("%s*(%s&&Dpt>%f&&Dpt<%f&&(Dgen==23344))","1",selmc.Data(),ptmin,ptmax));
  
  f->SetParLimits(4,-1000,1000);
  f->SetParLimits(10,0.001,0.05);
  f->SetParLimits(2,0.01,0.5);
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
  h->GetEntries();
  
  hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  f->ReleaseParameter(1);
  hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
  
  f->FixParameter(1,f->GetParameter(1));
  f->FixParameter(2,f->GetParameter(2));
  f->FixParameter(10,f->GetParameter(10));
  f->FixParameter(9,f->GetParameter(9));
  f->FixParameter(7,0);
  f->ReleaseParameter(8);
  f->SetParameter(8,setparam8);
  
  hMCSwapped->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSwapped->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSwapped->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSwapped->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
  
  f->FixParameter(7,hMCSignal->Integral(0,1000)/(hMCSwapped->Integral(0,1000)+hMCSignal->Integral(0,1000)));
  f->FixParameter(8,f->GetParameter(8));
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);
  f->ReleaseParameter(6);

  f->SetLineColor(kRed);
  
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  f->ReleaseParameter(1);
  f->SetParLimits(1,1.86,1.87);
  f->ReleaseParameter(11);
  f->SetParLimits(11,-1.,1.);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
  
  TF1* background = new TF1(Form("background%d",count),"[0]+[1]*x+[2]*x*x+[3]*x*x*x");
  background->SetParameter(0,f->GetParameter(3));
  background->SetParameter(1,f->GetParameter(4));
  background->SetParameter(2,f->GetParameter(5));
  background->SetParameter(3,f->GetParameter(6));
  background->SetLineColor(4);
  background->SetRange(minhisto,maxhisto);
  background->SetLineStyle(2);
  
  TF1* mass = new TF1(Form("fmass_%.0f_%.0f",ptmin,ptmax),"[0]*([3]*([4]*Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))))");
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
  
  TF1* massSwap = new TF1(Form("fmassSwap_%.0f_%.0f",ptmin,ptmax),"[0]*(1-[2])*Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))");
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
  
  h->SetXTitle("m_{#piK} (GeV/c^{2})");
  h->SetYTitle("Entries / (5 MeV/c^{2})");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  h->SetStats(0);
  h->Draw("e");

  background->Draw("same");   
  mass->SetRange(minhisto,maxhisto);
  mass->Draw("same");
  massSwap->SetRange(minhisto,maxhisto);
  massSwap->Draw("same");
  f->Draw("same");
  
  Double_t yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
  Double_t yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);
  
  std::cout<<"YIELD="<<yield<<std::endl;

  TLegend* leg = new TLegend(0.65,0.58,0.82,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(h,"Data","pl");
  leg->AddEntry(f,"Fit","l");
  leg->AddEntry(mass,"D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
  leg->AddEntry(massSwap,"K-#pi swapped","f");
  leg->AddEntry(background,"Combinatorial","l");
  leg->Draw("same");

  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();

  TLatex* texCol = new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();

  TLatex* tex;

  tex = new TLatex(0.22,0.78,Form("%.1f < p_{T} < %.1f GeV/c",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  TString texper="%";
  tex = new TLatex(0.22,0.71,Form("Centrality %.0f-%.0f%s",centMin,centMax,texper.Data()));//0.2612903,0.8425793
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.22,0.83,"|y| < 1.0");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  c->SaveAs(Form("plotFits/DMass_cent_%.0f_%.0f_pt_%.0f_%.0f_dphi_%d.pdf",centMin,centMax,ptmin,ptmax,count));
  
  return mass;
}
