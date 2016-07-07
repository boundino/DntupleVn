using namespace std;
#include "uti.h"
#include "setBranches.h"
#include "saveMassHisto.h"

int saveMassHisto(TString inputdata="", TString inputmc="",
                  TString outputfile="outfiles/MassHisto",
                  Float_t centmin=0., Float_t centmax=100.)
{
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TH1D* ah[nPtBins][4];
  TH1D* ahMCSignal[nPtBins];
  TH1D* ahMCSwapped[nPtBins];

  for(int i=0;i<nPtBins;i++)
    {
      for(int j=0;j<4;j++)
        {
          ah[i][j] = new TH1D(Form("h_%d_%s",i,tfend[j].Data()),"",nbinsmasshisto,minhisto,maxhisto);
          ah[i][j]->Sumw2();
        }
      ahMCSignal[i] = new TH1D(Form("hMCSignal_%d",i),"",nbinsmasshisto,minhisto,maxhisto);
      ahMCSignal[i]->Sumw2();
      ahMCSwapped[i] = new TH1D(Form("hMCSwapped_%d",i),"",nbinsmasshisto,minhisto,maxhisto);
      ahMCSwapped[i]->Sumw2();
    }

  TFile* inf = new TFile(inputdata.Data());
  TTree* nt = (TTree*)inf->Get("ntDkpi");
  TTree* ntHlt = (TTree*)inf->Get("ntHlt");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  TTree* ntSkim = (TTree*)inf->Get("ntSkim");

  SetRecoBranches(nt);
  SetHltBranches(ntHlt);
  SetHiBranches(ntHi,true);
  SetSkimBranches(ntSkim);

  cout<<endl;
  cout<<"--- Processing - DATA"<<endl;
  cout<<"--- Check the number of events for all trees"<<endl;
  cout<<nt->GetEntries()<<" "<<ntHlt->GetEntries()<<" "<<ntHi->GetEntries()<<" "<<ntSkim->GetEntries()<<endl;
  cout<<"--- Processing events"<<endl;
  Int_t nentries = nt->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      nt->GetEntry(i);
      ntHlt->GetEntry(i);
      ntHi->GetEntry(i);
      ntSkim->GetEntry(i);
      if(i%1000000==0) cout<<setiosflags(ios::left)<<"    "<<setw(8)<<i<<" / "<<nentries<<endl;
      if(!passtriggersel()||!passeventsel(centmin,centmax)) continue;
      for(int j=0;j<Dsize;j++)
        {
          Int_t iptbin = findptbin(j);
          if(iptbin<0) continue;
          if(!passcutsel(iptbin,j)) continue;
          Float_t* delphi = new Float_t[2];
          caldeltaphi(j,delphi);
          Int_t* ihistbin = new Int_t[2];
          if(delphi[0]>=0 && delphi[0]<(PI/4.)) ihistbin[0] = 0;//in-plain
          else ihistbin[0] = 1;//out-of-plain
          if(delphi[1]>=0 && delphi[1]<(PI/6.)) ihistbin[1] = 2;
          else ihistbin[1] = 3;
          ah[iptbin][ihistbin[0]]->Fill(Dmass[j]);
          ah[iptbin][ihistbin[1]]->Fill(Dmass[j]);
        }
    }

  //

  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("ntDkpi");
  TTree* ntMCHi = (TTree*)infMC->Get("ntHi");
  TTree* ntMCSkim = (TTree*)infMC->Get("ntSkim");

  SetRecoBranches(ntMC);
  SetHiBranches(ntMCHi,false);
  SetSkimBranches(ntMCSkim);

  cout<<endl;
  cout<<"--- Processing - MC: Prompt"<<endl;
  cout<<"--- Check the number of events for all trees"<<endl;
  cout<<ntMC->GetEntries()<<" "<<ntMCHi->GetEntries()<<" "<<ntMCSkim->GetEntries()<<endl;
  cout<<"--- Processing events"<<endl;
  nentries = ntMC->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      ntMC->GetEntry(i);
      ntMCHi->GetEntry(i);
      ntMCSkim->GetEntry(i);
      if(i%1000000==0) cout<<setiosflags(ios::left)<<"    "<<setw(8)<<i<<" / "<<nentries<<endl;
      if(!passeventsel(0,100)) continue;
      for(int j=0;j<Dsize;j++)
        {
          Int_t iptbin = findptbin(j);
          if(iptbin<0) continue;
          if(!passcutsel(iptbin,j)) continue;
          if(!passprompt(j)) continue;
          if(Dgen[j]==23333) ahMCSignal[iptbin]->Fill(Dmass[j]);
          else if(Dgen[j]==23344) ahMCSwapped[iptbin]->Fill(Dmass[j]);
        }
    }

  //
  cout<<endl;
  cout<<"--- Creating output files"<<endl;

  for(int i=0;i<nPtBins;i++)
    {
      TFile* outf = new TFile(Form("%s_cent_%.0f_%.0f_pt_%.0f_%.0f.root",outputfile.Data(),centmin,centmax,ptBins[i],ptBins[i+1]),"recreate");
      outf->cd();
      for(int j=0;j<4;j++)
        {      
          TH1D* h = (TH1D*)ah[i][j]->Clone(Form("h_%s",tfend[j].Data()));
          h->Write();
          delete h;
        }
      TH1D* hMCSignal = (TH1D*)ahMCSignal[i]->Clone("hMCSignal");;
      TH1D* hMCSwapped = (TH1D*)ahMCSwapped[i]->Clone("hMCSwapped");
      hMCSignal->Write();
      hMCSwapped->Write();
      delete hMCSwapped;
      delete hMCSignal;
      delete outf;
    }

  cout<<endl;
  return 0;
}

int main(int argc, char *argv[])
{
  if(argc==6)
    {
      saveMassHisto(argv[1],argv[2], argv[3], atof(argv[4]), atof(argv[5]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
