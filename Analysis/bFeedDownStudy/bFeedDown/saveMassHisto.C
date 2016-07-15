using namespace std;
#include "saveMassHisto.h"

int saveMassHisto(TString inputdata="", TString inputmc="", TString inputmcNP="",
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

  //
  initBins();
  //

  TH3D* hDataSigreg[3][2];
  TH3D* hDataSidbnd[3][2];
  TH3D* hDataPtMDca[3][2];

  for(int i=0;i<3;i++)
    {
      for(int j=0;j<2;j++)
        {
          if(i==2&&j==1) continue;
          hDataSigreg[i][j] = new TH3D(Form("hDataSigreg_%d_%d",i,j),"",nPtBins,ptBins,nDcaBins,dcaBins,nD0Bins,d0Bins);
          hDataSidbnd[i][j] = new TH3D(Form("hDataSidbnd_%d_%d",i,j),"",nPtBins,ptBins,nDcaBins,dcaBins,nD0Bins,d0Bins);          
          hDataPtMDca[i][j] = new TH3D(Form("hDataPtMDca_%d_%d",i,j),"",nPtBins,ptBins,nMassBins,massBins,nDcaBins,dcaBins);          
        }
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
          if(delphi[1]>=0 && delphi[1]<(PI/6.)) ihistbin[1] = 0;
          else ihistbin[1] = 1;
          hDataPtMDca[0][ihistbin[0]]->Fill(Dpt[j],Dmass[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]));
          hDataPtMDca[1][ihistbin[1]]->Fill(Dpt[j],Dmass[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]));
          hDataPtMDca[2][0]->Fill(Dpt[j],Dmass[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]));
          if(passsigreg(j))
            {
              hDataSigreg[0][ihistbin[0]]->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j]);              
              hDataSigreg[1][ihistbin[1]]->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j]);
              hDataSigreg[2][0]->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j]);
            }
          if(passsidbnd(j))
            {
              hDataSidbnd[0][ihistbin[0]]->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j]);              
              hDataSidbnd[1][ihistbin[1]]->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j]);
              hDataSidbnd[2][0]->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j]);
            }
        }
    }

  //
  TH3D* hMCSignalPr = new TH3D("hMCSignalPr","",nPtBins,ptBins,nDcaBins,dcaBins,nD0Bins,d0Bins);
  TH3D* hMCSignalNP = new TH3D("hMCSignalNP","",nPtBins,ptBins,nDcaBins,dcaBins,nD0Bins,d0Bins);
  TH3D* hMCPtMDcaSignal = new TH3D("hMCPtMDcaSignal","",nPtBins,ptBins,nMassBins,massBins,nDcaBins,dcaBins);
  TH3D* hMCPtMDcaSwapped = new TH3D("hMCPtMDcaSwapped","",nPtBins,ptBins,nMassBins,massBins,nDcaBins,dcaBins);

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
          if((Dgen[j]==23333||Dgen[j]==23344)&&passsigreg(j)) hMCSignalPr->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j],pthatweight);
          if(Dgen[j]==23333) hMCPtMDcaSignal->Fill(Dpt[j],Dmass[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]));
          if(Dgen[j]==23344) hMCPtMDcaSwapped->Fill(Dpt[j],Dmass[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]));
        }
    }

  //
  TFile* infMCNP = new TFile(inputmcNP.Data());
  TTree* ntMCNP = (TTree*)infMCNP->Get("ntDkpi");
  TTree* ntMCNPHi = (TTree*)infMCNP->Get("ntHi");
  TTree* ntMCNPSkim = (TTree*)infMCNP->Get("ntSkim");

  SetRecoBranches(ntMCNP);
  SetHiBranches(ntMCNPHi,false);
  SetSkimBranches(ntMCNPSkim);

  cout<<endl;
  cout<<"--- Processing - MC: Non-prompt"<<endl;
  cout<<"--- Check the number of events for all trees"<<endl;
  cout<<ntMCNP->GetEntries()<<" "<<ntMCNPHi->GetEntries()<<" "<<ntMCNPSkim->GetEntries()<<endl;
  cout<<"--- Processing events"<<endl;
  nentries = ntMCNP->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      ntMCNP->GetEntry(i);
      ntMCNPHi->GetEntry(i);
      ntMCNPSkim->GetEntry(i);
      if(i%1000000==0) cout<<setiosflags(ios::left)<<"    "<<setw(8)<<i<<" / "<<nentries<<endl;
      if(!passeventsel(0,100)) continue;
      for(int j=0;j<Dsize;j++)
        {
          Int_t iptbin = findptbin(j);
          if(iptbin<0) continue;
          if(!passcutsel(iptbin,j)) continue;
          if(passprompt(j)) continue;
          if((Dgen[j]==23333||Dgen[j]==23344)&&passsigreg(j)) hMCSignalNP->Fill(Dpt[j],DsvpvDistance[j]*TMath::Sin(Dalpha[j]),DsvpvDistance[j]/DsvpvDisErr[j],pthatweight);
        }
    }

  cout<<endl;
  cout<<"--- Creating output files"<<endl;

  for(int i=0;i<3;i++)
    {
      for(int j=0;j<2;j++)
        {
          if(i==2&&j==1) continue;

          TH3D* ohDataSigreg = (TH3D*)hDataSigreg[i][j]->Clone("hDataSigreg");
          TH3D* ohDataSidbnd = (TH3D*)hDataSidbnd[i][j]->Clone("hDataSidbnd");
          TH3D* ohDataPtMDca = (TH3D*)hDataPtMDca[i][j]->Clone("hDataPtMDca");
          TH3D* ohMCSignalPr = (TH3D*)hMCSignalPr->Clone("hMCSignalPr");
          TH3D* ohMCSignalNP = (TH3D*)hMCSignalNP->Clone("hMCSignalNP");
          TH3D* ohMCPtMDcaSignal = (TH3D*)hMCPtMDcaSignal->Clone("hMCPtMDcaSignal");
          TH3D* ohMCPtMDcaSwapped = (TH3D*)hMCPtMDcaSwapped->Clone("hMCPtMDcaSwapped");

          TFile* outf = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputfile.Data(),centmin,centmax,tfname[i][j].Data()),"recreate");
          outf->cd();
          ohDataSigreg->Write();
          ohDataSidbnd->Write();
          ohDataPtMDca->Write();
          ohMCSignalPr->Write();
          ohMCSignalNP->Write();
          ohMCPtMDcaSignal->Write();
          ohMCPtMDcaSwapped->Write();
          delete ohMCPtMDcaSwapped;
          delete ohMCPtMDcaSignal;
          delete ohMCSignalNP;
          delete ohMCSignalPr;
          delete ohDataPtMDca;
          delete ohDataSidbnd;
          delete ohDataSigreg;
          delete outf;
        }
    }

  cout<<endl;
  return 0;
}

int main(int argc, char *argv[])
{
  if(argc==7)
    {
      saveMassHisto(argv[1],argv[2], argv[3], argv[4],atof(argv[5]), atof(argv[6]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
