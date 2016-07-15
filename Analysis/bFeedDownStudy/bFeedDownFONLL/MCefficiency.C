using namespace std;
#include "../include/uti.h"
#include "../include/parameters.h"

void MCefficiency(TString inputmc, TString outputfile, TString tfend, TString selmcgen, TString cut, TString weight, Float_t centmin, Float_t centmax)
{
  Float_t hiBinMin,hiBinMax;
  hiBinMin = centmin*2;
  hiBinMax = centmax*2;

  selmcgen = Form("%s&&hiBin>=%f&&hiBin<%f",selmcgen.Data(),hiBinMin,hiBinMax);
  cut = Form("%s&&hiBin>=%f&&hiBin<%f",cut.Data(),hiBinMin,hiBinMax);

  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("ntDkpi");
  ntMC->AddFriend("ntHi");
  ntMC->AddFriend("ntSkim");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  ntGen->AddFriend("ntHi");
  ntGen->AddFriend("ntSkim");

  TH1D* hPtMC = new TH1D("hPtMC","",nPtBins,ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",nPtBins,ptBins);

  ntMC->Project("hPtMC","Dpt",TCut(weight)*(TCut(cut.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMC);
  ntGen->Project("hPtGen","Gpt",TCut(weight)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGen);
  TH1D* hEff = (TH1D*)hPtMC->Clone("hEff");
  hEff->Divide(hPtGen);

  TFile* fout=new TFile(Form("%s_cent_%.0f_%.0f_%s.root",outputfile.Data(),centmin,centmax,tfend.Data()),"recreate");
  fout->cd();
  hEff->Write();
  fout->Close();
}

int main(int argc, char *argv[])
{
  if(argc==9)
    {
      MCefficiency(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],atof(argv[7]),atof(argv[8]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
