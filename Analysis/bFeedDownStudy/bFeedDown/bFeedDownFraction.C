using namespace std;
#include "../include/parameters.h"
#include "bFeedDownFraction.h"

void bFeedDownFraction(TString inputfile="", TString tfend="", Float_t centmin=0., Float_t centmax=100.)
{
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetTitleX(.0f);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* c4 = new TCanvas("c4","",800,600);
  c4->Divide(2,2);
  TCanvas* c2 = new TCanvas("c2","",400,600);
  c2->Divide(1,2);
  TCanvas* c1 = new TCanvas("c1","",600,400);
  TCanvas* c15 = new TCanvas("c15","",810,1000);
  c15->Divide(3,5);

  TFile* inf = new TFile(Form("%s_cent_%.0f_%.0f_%s.root",inputfile.Data(),centmin,centmax,tfend.Data()));
  TH3D* hData = (TH3D*)inf->Get("hDataSigreg");
  TH3D* hSideband = (TH3D*)inf->Get("hDataSidbnd");
  TH3D* hPtMD0Dca = (TH3D*)inf->Get("hDataPtMDca");
  TH3D* hMCPSignal = (TH3D*)inf->Get("hMCSignalPr");
  TH3D* hMCNPSignal = (TH3D*)inf->Get("hMCSignalNP");
  TH3D* hPtMD0DcaMCPSignal =(TH3D*)inf->Get("hMCPtMDcaSignal");
  TH3D* hPtMD0DcaMCPSwapped = (TH3D*)inf->Get("hMCPtMDcaSwapped");

  hData->Sumw2();
  hSideband->Sumw2();
  hPtMD0Dca->Sumw2();
  hMCPSignal->Sumw2();
  hMCNPSignal->Sumw2();
  hPtMD0DcaMCPSignal->Sumw2();
  hPtMD0DcaMCPSwapped->Sumw2();

  float invMSideBandSwitchPt = 15;

  float pts[nPtBins];
  float ptErrors[nPtBins];
  float promptFraction[nPtBins];
  float promptFractionError[nPtBins];
  float promptFractionErrorDataOnly[nPtBins];
  float promptFraction_DCA[nPtBins];
  float promptFractionError_DCA[nPtBins];
  float totalYield[nPtBins];
  float totalYieldError[nPtBins];
  float bToDYield[nPtBins];
  float bToDYieldError[nPtBins];
  float promptDYield[nPtBins];
  float promptDYieldError[nPtBins];

  initBins();

  TLatex* texCms = new TLatex(0.14,0.95, "#scale[1.25]{CMS} #bf{#it{Preliminary}}");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.06);
  texCms->SetTextFont(42);
  TLatex* texCol = new TLatex(0.96,0.95, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.06);
  texCol->SetTextFont(42);

  for(int i=1;i<nPtBins+1;i++)
    {
      pts[i-1] = 0.5*(ptBins[i-1]+ptBins[i]);
      ptErrors[i-1] = 0.5*(ptBins[i]-ptBins[i-1]);
      float ptLow = ptBins[i-1];
      float ptHigh = ptBins[i];
      cout<<endl;
      cout<<"  --- Processing pT range: "<<ptLow<<" - "<<ptHigh<<endl;

      TLatex* texPtY = new TLatex(0.94,0.82,Form("%.1f < p_{T} < %.1f GeV/c      |y| < 1.0",ptLow,ptHigh));
      texPtY->SetNDC();
      texPtY->SetTextAlign(32);
      texPtY->SetTextFont(42);
      texPtY->SetTextSize(0.06);
      texPtY->SetLineWidth(2);

      TLatex* texPt = new TLatex(0.18,0.82,Form("%.1f < p_{T} < %.1f GeV/c",ptLow,ptHigh));
      texPt->SetNDC();
      texPt->SetTextFont(42);
      texPt->SetTextSize(0.06);
      texPt->SetLineWidth(2);

      TLatex* texY = new TLatex(0.18,0.74,Form("|y| < 1.0"));
      texY->SetNDC();
      texY->SetTextFont(42);
      texY->SetTextSize(0.06);
      texY->SetLineWidth(2);

      c2->cd(1);

      hPtMD0Dca->GetZaxis()->SetRange(1,100);///////////////////////////////
      hPtMD0Dca->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);
      hPtMD0DcaMCPSignal->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);
      hPtMD0DcaMCPSwapped->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);
      TH1D* hMData = (TH1D*)hPtMD0Dca->Project3D("y")->Clone(Form("hM_%1.1f_%1.1f", ptLow, ptHigh));
      TH1D* hMMCSignal = (TH1D*)hPtMD0DcaMCPSignal->Project3D("y");
      TH1D* hMMCSwapped = (TH1D*)hPtMD0DcaMCPSwapped->Project3D("y");

      setColorTitleLabel(hMData, kBlack);
      TF1* fMass = fitMass(hMData, hMMCSignal, hMMCSwapped);

      texCms->Draw();
      texCol->Draw();
      texPt->Draw();
      texY->Draw();

      TF1* background = new TF1("fBackground","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
      background->SetParameter(0,fMass->GetParameter(3));
      background->SetParameter(1,fMass->GetParameter(4));
      background->SetParameter(2,fMass->GetParameter(5));
      background->SetParameter(3,fMass->GetParameter(6));

      float massD = 1.8649;
      float massSignal1 = massD-0.025;
      float massSignal2 = massD+0.025;
      float massSideBand1 = massD-0.1;
      float massSideBand2 = massD-0.075;
      float massSideBand3 = massD+0.075;
      float massSideBand4 = massD+0.1;

      float scaleSideBand = background->Integral(massSignal1, massSignal2)/(background->Integral(massSideBand1, massSideBand2)+background->Integral(massSideBand3, massSideBand4));
      cout<<"scaleSideBand: "<<scaleSideBand<<endl;

      TLatex* texScale = new TLatex(0.18,0.66,Form("side band scale: %1.3f", scaleSideBand));
      texScale->SetNDC();
      texScale->SetTextFont(42);
      texScale->SetTextSize(0.06);
      texScale->SetLineWidth(2);
      texScale->Draw();

      TLine* lineSignal1 = new TLine(massSignal1, 0, massSignal1, background->Eval(massSignal1));  lineSignal1->SetLineWidth(2);  lineSignal1->SetLineStyle(5);  lineSignal1->SetLineColor(14);
      TLine* lineSignal2 = new TLine(massSignal2, 0, massSignal2, background->Eval(massSignal2));  lineSignal2->SetLineWidth(2);  lineSignal2->SetLineStyle(5);  lineSignal2->SetLineColor(14);
      TLine* lineSideBand1 = new TLine(massSideBand1, 0, massSideBand1, background->Eval(massSideBand1));  lineSideBand1->SetLineWidth(2);  lineSideBand1->SetLineStyle(5);  lineSideBand1->SetLineColor(14);
      TLine* lineSideBand2 = new TLine(massSideBand2, 0, massSideBand2, background->Eval(massSideBand2));  lineSideBand2->SetLineWidth(2);  lineSideBand2->SetLineStyle(5);  lineSideBand2->SetLineColor(14);
      TLine* lineSideBand3 = new TLine(massSideBand3, 0, massSideBand3, background->Eval(massSideBand3));  lineSideBand3->SetLineWidth(2);  lineSideBand3->SetLineStyle(5);  lineSideBand3->SetLineColor(14);
      TLine* lineSideBand4 = new TLine(massSideBand4, 0, massSideBand4, background->Eval(massSideBand4));  lineSideBand4->SetLineWidth(2);  lineSideBand4->SetLineStyle(5);  lineSideBand4->SetLineColor(14);
      lineSignal1->Draw();
      lineSignal2->Draw();
      lineSideBand1->Draw();
      lineSideBand2->Draw();
      lineSideBand3->Draw();
      lineSideBand4->Draw();

      c2->cd(2);
      gPad->SetLogy();

      hData->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);
      hSideband->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);
      hMCPSignal->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);
      hMCNPSignal->GetXaxis()->SetRangeUser(ptLow+0.001,ptHigh-0.001);

      TH1D* hD0DcaData0 = (TH1D*)hData->Project3D("y")->Clone("hD0DcaData0");
      TH1D* hD0DcaSideband = (TH1D*)hSideband->Project3D("y")->Clone("hD0DcaSideband");
      TH1D* hD0DcaMCPSignal0 = (TH1D*)hMCPSignal->Project3D("y")->Clone("hD0DcaMCPSignal0");
      TH1D* hD0DcaMCNPSignal0 = (TH1D*)hMCNPSignal->Project3D("y")->Clone("hD0DcaMCNPSignal0");

      float integralRawYieldMCP = hD0DcaMCPSignal0->Integral();
      float integralRawYieldMCNP = hD0DcaMCNPSignal0->Integral();
      cout<<"integralRawYieldMCP: "<<integralRawYieldMCP<<endl;
      cout<<"integralRawYieldMCNP: "<<integralRawYieldMCNP<<endl;

      hD0DcaMCPSignal = hD0DcaMCPSignal0;
      hD0DcaMCNPSignal = hD0DcaMCNPSignal0;

      divideBinWidth(hD0DcaData0);
      divideBinWidth(hD0DcaSideband);
      setColorTitleLabel(hD0DcaData0, 1);
      hD0DcaData0->GetXaxis()->SetRangeUser(0,0.07);
      hD0DcaData0->GetYaxis()->SetTitle("counts per cm");

      TH1D* hD0DcaSideband0 = (TH1D*)hD0DcaSideband->Clone("hD0DcaSideband0");
      hD0DcaSideband->Scale(scaleSideBand);

      TH1D* hD0DcaDataSubSideBand = (TH1D*)hD0DcaData0->Clone("hD0DcaDataSubSideBand");
      hD0DcaDataSubSideBand->Add(hD0DcaSideband,-1);

      hD0DcaData0->SetMarkerSize(0.6);
      hD0DcaData0->SetTitle(";D^{0} DCA (cm);counts per cm");
      hD0DcaData0->Draw();
      hD0DcaSideband->Draw("hsame");
      hD0DcaSideband0->SetLineStyle(2);
      hD0DcaSideband0->Draw("hsame");

      TLegend* leg1 = new TLegend(0.44,0.6,0.90,0.76,NULL,"brNDC");
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.06);
      leg1->SetTextFont(42);
      leg1->SetFillStyle(0);
      leg1->AddEntry(hD0DcaData0,"D^{0} candidate","pl");
      leg1->AddEntry(hD0DcaSideband,"side band","l");
      leg1->AddEntry(hD0DcaSideband0,"side band unscaled","l");
      leg1->Draw("same");

      texCms->Draw();
      texCol->Draw();
      texPtY->Draw();

      c2->SaveAs(Form("plots/PbPb_cent_%.0f_%.0f_%s_pt_%.1f_%.1f_sideBand.pdf",centmin,centmax,tfend.Data(),ptLow,ptHigh));

      c15->cd(1);

      fitMass(hMData, hMMCSignal, hMMCSwapped);

      texPt->Draw();
      texY->Draw();

      TH1D* hD0DcaDataFit = new TH1D("hD0DcaDataFit", ";D^{0} DCA (cm);dN / d(D^{0} DCA) (cm^{-1})", nDcaBins, dcaBins);

      for(int j=1; j<=14; j++)
        {
          c15->cd(j+1);
          hPtMD0Dca->GetZaxis()->SetRange(j,j);
          float D0DcaLow = hPtMD0Dca->GetZaxis()->GetBinLowEdge(j);
          float D0DcaHigh = hPtMD0Dca->GetZaxis()->GetBinUpEdge(j);
          TH1D* hMData_D0Dca = (TH1D*)hPtMD0Dca->Project3D("y")->Clone(Form("hM_pt_%1.1f_%1.1f_D0Dca_%1.4f_%1.4f", ptLow, ptHigh, D0DcaLow, D0DcaHigh));
          setColorTitleLabel(hMData_D0Dca,kBlack);
          fMass = fitMass(hMData_D0Dca, hMMCSignal, hMMCSwapped);

          float yield = fMass->GetParameter(0)*fMass->GetParameter(7)/hMData_D0Dca->GetBinWidth(1);
          float yieldError = fMass->GetParError(0)*fMass->GetParameter(7)/hMData_D0Dca->GetBinWidth(1);

          hD0DcaDataFit->SetBinContent(j, yield);
          hD0DcaDataFit->SetBinError(j, yieldError);

          TLatex* texD0Dca = new TLatex(0.18,0.82,Form("D^{0} DCA: %1.4f - %1.4f",D0DcaLow,D0DcaHigh));
          texD0Dca->SetNDC();
          texD0Dca->SetTextFont(42);
          texD0Dca->SetTextSize(0.06);
          texD0Dca->SetLineWidth(2);
          texD0Dca->Draw();

          TLatex* texYield = new TLatex(0.18,0.74,Form("D^{0} yield: %1.0f #pm %1.0f",yield,yieldError));
          texYield->SetNDC();
          texYield->SetTextFont(42);
          texYield->SetTextSize(0.06);
          texYield->SetLineWidth(2);
          texYield->Draw();

        }

      c15->SaveAs(Form("plots/PbPb_cent_%.0f_%.0f_%s_pt_%.1f_%.1f_invMassFit.pdf",centmin,centmax,tfend.Data(),ptLow,ptHigh));

      divideBinWidth(hD0DcaDataFit);

      c4->cd(1);
      gPad->SetLogy();

      normalize(hD0DcaMCPSignal);
      setColorTitleLabel(hD0DcaMCPSignal, 2);
      hD0DcaMCPSignal->GetXaxis()->SetRangeUser(0,0.07);

      normalize(hD0DcaMCNPSignal);
      setColorTitleLabel(hD0DcaMCNPSignal, 4);
      hD0DcaMCNPSignal->GetXaxis()->SetRangeUser(0,0.07);
      hD0DcaMCNPSignal->GetYaxis()->SetTitle("dN / d(D^{0} DCA) (cm^{-1})");
      hD0DcaMCNPSignal->GetXaxis()->SetTitle("D^{0} DCA (cm)");
      hD0DcaMCNPSignal->SetMaximum(hD0DcaMCPSignal->GetMaximum()*3.);

      hD0DcaMCNPSignal->Draw("");
      hD0DcaMCPSignal->Draw("same");

      TLegend* leg2 = new TLegend(0.52,0.72,0.90,0.88,NULL,"brNDC");
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.06);
      leg2->SetTextFont(42);
      leg2->SetFillStyle(0);
      leg2->AddEntry(hD0DcaMCPSignal,"MC Prompt D^{0}","pl");
      leg2->AddEntry(hD0DcaMCNPSignal,"MC Non-prompt D^{0}","pl");
      leg2->Draw("same");

      c4->cd(2);
      gPad->SetLogy();

      TH1D* hD0DcaData = hD0DcaDataFit;
      if(pts[i-1]>invMSideBandSwitchPt) hD0DcaData = hD0DcaDataSubSideBand;

      setColorTitleLabel(hD0DcaData, 1);

      double integralTotalYield = hD0DcaData->Integral(1,hD0DcaData->GetXaxis()->GetNbins(),"width");
      cout<<"integralTotalYield: "<<integralTotalYield<<endl;

      TF1* fMix = new TF1("fMix",&funMix, 0., 0.5, 2);
      fMix->SetParameters(integralTotalYield,0.9);
      fMix->SetParLimits(0,0,2*integralTotalYield);
      fMix->SetParLimits(1,0,1);

      fMix->SetLineColor(2);
      fMix->SetFillColor(kRed-9);
      fMix->SetFillStyle(1001);

      cout<<"   -- Fitting DCA"<<endl;

      float fitRangeL = 0;
      float fitRangeH = 0.07;

      hD0DcaData->GetXaxis()->SetRangeUser(0,0.07);
      hD0DcaData->Draw();
      int fitStatus = 1;
      TFitResultPtr fitResult;
      double fitPrecision = 1.e-6;
      while(fitStatus)
        {
          TFitter::SetPrecision(fitPrecision);
          fMix->SetParameters(integralTotalYield,0.9);
          fMix->SetParError(0,0.1*integralTotalYield);
          fMix->SetParError(1,0.1);
          fitResult = hD0DcaData->Fit("fMix","E SNQ0", "", fitRangeL, fitRangeH);
          fitStatus = fitResult->Status();
          cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<endl;
          if(fitStatus)
            fitPrecision *= 10;
        }
      fMix->SetParameters(integralTotalYield,0.9);
      fMix->SetParError(0,0.1*integralTotalYield);
      fMix->SetParError(1,0.1);
      fMix->SetNpx(10000);
      fitResult = hD0DcaData->Fit("fMix","E S0", "", fitRangeL, fitRangeH);
      hD0DcaData->GetFunction("fMix")->Draw("flsame");
      fitStatus = fitResult->Status();
      cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<endl;

      TF1* fNP = new TF1("fNP",&funNonPrompt, 0., 0.5, 2);
      fNP->SetParameters(fMix->GetParameter(0),fMix->GetParameter(1));
      fNP->SetRange(fitRangeL,fitRangeH);
      fNP->SetLineColor(4);
      fNP->SetFillStyle(1001);
      fNP->SetFillColor(kBlue-9);
      fNP->SetNpx(10000);
      fNP->Draw("same");

      hD0DcaData->Draw("same");

      //      cout<<"NP integral fraction: "<<fNP->Integral(fitRangeL,fitRangeH,1.e-7)/fMix->Integral(fitRangeL,fitRangeH,1.e-7)<<endl;
      //cout<<"NP integral fraction: "<<fNP->Integral(fitRangeL,fitRangeH,fNP->GetParameters(),1.e-7)/fMix->Integral(fitRangeL,fitRangeH,fNP->GetParameters(),1.e-7)<<endl;
      cout<<"prompt fraction (integral): "<<1-fNP->Integral(fitRangeL,fitRangeH)/fMix->Integral(fitRangeL,fitRangeH)<<endl;
      cout<<"prompt fraction (real data statistic error only): "<<fMix->GetParameter(1)<<" +- "<<fMix->GetParError(1)<<endl;
      cout<<"chi2 / NDF: "<<fitResult->Chi2()<<" / "<<fitResult->Ndf()<<endl;
      cout<<"total yield: "<<integralTotalYield<<" (integral) vs. "<<fMix->GetParameter(0)<<" +- "<<fMix->GetParError(0)<<" (fit)"<<endl;

      promptFraction[i-1] = fMix->GetParameter(1);
      promptFractionErrorDataOnly[i-1] = fMix->GetParError(1);
      promptFraction_DCA[i-1] = 1-fNP->Integral(0,0.008)/fMix->Integral(0,0.008);
      totalYield[i-1] = fMix->GetParameter(0);
      totalYieldError[i-1] = fMix->GetParError(0);

      texCms->Draw();
      texCol->Draw();
      texPtY->Draw();

      TLatex* texRatio = new TLatex(0.47,0.73,Form("Prompt frac. = %.1f #pm %.1f %%",100*fMix->GetParameter(1),100*fMix->GetParError(1)));
      texRatio->SetNDC();
      texRatio->SetTextFont(42);
      texRatio->SetTextSize(0.06);
      texRatio->SetLineWidth(2);
      texRatio->Draw();

      TLegend* leg4 = new TLegend(0.56,0.46,0.90,0.7);
      leg4->SetBorderSize(0);
      leg4->SetTextSize(0.06);
      leg4->SetTextFont(42);
      leg4->SetFillStyle(0);
      leg4->AddEntry(hD0DcaData,"Data","pl");
      leg4->AddEntry(fMix,"Prompt D^{0}","f");
      leg4->AddEntry(fNP,"Non-Prompt D^{0}","f");
      leg4->Draw("same");

      //smear MC smaple with the error, to simulate the MC statistic error effect.
      c4->cd(3);
      gPad->SetLogy();

      hD0DcaMCPSignal = (TH1D*)hD0DcaMCPSignal0->Clone("hMCPSignal");
      hD0DcaMCNPSignal = (TH1D*)hD0DcaMCNPSignal0->Clone("hMCNPSignal");

      TH1D* hPromptRatio = new TH1D("hPromptRatio", "hPromptRatio", 100, 0, 1);
      setColorTitleLabel(hPromptRatio, 1);

      int nSmear = 1000;

      for(int j=0; j<nSmear; j++)
        {
          RandomSmear(hD0DcaMCPSignal0, hD0DcaMCPSignal);
          RandomSmear(hD0DcaMCNPSignal0, hD0DcaMCNPSignal);

          hD0DcaData->Fit("fMix","E QN0");

          hPromptRatio->Fill(fMix->GetParameter(1));
        }

      hPromptRatio->GetXaxis()->SetTitle("prompt ratio");
      hPromptRatio->GetYaxis()->SetTitle("counts");
      hPromptRatio->SetMarkerStyle(20);
      hPromptRatio->SetStats(0);
      hPromptRatio->Draw("e");
      hPromptRatio->Fit("gaus");

      TLatex* texGaussMean = new TLatex(0.27,0.83,Form("#mu: %.3f",hPromptRatio->GetFunction("gaus")->GetParameter(1)));
      texGaussMean->SetNDC();
      texGaussMean->SetTextFont(42);
      texGaussMean->SetTextSize(0.06);
      texGaussMean->SetLineWidth(2);
      texGaussMean->Draw();

      TLatex* texGaussSigma = new TLatex(0.27,0.73,Form("#sigma: %.3f",hPromptRatio->GetFunction("gaus")->GetParameter(2)));
      texGaussSigma->SetNDC();
      texGaussSigma->SetTextFont(42);
      texGaussSigma->SetTextSize(0.06);
      texGaussSigma->SetLineWidth(2);
      texGaussSigma->Draw();

      float promptFractionErrorMc = hPromptRatio->GetFunction("gaus")->GetParameter(2);
      promptFractionError[i-1] = sqrt(pow(promptFractionErrorDataOnly[i-1],2)+pow(promptFractionErrorMc,2));
      promptFractionError_DCA[i-1] = promptFraction_DCA[i-1]*promptFractionError[i-1]/promptFraction[i-1];
      cout<<"prompt fraction: "<<promptFraction[i-1]<<" +- "<<promptFractionError[i-1]<<" (+- "<<promptFractionErrorDataOnly[i-1]<<" +- "<<promptFractionErrorMc<<" )"<<endl;

      bToDYield[i-1] = totalYield[i-1]*(1-promptFraction[i-1]);
      bToDYieldError[i-1] = bToDYield[i-1]*sqrt(pow(promptFractionError[i-1]/(1-promptFraction[i-1]),2) + pow(totalYieldError[i-1]/totalYield[i-1],2));
      cout<<"B to D yield: "<<bToDYield[i-1]<<" +- "<<bToDYieldError[i-1]<<endl;

      promptDYield[i-1] = totalYield[i-1]*promptFraction[i-1];
      promptDYieldError[i-1] = promptDYield[i-1]*sqrt(pow(promptFractionError[i-1]/promptFraction[i-1],2) + pow(totalYieldError[i-1]/totalYield[i-1],2));
      cout<<"prompt D yield: "<<promptDYield[i-1]<<" +- "<<promptDYieldError[i-1]<<endl;

      float promptMCScale = promptDYield[i-1]/integralRawYieldMCP;
      float nonPromptMCScale = bToDYield[i-1]/integralRawYieldMCNP;

      cout<<"promptMCScale: "<<promptMCScale<<endl;
      cout<<"nonPromptMCScale: "<<nonPromptMCScale<<endl;

      //restore original unsmeared histograms before saving plots
      delete hD0DcaMCPSignal;
      delete hD0DcaMCNPSignal;
      hD0DcaMCPSignal = hD0DcaMCPSignal0;
      hD0DcaMCNPSignal = hD0DcaMCNPSignal0;
      hD0DcaData->Fit("fMix","E QN0");

      c4->cd(4);

      TH1D* hD0DcaDataOverFit = (TH1D*)hD0DcaData->Clone("hD0DcaDataOverFit");
      hD0DcaDataOverFit->Divide(fMix);
      hD0DcaDataOverFit->GetYaxis()->SetTitle("data / fit");
      hD0DcaDataOverFit->GetYaxis()->SetRangeUser(0,5);
      hD0DcaDataOverFit->GetXaxis()->SetRangeUser(0,0.07);
      setColorTitleLabel(hD0DcaDataOverFit, 1);
      hD0DcaDataOverFit->Draw("e");

      TF1* fLine1 = new TF1("fLine1", "1", 0,1);
      fLine1->Draw("same");
      hD0DcaDataOverFit->Draw("esame");

      c4->SaveAs(Form("plots/PbPb_cent_%.0f_%.0f_%s_pt_%.1f_%.1f_fit.pdf",centmin,centmax,tfend.Data(),ptLow,ptHigh));

      delete hD0DcaMCPSignal;
      delete hD0DcaMCNPSignal;

    }
  c1->cd();

  TH1D* hStupidJie = new TH1D("hStupidJie", "", 40, 0, 40);
  hStupidJie->GetYaxis()->SetRangeUser(0,1);
  hStupidJie->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hStupidJie->GetYaxis()->SetTitle("prompt fraction");
  hStupidJie->SetStats(0);
  hStupidJie->Draw();
  TGraphErrors* grFraction = new TGraphErrors(nPtBins, pts, promptFraction, ptErrors, promptFractionError);
  grFraction->SetName("grPromptFraction");
  grFraction->SetMarkerStyle(20);
  grFraction->Draw("psame");

  TGraphErrors* grFraction2 = new TGraphErrors(nPtBins, pts, promptFraction, ptErrors, promptFractionErrorDataOnly);
  grFraction2->SetName("grPromptFractionErrorFromRealDataOnly");
  grFraction2->SetMarkerStyle(20);
  grFraction2->SetMarkerColor(4);
  grFraction2->SetLineColor(4);
  grFraction2->Draw("psame");

  TGraphErrors* grFraction_DCA = new TGraphErrors(nPtBins, pts, promptFraction_DCA, ptErrors, promptFractionError_DCA);
  grFraction_DCA->SetName("grPromptFraction_DCA");
  grFraction_DCA->SetMarkerStyle(20);
  grFraction_DCA->Draw("psame");

  TLegend* leg = new TLegend(0.2, 0.3, 0.9, 0.5);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(grFraction, "statistic error from both real data and MC", "l");
  leg->AddEntry(grFraction2, "statistic error from real data only", "l");
  leg->Draw();

  c1->SaveAs(Form("plotsResult/promptFraction_cent_%.0f_%.0f_%s.pdf",centmin,centmax,tfend.Data()));

  c1->SetLogy();

  TH1D* hBtoDRawYield = new TH1D("hBtoDRawYield", ";p_{T} (GeV/c);counts per GeV/c", nPtBins, ptBins);
  for(int i=1; i<=nPtBins; i++)
    {
      if(bToDYield[i-1] <= 0) continue;
      hBtoDRawYield->SetBinContent(i, bToDYield[i-1]);
      hBtoDRawYield->SetBinError(i, bToDYieldError[i-1]);
    }
  divideBinWidth(hBtoDRawYield);
  setColorTitleLabel(hBtoDRawYield, 1);
  c1->SetBottomMargin(0.14);
  hBtoDRawYield->Draw("p");

  c1->SaveAs(Form("plotsResult/BtoD_cent_%.0f_%.0f_%s.pdf",centmin,centmax,tfend.Data()));

  TH1D* hPromptDRawYield = new TH1D("hPromptDRawYield", ";p_{T} (GeV/c);counts per GeV/c", nPtBins, ptBins);
  for(int i=1; i<=nPtBins; i++)
    {
      if(promptDYield[i-1] <= 0) continue;
      hPromptDRawYield->SetBinContent(i, promptDYield[i-1]);
      hPromptDRawYield->SetBinError(i, promptDYieldError[i-1]);
    }
  divideBinWidth(hPromptDRawYield);
  setColorTitleLabel(hPromptDRawYield, 1);
  c1->SetBottomMargin(0.14);
  hPromptDRawYield->Draw("p");

  c1->SaveAs(Form("plotsResult/promptD_cent_%.0f_%.0f_%s.pdf",centmin,centmax,tfend.Data()));

  TFile* fOut = new TFile(Form("outfilesResult/bFeedDownResult_cent_%.0f_%.0f_%s.root",centmin,centmax,tfend.Data()), "recreate");
  grFraction->Write();
  grFraction2->Write();
  grFraction_DCA->Write();
  hBtoDRawYield->Write();
  hPromptDRawYield->Write();
  fOut->Close();
}

int main(int argc, char *argv[])
{
  if(argc==5)
    {
      bFeedDownFraction(argv[1],argv[2], atof(argv[3]), atof(argv[4]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
