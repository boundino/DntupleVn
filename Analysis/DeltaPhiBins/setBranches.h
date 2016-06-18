#ifndef _ANA_DELTAPHIBINS_SETBRANCHES_H_
#define _ANA_DELTAPHIBINS_SETBRANCHES_H_

#include "uti.h"
#define MAX_XB        20000
#define kMaxEvtPlanes 200

Int_t     Dsize;
Float_t   PVz;
Float_t   Dmass[MAX_XB];
Float_t   Dpt[MAX_XB];
Float_t   Deta[MAX_XB];
Float_t   Dphi[MAX_XB];
Float_t   Dy[MAX_XB];
Float_t   Dchi2ndf[MAX_XB];
Float_t   Dchi2cl[MAX_XB];
Float_t   Dalpha[MAX_XB];
Float_t   DsvpvDistance[MAX_XB];
Float_t   DsvpvDisErr[MAX_XB];
Float_t   DlxyBS[MAX_XB];
Float_t   DlxyBSErr[MAX_XB];
Float_t   Dtrk1Pt[MAX_XB];
Float_t   Dtrk2Pt[MAX_XB];
Float_t   Dtrk1Eta[MAX_XB];
Float_t   Dtrk2Eta[MAX_XB];
Float_t   Dtrk1Phi[MAX_XB];
Float_t   Dtrk2Phi[MAX_XB];
Float_t   Dtrk1PtErr[MAX_XB];
Float_t   Dtrk2PtErr[MAX_XB];
Float_t   Dtrk1PixelHit[MAX_XB];
Float_t   Dtrk2PixelHit[MAX_XB];
Float_t   Dtrk1StripHit[MAX_XB];
Float_t   Dtrk2StripHit[MAX_XB];
Float_t   Dtrk1nStripLayer[MAX_XB];
Float_t   Dtrk2nStripLayer[MAX_XB];
Float_t   Dtrk1nPixelLayer[MAX_XB];
Float_t   Dtrk2nPixelLayer[MAX_XB];
Float_t   Dtrk1Chi2ndf[MAX_XB];
Float_t   Dtrk2Chi2ndf[MAX_XB];
Bool_t    Dtrk1highPurity[MAX_XB];
Bool_t    Dtrk2highPurity[MAX_XB];
Float_t   Dgen[MAX_XB];

void SetRecoBranches(TTree* nt)
{
  nt->SetBranchAddress("Dsize",&Dsize);
  nt->SetBranchAddress("PVz",&PVz);
  nt->SetBranchAddress("Dmass",Dmass);
  nt->SetBranchAddress("Dpt",Dpt);
  nt->SetBranchAddress("Deta",Deta);
  nt->SetBranchAddress("Dphi",Dphi);
  nt->SetBranchAddress("Dy",Dy);
  nt->SetBranchAddress("Dchi2ndf",Dchi2ndf);
  nt->SetBranchAddress("Dchi2cl",Dchi2cl);
  nt->SetBranchAddress("Dalpha",Dalpha);
  nt->SetBranchAddress("DsvpvDistance",DsvpvDistance);
  nt->SetBranchAddress("DsvpvDisErr",DsvpvDisErr);
  nt->SetBranchAddress("DlxyBS",DlxyBS);
  nt->SetBranchAddress("DlxyBSErr",DlxyBSErr);
  nt->SetBranchAddress("Dtrk1Pt",Dtrk1Pt);
  nt->SetBranchAddress("Dtrk2Pt",Dtrk2Pt);
  nt->SetBranchAddress("Dtrk1Eta",Dtrk1Eta);
  nt->SetBranchAddress("Dtrk2Eta",Dtrk2Eta);
  nt->SetBranchAddress("Dtrk1Phi",Dtrk1Phi);
  nt->SetBranchAddress("Dtrk2Phi",Dtrk2Phi);
  nt->SetBranchAddress("Dtrk1PtErr",Dtrk1PtErr);
  nt->SetBranchAddress("Dtrk2PtErr",Dtrk2PtErr);
  nt->SetBranchAddress("Dtrk1PixelHit",Dtrk1PixelHit);
  nt->SetBranchAddress("Dtrk2PixelHit",Dtrk2PixelHit);
  nt->SetBranchAddress("Dtrk1StripHit",Dtrk1StripHit);
  nt->SetBranchAddress("Dtrk2StripHit",Dtrk2StripHit);
  nt->SetBranchAddress("Dtrk1nStripLayer",Dtrk1nStripLayer);
  nt->SetBranchAddress("Dtrk2nStripLayer",Dtrk2nStripLayer);
  nt->SetBranchAddress("Dtrk1nPixelLayer",Dtrk1nPixelLayer);
  nt->SetBranchAddress("Dtrk2nPixelLayer",Dtrk2nPixelLayer);
  nt->SetBranchAddress("Dtrk1Chi2ndf",Dtrk1Chi2ndf);
  nt->SetBranchAddress("Dtrk2Chi2ndf",Dtrk2Chi2ndf);
  nt->SetBranchAddress("Dtrk1highPurity",Dtrk1highPurity);
  nt->SetBranchAddress("Dtrk2highPurity",Dtrk2highPurity);
  nt->SetBranchAddress("Dgen",Dgen);
}

Int_t     HLT_HIL1MinimumBiasHF1AND_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part1_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part2_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part3_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1;
Int_t     HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1;
void SetHltBranches(TTree* nt)
{
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&HLT_HIL1MinimumBiasHF1AND_v1);
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part1_v1",&HLT_HIL1MinimumBiasHF2AND_part1_v1);
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part2_v1",&HLT_HIL1MinimumBiasHF2AND_part2_v1);
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part3_v1",&HLT_HIL1MinimumBiasHF2AND_part3_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1",&HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1);
}

Int_t     hiBin;
Int_t     hiNevtPlane;
Float_t   hiEvtPlanes[kMaxEvtPlanes];
Float_t   hiEvtPlanesqx[kMaxEvtPlanes];
Float_t   hiEvtPlanesqy[kMaxEvtPlanes];
Float_t   pthatweight;
void SetHiBranches(TTree* nt, Bool_t isData)
{
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("hiNevtPlane",&hiNevtPlane);
  nt->SetBranchAddress("hiEvtPlanes",hiEvtPlanes);
  nt->SetBranchAddress("hiEvtPlanesqx",hiEvtPlanesqx);
  nt->SetBranchAddress("hiEvtPlanesqy",hiEvtPlanesqy);
  if(!isData) nt->SetBranchAddress("pthatweight",&pthatweight);
}

Int_t     pcollisionEventSelection;
Int_t     pprimaryVertexFilter;
Int_t     phfCoincFilter3;
Int_t     pclusterCompatibilityFilter;
void SetSkimBranches(TTree* nt)
{
  nt->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  nt->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  nt->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
  nt->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
}

#endif
