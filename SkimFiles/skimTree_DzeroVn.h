#ifndef _SKIM__SKIMTREEDZEROVN_H_
#define _SKIM__SKIMTREEDZEROVN_H_

#include "uti.h"

void SelectRecoBranches(TTree* nt)
{
  nt->SetBranchStatus("DtktkRes*",0);
  nt->SetBranchStatus("Dtrk3*",0);
  nt->SetBranchStatus("Dtrk4*",0);
  nt->SetBranchStatus("DRestrk*",0);
}

Int_t HLT_HIL1MinimumBiasHF2AND_part1_v1;
Int_t HLT_HIL1MinimumBiasHF2AND_part2_v1;
Int_t HLT_HIL1MinimumBiasHF2AND_part3_v1;
Int_t HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1;
Int_t HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1;
Int_t HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1;
Int_t HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1;
void SetHltBranchAddress(TTree* nt)
{
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part1_v1",&HLT_HIL1MinimumBiasHF2AND_part1_v1);
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part2_v1",&HLT_HIL1MinimumBiasHF2AND_part2_v1);
  nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part3_v1",&HLT_HIL1MinimumBiasHF2AND_part3_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1);
  nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1);
}

Int_t hiBin;
void SetHiBranchAddress(TTree* nt)
{
  nt->SetBranchAddress("hiBin",&hiBin);
}

bool hltselection(int PDno)
{
  if(PDno==1 && (HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1 && !(HLT_HIL1MinimumBiasHF2AND_part1_v1||HLT_HIL1MinimumBiasHF2AND_part2_v1||HLT_HIL1MinimumBiasHF2AND_part3_v1) && !(HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1||HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1||HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1))) return true;
  else if(PDno==2 && (HLT_HIL1MinimumBiasHF2AND_part1_v1)) return true;
  else if(PDno==3 && (HLT_HIL1MinimumBiasHF2AND_part2_v1)) return true;
  else if(PDno==4 && (HLT_HIL1MinimumBiasHF2AND_part3_v1)) return true;
  else if(PDno==5 && (HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1 && !(HLT_HIL1MinimumBiasHF2AND_part1_v1||HLT_HIL1MinimumBiasHF2AND_part2_v1||HLT_HIL1MinimumBiasHF2AND_part3_v1))) return true;
  else if(PDno==6 && (HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1 && !(HLT_HIL1MinimumBiasHF2AND_part1_v1||HLT_HIL1MinimumBiasHF2AND_part2_v1||HLT_HIL1MinimumBiasHF2AND_part3_v1))) return true;
  else if(PDno==7 && (HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1 && !(HLT_HIL1MinimumBiasHF2AND_part1_v1||HLT_HIL1MinimumBiasHF2AND_part2_v1||HLT_HIL1MinimumBiasHF2AND_part3_v1))) return true;
  else return false;
}

#endif
