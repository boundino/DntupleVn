#!/bin/bash

CENTMIN=10
CENTMAX=30

DO_FONLL=1



##

INPUTDATA="/data/wangj/Data2015/Dntuple/DVn/Dntuple_crab_PbPb_HIMinimumBias1-7_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016_skim_cent${CENTMIN}${CENTMAX}.root"
INPUTMC="/data/HeavyFlavourRun2/MC2015/Dntuple/PbPb/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
INPUTMCNP="/data/HeavyFlavourRun2/MC2015/Dntuple/PbPb/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_nonprompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"

INPUTFONLLP="fonlls/FONLL_pp_promptDzero_5TeV_y1.dat"
INPUTFONLLNP="fonlls/FONLL_pp_nonpromptDzero_5TeV_y1.dat"
OUTPUTFONLLP="outfiles/FONLL_pp_promptDzero_5TeV_y1.root"
OUTPUTFONLLNP="outfiles/FONLL_pp_nonpromptDzero_5TeV_y1.root"

##

if [ $DO_FONLL -eq 1 ]; then 
g++ fonllDsigmadpt.C $(root-config --cflags --libs) -g -o fonllDsigmadpt.exe 
./fonllDsigmadpt.exe "$INPUTFONLLP" "$OUTPUTFONLLP"
./fonllDsigmadpt.exe "$INPUTFONLLNP" "$OUTPUTFONLLNP"
rm fonllDsigmadpt.exe
fi



