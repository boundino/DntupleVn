#!/bin/bash

CENTMIN=30
CENTMAX=50
TFEND=("inclusive" "v2_inpl" "v2_outpl" "v3_inpl" "v3_outpl")
#

DO_SAVEMASSHISTO=1
DO_BFEEDDOWNFRACTION=1
DO_PLOTFRACTIONS=1

##

INPUTDATA="/data/wangj/Data2015/Dntuple/DVn/Dntuple_crab_PbPb_HIMinimumBias1-7_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016_skim_cent${CENTMIN}${CENTMAX}.root"
INPUTMC="/data/HeavyFlavourRun2/MC2015/Dntuple/PbPb/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
INPUTMCNP="/data/HeavyFlavourRun2/MC2015/Dntuple/PbPb/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_nonprompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
OUTPUT_SAVEMASSHISTO="outfiles/MassHisto"

##

if [ $DO_SAVEMASSHISTO -eq 1 ]; then 
g++ saveMassHisto.C $(root-config --cflags --libs) -g -o saveMassHisto.exe 
./saveMassHisto.exe "$INPUTDATA" "$INPUTMC" "$INPUTMCNP" "$OUTPUT_SAVEMASSHISTO" "$CENTMIN" "$CENTMAX"
rm saveMassHisto.exe
fi

if [ $DO_BFEEDDOWNFRACTION -eq 1 ]; then 
g++ bFeedDownFraction.C $(root-config --cflags --libs) -g -o bFeedDownFraction.exe 
for tend in ${TFEND[@]}
do
    ./bFeedDownFraction.exe "$OUTPUT_SAVEMASSHISTO" "$tend" "$CENTMIN" "$CENTMAX"
done
rm bFeedDownFraction.exe
fi

if [ $DO_PLOTFRACTIONS -eq 1 ]; then 
g++ plotFractions.C $(root-config --cflags --libs) -g -o plotFractions.exe 
./plotFractions.exe "$CENTMIN" "$CENTMAX"
rm plotFractions.exe
fi