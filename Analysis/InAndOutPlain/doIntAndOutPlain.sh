#!/bin/bash

CENTMIN=30
CENTMAX=50
ISPROMPT=1

#

DO_SAVEMASSHISTO=0
DO_FITMASSHISTO=0
DO_CALVN=0
##

INPUTDATA="/data/wangj/Data2015/Dntuple/DVn/Dntuple_crab_PbPb_HIMinimumBias1-7_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016_skim_cent${CENTMIN}${CENTMAX}.root"
INPUTMC="/data/jisun/PbPb2015/Dntuple_crab_PbPb_Pythia8_prompt_D0pt0p0_AllPthat_Hydjet_MB_tkpt0p5eta1p5_Dy1p1_06092016_Pthatweight.root"
INPUT_BFEEDDOWN="../bFeedDown/outfilesResult/bFeedDownResult"
OUTPUT_SAVEMASSHISTO="outfiles/MassHisto"
OUTPUT_FITMASSHISTO="outfiles/YieldHisto"
OUTPUT_CALVN="outfiles/VnPtHisto"

##

if [ $DO_SAVEMASSHISTO -eq 1 ]; then 
g++ saveMassHisto.C $(root-config --cflags --libs) -g -o saveMassHisto.exe 
./saveMassHisto.exe "$INPUTDATA" "$INPUTMC" "$OUTPUT_SAVEMASSHISTO" "$CENTMIN" "$CENTMAX"
rm saveMassHisto.exe
fi

if [ $DO_FITMASSHISTO -eq 1 ]; then 
g++ fitMassHisto.C $(root-config --cflags --libs) -g -o fitMassHisto.exe 
./fitMassHisto.exe "$OUTPUT_SAVEMASSHISTO" "$OUTPUT_FITMASSHISTO" "$CENTMIN" "$CENTMAX"
rm fitMassHisto.exe
fi

if [ $DO_CALVN -eq 1 ]; then 
g++ calVn.C $(root-config --cflags --libs) -g -o calVn.exe 
./calVn.exe "$OUTPUT_FITMASSHISTO" "$OUTPUT_CALVN" "$CENTMIN" "$CENTMAX" "$ISPROMPT" "$INPUT_BFEEDDOWN"
rm calVn.exe
fi

###

DO_PLOTVN=1

if [ $DO_PLOTVN -eq 1 ]; then
g++ plotVn.C $(root-config --cflags --libs) -g -o plotVn.exe
./plotVn.exe "$OUTPUT_CALVN" "$CENTMIN" "$CENTMAX" 
rm plotVn.exe
fi