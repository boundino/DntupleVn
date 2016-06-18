#!/bin/bash

DO_SAVEMASSHISTO=0
DO_FITMASSHISTO=0
DO_FITDELTAPHI=1

INPUTDATA="/data/wangj/Data2015/Dntuple/DVn/Dntuple_crab_PbPb_HIMinimumBias1-7_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016_skim_cent3050.root"
INPUTMC="/data/jisun/PbPb2015/Dntuple_crab_PbPb_Pythia8_prompt_D0pt0p0_AllPthat_Hydjet_MB_tkpt0p5eta1p5_Dy1p1_06092016_Pthatweight.root"
OUTPUT_SAVEMASSHISTO="outfiles/MassHisto"
OUTPUT_FITMASSHISTO="outfiles/DeltaPhiHisto"

if [ $DO_SAVEMASSHISTO -eq 1 ]; then 
g++ saveMassHisto.C $(root-config --cflags --libs) -g -o saveMassHisto.exe 
./saveMassHisto.exe "$INPUTDATA" "$INPUTMC" "$OUTPUT_SAVEMASSHISTO" 30 50
rm saveMassHisto.exe
fi

if [ $DO_FITMASSHISTO -eq 1 ]; then 
g++ fitMassHisto.C $(root-config --cflags --libs) -g -o fitMassHisto.exe 
./fitMassHisto.exe "$OUTPUT_SAVEMASSHISTO" "$OUTPUT_FITMASSHISTO" 30 50
rm fitMassHisto.exe
fi

if [ $DO_FITDELTAPHI -eq 1 ]; then 
g++ fitDeltaPhi.C $(root-config --cflags --libs) -g -o fitDeltaPhi.exe 
./fitDeltaPhi.exe "$OUTPUT_FITMASSHISTO" 30 50
rm fitDeltaPhi.exe
fi