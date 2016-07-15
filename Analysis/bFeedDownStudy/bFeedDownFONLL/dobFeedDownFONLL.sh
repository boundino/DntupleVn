#!/bin/bash

CENTMIN=30
CENTMAX=50

IS_DCA=1
#
DO_FONLL=0
DO_MCEFFICIENCY=0
DO_PLOTPNNP=0
DO_BFEEDDOWNFONLL=0
#
DO_PLOTFRACTIONS=1

##

INPUTDATA="/data/wangj/Data2015/Dntuple/DVn/Dntuple_crab_PbPb_HIMinimumBias1-7_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016_skim_cent${CENTMIN}${CENTMAX}.root"
INPUTMCP="/data/HeavyFlavourRun2/MC2015/Dntuple/PbPb/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
INPUTMCNP="/data/HeavyFlavourRun2/MC2015/Dntuple/PbPb/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_nonprompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"

SELGEN="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))"
CUT="pclusterCompatibilityFilter&&pprimaryVertexFilter&&phfCoincFilter3&&abs(PVz)<15&&Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>1.0&&Dtrk2Pt>1.0&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3&&abs(Dtrk1Eta)<1.5&&abs(Dtrk2Eta)<1.5&&((DlxyBS/DlxyBSErr)>2.5&&Dalpha<0.12&&((Dpt>1&&Dpt<2&&(DsvpvDistance/DsvpvDisErr)>6.0&&Dchi2cl>0.25)||(Dpt>2&&Dpt<4&&(DsvpvDistance/DsvpvDisErr)>5.86&&Dchi2cl>0.224)||(Dpt>4&&Dpt<5&&(DsvpvDistance/DsvpvDisErr)>5.46&&Dchi2cl>0.196)||(Dpt>5&&Dpt<6&&(DsvpvDistance/DsvpvDisErr)>4.86&&Dchi2cl>0.170)||(Dpt>6&&Dpt<8&&(DsvpvDistance/DsvpvDisErr)>4.54&&Dchi2cl>0.125)||(Dpt>8&&Dpt<10&&(DsvpvDistance/DsvpvDisErr)>4.42&&Dchi2cl>0.091)||(Dpt>10&&Dpt<15&&(DsvpvDistance/DsvpvDisErr)>4.06&&Dchi2cl>0.069)||(Dpt>15&&Dpt<20&&(DsvpvDistance/DsvpvDisErr)>3.71&&Dchi2cl>0.056)||(Dpt>20&&Dpt<25&&(DsvpvDistance/DsvpvDisErr)>3.25&&Dchi2cl>0.054)||(Dpt>25&&(DsvpvDistance/DsvpvDisErr)>2.97&&Dchi2cl>0.050)))"
WEIGHT="(pthatweight)*(maxDgenpt<pthat/1.2)*(6.14981+hiBin*(-0.156513)+hiBin*hiBin*(0.00149127)+hiBin*hiBin*hiBin*(-6.29087e-06)+hiBin*hiBin*hiBin*hiBin*(9.90029e-09))*(-0.00600791+maxDgenpt*(0.0838585)+maxDgenpt*maxDgenpt*(-0.00991096)+maxDgenpt*maxDgenpt*maxDgenpt*(0.000496019)+maxDgenpt*maxDgenpt*maxDgenpt*maxDgenpt*(-8.50065e-06))"

#

TFEND="inclusive"
CUTEND=""
if [ $IS_DCA -eq 1 ]; then
TFEND="DCA"
CUTEND="&&(DsvpvDistance*TMath::Sin(Dalpha)>0.008)"
fi

#

INPUTFONLLP="fonlls/FONLL_pp_promptDzero_5TeV_y1.dat"
INPUTFONLLNP="fonlls/FONLL_pp_nonpromptDzero_5TeV_y1.dat"
OUTPUTFONLLP="outfiles/FONLL_pp_promptDzero_5TeV_y1.root"
OUTPUTFONLLNP="outfiles/FONLL_pp_nonpromptDzero_5TeV_y1.root"

OUTPUTMCEFFP="outfiles/MCefficiencyP"
OUTPUTMCEFFNP="outfiles/MCefficiencyNP"

OUTPUTFRACTION="outfilesResult/bFeedDownFONLL"
OUTPUTFRACTIONcompare="../bFeedDown/outfilesResult/bFeedDownResult"

##

if [ $DO_FONLL -eq 1 ]; then 
g++ fonllDsigmadpt.C $(root-config --cflags --libs) -g -o fonllDsigmadpt.exe 
./fonllDsigmadpt.exe "$INPUTFONLLP" "$OUTPUTFONLLP"
./fonllDsigmadpt.exe "$INPUTFONLLNP" "$OUTPUTFONLLNP"
rm fonllDsigmadpt.exe
fi

if [ $DO_MCEFFICIENCY -eq 1 ]; then      
g++ MCefficiency.C $(root-config --cflags --libs) -g -o MCefficiency.exe 
./MCefficiency.exe "$INPUTMCP" "$OUTPUTMCEFFP" "$TFEND" "$SELGEN" "${CUT}${CUTEND}" "$WEIGHT" "$CENTMIN" "$CENTMAX"
./MCefficiency.exe "$INPUTMCNP" "$OUTPUTMCEFFNP" "$TFEND" "$SELGEN" "${CUT}${CUTEND}" "$WEIGHT" "$CENTMIN" "$CENTMAX"
rm MCefficiency.exe
fi
if [ $DO_PLOTPNNP -eq 1 ]; then      
g++ plotPnNP.C $(root-config --cflags --libs) -g -o plotPnNP.exe 
./plotPnNP.exe "$OUTPUTMCEFFP" "$OUTPUTMCEFFNP" "$TFEND" "$CENTMIN" "$CENTMAX"
rm plotPnNP.exe
fi

if [ $DO_BFEEDDOWNFONLL -eq 1 ]; then      
g++ bFeedDownFONLL.C $(root-config --cflags --libs) -g -o bFeedDownFONLL.exe 
./bFeedDownFONLL.exe "$OUTPUTFONLLP" "$OUTPUTFONLLNP" "$OUTPUTMCEFFP" "$OUTPUTMCEFFNP" "$OUTPUTFRACTION" "$TFEND" "$CENTMIN" "$CENTMAX"
rm bFeedDownFONLL.exe
fi

if [ $DO_PLOTFRACTIONS -eq 1 ]; then      
g++ plotFractions.C $(root-config --cflags --libs) -g -o plotFractions.exe 
./plotFractions.exe "$OUTPUTFRACTION" "$OUTPUTFRACTIONcompare" "$CENTMIN" "$CENTMAX"
rm plotFractions.exe
fi
