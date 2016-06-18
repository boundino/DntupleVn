PDNUMBER=1

INPUTFILE="/data/wangj/Data2015/Dntuple/vN/Dntuple_crab_PbPb_HIMinimumBias${PDNUMBER}_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016.root"
OUTPUTFILE="/data/wangj/Data2015/Dntuple/vN/Dntuple_crab_PbPb_HIMinimumBias${PDNUMBER}_tkpt0p5eta1p5_Dy1p1_EvtPlaneCalibration_v2v3etagap_05142016_skim.root"

g++ skimTree_DzeroVn.C $(root-config --cflags --libs) -g -o skimTree_DzeroVn.exe 
./skimTree_DzeroVn.exe "$INPUTFILE" "$OUTPUTFILE" "$PDNUMBER" 30 50
rm skimTree_DzeroVn.exe
