outDir="Out_$(date +"%d-%m-%Y_%H-%M")"
mkdir $outDir

###########################   Signal  #########################
./rootcom mutau_analyzer analyze_mutau

./MakeCondorFiles.csh analyze_mutau root://cmsxrootd.hep.wisc.edu//store/user/jmadhusu/with_boosted_taus/2017/zprimeBaryonic/Signal_Zpbaryonic2017_01_2.root Zpbaryonic2017_2.root -1 1000 2017 MC Zpbaryonic2017_2 $outDir
./MakeCondorFiles.csh analyze_mutau root://cmsxrootd.hep.wisc.edu//store/user/jmadhusu/with_boosted_taus/2017/zprimeBaryonic/Signal_Zpbaryonic2017_01_7.root Zpbaryonic2017_7.root -1 1000 2017 MC Zpbaryonic2017_7 $outDir
./MakeCondorFiles.csh analyze_mutau root://cmsxrootd.hep.wisc.edu//store/user/jmadhusu/with_boosted_taus/2017/zprimeBaryonic/Signal_Zpbaryonic2017_01_10.root Zpbaryonic2017_10.root -1 1000 2017 MC Zpbaryonic2017_10 $outDir
./MakeCondorFiles.csh analyze_mutau root://cmsxrootd.hep.wisc.edu//store/user/jmadhusu/with_boostedtau/2017_skimmed/mutau/DYJetsToLL_M-50_TuneCP5.root DYJetsToLL_M-50_TuneCP5.root -1 1000 2017 MC DYJetsToLL_M-50_TuneCP5 $outDir