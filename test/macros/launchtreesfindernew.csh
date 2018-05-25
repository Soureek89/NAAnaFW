#cd $HOME
cd ../../../../../../CMSSW_8_1_0_pre12/src/
cmsenv
cd -
#python treesFinder.py -g gfal -F txt -f samples_mc_edmntuples_Summer16.txt -x xrootd.ba.infn.it -o files/mc -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2017/Apr/25Apr/ -V 2 > & mcfinder.log &
#python treesFinder.py -g gfal -F txt -f samples_mc_edmntuples_sd_Summer16.txt -x xrootd.ba.infn.it -o files/mcsd -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2017/Apr/26Apr/ -V 2 > & mcfinder.log &
python treesFinder.py -g gfal -F txt -f samples_data_mu_edmntuples_Summer16.txt -D True -o files/data/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2017/Apr/25Apr/ -V 2 > & datafinder.log &
#cmsenv
