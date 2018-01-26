#cd $HOME
cd ../../../../../../CMSSW_8_1_0_pre12/src/
cmsenv
cd -

#python treesFinder.py -g gfal -F txt -f samples_edmntuple_mc_full_Summer16.tex -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mc/ -p /cms/store/user/oiorio/SingleTop/2017/July/13July/ -V 2 > & mcfinderpisa.log &


#python treesFinder.py -g gfal -F txt -f samples_data_mu_pt1_ReReReco_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/datamu/ -D true -p /cms/store/user/oiorio/SingleTop/2017/Julsy/13July/SingleMuon/ -V 2 > & datafinderpisamu1.log &

#python treesFinder.py -g gfal -F txt -f samples_data_mu_pt2_ReReReco_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/datamu/ -D true -p /cms/store/user/oiorio/SingleTop/2017/Julsy/13July/SingleMuon/ -V 2 > & datafinderpisamu2.log &

#python treesFinder.py -g gfal -F txt -f samples_data_ele_ReReReco_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/dataele/ -D true -p /cms/store/user/oiorio/SingleTop/2017/Julsy/13July/SingleElectron/ -V 2 > & datafinderpisael.log &

#python treesFinder.py -g gfal -F txt -f ../samples_mc_edmntuples_tch_sd_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mcsd/ -p /cms/store/user/oiorio/SingleTop/2017/Juls/25July/ -V 2 > & datafinderpisatchsd.log &

python treesFinder.py -g gfal -F txt -f ../samples_mc_edmntuples_tch_sd_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mcsd/ -p /cms/store/user/oiorio/SingleTop/2017/Aug/27Aug/ -V 2 > & datafinderpisatchpsd.log &
#python treesFinder.py -g gfal -F txt -f ../samples_mc_edmntuples_tt_sd_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mcsd/  -p /cms/store/user/oiorio/SingleTop/2017/Juls/25July/ -V 2 > & datafinderpisattsd.log &


#python treesFinder.py -g gfal -F txt -f edmntupples_syst_mc.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mcsyst/ -p /cms/store/user/oiorio/SingleTop/2017/July/26uly/ -V 2 > & datafinderpisasyst.log &
#python treesFinder.py -g gfal -F txt -f edmntupples_syst_one_mc.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mcsyst/ -p /cms/store/user/oiorio/SingleTop/2017/July/26uly/ -V 2 > & datafinderpisasyst.log &


#python treesFinder.py -g gfal -F txt -f samples_mc_edmntuples_sd_loc_Summer16.txt -s stormfe1.pi.infn.it:8443/srm/managerv2 -x cms-xrd-global.cern.ch -o files/mcsd/ -p /cms/store/user/oiorio/SingleTop/2017/June/17June/ -V 2 > & mcfinderpisasdv2.log &



