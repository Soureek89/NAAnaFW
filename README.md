## Full description of the current package. Point 1 is the info for compilation, points 2 through 4 are the detailed info for the analysis. A shortened version of the analysis steps is reported at the end in part 5.


## Part 1: Environment setup and compilation  ##

-- setenv SCRAM_ARCH slc6_amd64_gcc530

-- cd CMSSW_8_0_26_patch1/src

-- cmsenv

-- source /cvmfs/cms.cern.ch/crab3/crab.csh

-- cmsenv

-- mkdir Analysis/NAAnaFW

-- git clone https://github.com/oiorio/NAAnaFW.git Analysis/NAAnaFW

-- scram b -j 10 > & compilation.log &

## Part 2: Tree running  ##
Everything is in the NAAnaFW/test folder:

-- cd test 

##MC:
nohup cmsRun topplusdmTrees_skim_cfg.py
 #runs a skimmed sample, just the nominal trees (no jes/jer), preselection of >=1 tight lepton, >=1 pfLoose jet with pt >40 GeV

##Data:
nohup cmsRun topplusdmTrees_skim_cfg.py isData=True

##Finding files for trees:
cd macros

python treesFinder.py -g gfal -f samples_HLTv1_muons.txt -s storage01.lcg.cscs.ch:8443/srm/managerv2 -o samples/mc/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2016/Oct/ -o files/

NOTA BENE: CURRENTLY THIS DOES NOT WORK ON LXPLUS AND CMSSW_8_0_X BECAUSE OF MISSING lcg-ls COMMANDS! WILL NEED TO USE gfal!!!
NOTA EVEN MORE BENE: Unfortunately, gfal does not work on CMSSW_8_0_X, so you'll have to migrate to 81X to make this script work:
TO MAKE IT WORK PLEASE DO 

-- cd $HOME 

-- cmsrel CMSSW_8_1_0_pre12 

-- cd CMSSW_8_1_0_pre12/src

-- cmsenv

This one fetches recursively all directories in a path [option -p] in the storage element [option -s] looking for samples indicated in 
a txt file [option -f], with some parsing options (veto  a word, look for a particular date etc. ). 
The output are python files including the list of samples, as list of strings with the xrootd paths to all.

Syntax:

-- python treesFinder.py [-f crabFilesToFetch] [-s storageElement] [-p basePathOnTheSE] [-o outputDirToStoreFiles] [-d dateInTheYYMMDD_HHMMSSCrabFormat] [-v VetoWord] [-V verboseLevelFrom0To2]

Can also use the built-in help function:

-- python treesFinder.py --help 


In practice you create a folder for example with

-- cd ../../../../
-- cmsrel CMSSW_8_1_0_pre12

and launch with the script in:

-- cd Analysis/NAAnaFWtest/macros/
-- voms-proxy-init --voms cms
-- source launchtreesfinder.csh

## Part 3: Running the detailed selection and analysis  

Everything is stored in the NAAnaFW/bin folder: cd bin

To rename the output to the standard naming of the analysis please follow the convention found at: 

script_rename.py -s files/ -l files/renamed/

The SingleTopAnalysis.cpp contains the event selection and the systematics application.

Before using it do the 

source preparesetup.csh

to setup the grid certificates to run on the batch queue. 

The SingleTopAnalysis.cpp be used with the following python script:

python new_singletop.py --t3batch -f trees/mc/renamed/ -P ST,TT,VJ,VV -S 10

This launches the analysis on batch queues at Cern by splitting it in groups of 10 root files. Other options for launching the analysis are are: 

-- The "t3se" option will run it interactively on t3 Storage Element

-- The "local" will run the job locally

Examples are therefore

1- python new_singletop.py --t3batch -f trees/mc/renamed/ -P ST,TT,VJ,VV -S 10

2- python new_singletop.py -m t3se -P ... -> launches interactively taking files from the grid

3- python new_singletop.py  -m local -l localfiles/ST  -P... ->launches interactively on the files in the local folder -l blabla

The -P option will take the input processes specified in the samplesST.py, samplesTT.py, samplesVJ.py etc files, also separated by comma.
the results are put in the ./res folder, unless differently specified, in the form of root files like ST_T_muon.root, etc.

If the -S NJobs option is specified, the output will be split in parts ST_T_partN_muon.root.

The parts can be merged with the

-- merge_res.py -l ./res -P ST,TT,VJ,VV --rm True

- which will merge them together for ALL systematics available and, if with the --rm True option is specified, will remove all the separate *_partN_* files.

## Part 4: Making the plots
Everything is stored in the NAAnaFW/bin folder.

Once trees are merged together the makeplot.py can be used.  It allows to perform several operations on the trees. Let's see how it works. NOTE: for all the following operation the trees must be put in the directory trees and the subdirectories: muon, muonantiiso, electron and electronantiiso.

The first operation is the normalitation to the luminosity with the command:

-- python makeplot.py --lumi -L muon,muonantiiso,electron,electronantiiso

The -L option stands for the lepton channel to work out 
The output will be put in the directory trees_lumi which will contain the same four subdirectories of trees.

The second operation is the merging of the components of the samples:

-- python makeplot.py --mertree -L muon,muonantiiso,electron,electronantiiso

These trees can be now used for the trainings of the BDTs: go into the MVA folder and choose the training to perfom in the train_trees.py script; then use the command:

-- python train_trees.py

After all the trainings are be performed move again in bin: cd ..

At this point the MVAs variables can be included into the trees: go in makeplot.py and select the trainings that must be included into the trees adding them at trainings_njmt, for muon, and trainings_njmt, for electron. Then use the command:

-- python makeplot.py --mva -L muon,muonantiiso,electron,electronantiiso

All the previous operations have to be performed only once.
Now it's time to plot. Selecting the variables of interest in makeplot.py one can plot them using the command:

-- python makeplot.py -T 2j1t -C "mtw>50&&etajprime<2.4" -S "noSyst" --plot -L electronantiiso

!!!NOTE!!! For this case the electronantiiso choice has to be worked out alone beacuse the extra cut "mlb>30" must be required.

Other options that can be used are:

-T njmt that select the topology of the event; default perform all the topologies;

-C "selection" to plot the variables with the additional cut selection;

-S "syst" select the systematic to plot: the choices are -S "noSyst" to plot only nominal, -S "syst" to plot only the systematic syst or nothing for the default plotting of all the systematics.   

Then it's possible to make the stacks with the command:

-- python makeplot.py -T 2j1t -C "mtw>50&&etajprime<2.4" -L electron,muon

with this command the DDQCD sample is also produced.

To make the 2D histos of the two BDTs for each region you have to use the command:

-- python makeplot.py --h2D

and then you can unroll them in 1D histos with the command:

-- python makeplot.py --h1D

At the end use the command:

-- python fit_setup.py 

to adjust all the systematics and the plots for the final fit.


## Part 5: Statistical inference
The fit framework can be downloaded following the instructions that you can find here:

https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit

TO be continued when we will commit the repository on GitHub.

Go to the folder stat where you can find the script for running the combine tool. 
The fit to the QCD scale factor can be performed running the script DDQCD_fit.csh wih the command: 

-- source DDQCD_fit.csh

The final fit can be performed running the script fit_muele.csh with the command:

-- source fit_muele.csh

## Summary of the analyis steps:
- 0 Compile following the instructions in Part 1.
*
- 1 Retrieve the files location with:
cd test/macros/
python treesFinder.py -g gfal -F txt -f samples_HLTv1_muons.txt -o ../../bin/files/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2016/Oct/31Oct/ -d 161031 -V 2 > & mcfinder.log &
python treesFinder.py -g gfal -F txt -f samples_data_mu_edmndtuple.txt -D True -o ../../bin/files/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/ser/oiorio/ttDM/trees/2016/Oct/12Oct/SingleMuon/ -V 2 > & datafinder.log &

following instructions here:

test/macros/launchtreesfinder.csh

- 2 Rename the files location according to the convention we follow:

cd ../../bin

python script_rename.py  -s files/ -d files/renamed/

- 3 Launch the batch queues, splitting it wherever needed:

if need be change the xrd to one of your liking (it can improve drastically the speed of your access to the file via xrootd):

python script_replacexrd.py -f files/renamed/ -o files/final/ -x xrootd.ba.infn.it -P ST,VV,VJ,TT,SingleMuon

python new_singletop.py --t3batch -f files/final/ -P ST,TT,VJ,VV -S 10

python new_singletop.py --t3batch -f files/final/ -P SingleMuon -d DATA -S 10

- 4 merge the batch queue result

merge_res.py -l ./res -P ST,TT,VJ,VV --rm True

- 5 make the plots with the makePlot

- 6 Statistical inference with necromantic magic.
