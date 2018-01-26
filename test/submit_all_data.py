#!/usr/bin/env python
#usage: python submit_all_data.py -c topplusdmTrees_skim_cfg.py -f samples.txt -d crab_workdir
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR ')

    parser = OptionParser(usage=usage)    
    parser.add_option("-c", "--config", dest="config",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-d", "--dir", dest="dir",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    parser.add_option("-n", "--dryrun", dest="dryrun",
        help=("Dry run"),
        metavar="DRYRUN")
    (options, args) = parser.parse_args()


    if options.config == None or options.dir == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()

    from WMCore.Configuration import Configuration
    config = Configuration()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    config.section_("General")
    config.General.workArea = options.dir
   
    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.config
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.pyCfgParams = ["isData=False", "changeJECs=True"]
    jecversion="V4"
    jecera="Summer16_23Sep2016"
#    jecversion="V2"
#    jecera="Spring16_23Sep2016"
    config.JobType.inputFiles = [jecera+"BCD"+jecversion+"_DATA_UncertaintySources_AK4PFchs.txt" ,
                                 jecera+jecversion+"_MC_L2L3Residual_AK4PFchs.txt",#New MC
                                 jecera+jecversion+"_MC_L3Absolute_AK4PFchs.txt",
                                 jecera+jecversion+"_MC_L1FastJet_AK4PFchs.txt",
                                 jecera+jecversion+"_MC_L2Relative_AK4PFchs.txt",
                                 jecera+"BCD"+jecversion+"_DATA_L1FastJet_AK4PFchs.txt",#RunBCD
                                 jecera+"BCD"+jecversion+"_DATA_L1RC_AK4PFchs.txt",
                                 jecera+"BCD"+jecversion+"_DATA_L2Relative_AK4PFchs.txt",
                                 jecera+"BCD"+jecversion+"_DATA_L3Absolute_AK4PFchs.txt",
                                 jecera+"BCD"+jecversion+"_DATA_L2L3Residual_AK4PFchs.txt",
                                 jecera+"EF"+jecversion+"_DATA_UncertaintySources_AK4PFchs.txt",#RunE
                                 jecera+"EF"+jecversion+"_DATA_L1FastJet_AK4PFchs.txt",
                                 jecera+"EF"+jecversion+"_DATA_L1RC_AK4PFchs.txt",
                                 jecera+"EF"+jecversion+"_DATA_L2Relative_AK4PFchs.txt",
                                 jecera+"EF"+jecversion+"_DATA_L3Absolute_AK4PFchs.txt",
                                 jecera+"EF"+jecversion+"_DATA_L2L3Residual_AK4PFchs.txt",
                                 jecera+"G"+jecversion+"_DATA_UncertaintySources_AK4PFchs.txt",#RunF
                                 jecera+"G"+jecversion+"_DATA_L1FastJet_AK4PFchs.txt",
                                 jecera+"G"+jecversion+"_DATA_L1RC_AK4PFchs.txt",
                                 jecera+"G"+jecversion+"_DATA_L2Relative_AK4PFchs.txt",
                                 jecera+"G"+jecversion+"_DATA_L3Absolute_AK4PFchs.txt",
                                 jecera+"G"+jecversion+"_DATA_L2L3Residual_AK4PFchs.txt",
                                 jecera+"H"+jecversion+"_DATA_UncertaintySources_AK4PFchs.txt",#RunGH
                                 jecera+"H"+jecversion+"_DATA_L1FastJet_AK4PFchs.txt",
                                 jecera+"H"+jecversion+"_DATA_L1RC_AK4PFchs.txt",
                                 jecera+"H"+jecversion+"_DATA_L2Relative_AK4PFchs.txt",
                                 jecera+"H"+jecversion+"_DATA_L3Absolute_AK4PFchs.txt",
                                 jecera+"H"+jecversion+"_DATA_L2L3Residual_AK4PFchs.txt",
                                 jecera+jecversion+"_MC_UncertaintySources_AK4PFchs.txt",
                                 "cMVAv2_Moriond17_B_H.csv",
                                 "CSVv2_Moriond17_B_H.csv",
                                 "btagging_cmva.root",
                                 ]

    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.ignoreLocality = True
    #    config.Data.outLFNDirBase = '/store/user/oiorio/ttDM/samples/2016/Oct/'
    #    config.Data.outLFNDirBase = '/store/user/oiorio/ttDM/trees/2016/Nov/22Nov/'
    config.Data.outLFNDirBase = '/store/user/oiorio/SingleTop/2017/Julsy/13July/'
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
    #    config.Data.totalUnits = -1
    config.Data.unitsPerJob = 10
    config.Data.publication = False

    config.section_("Site")
    config.Site.storageSite = 'T2_IT_Pisa'
    #config.Site.whitelist = ['T2_CH_CERN','T2_IT_*','T2_DE_*','T2_CH_*']


    print 'Using config ' + options.config
    print 'Writing to directory ' + options.dir

    def getJECV(name):
        if "Run2016E" in name: return "DATA2016EF"
        if "Run2016F" in name: return "DATA2016EF"
        if "Run2016G" in name: return "DATA2016G"
        if "Run2016H" in name: return "DATA2016H"
        else: return "DATA2016BCD"
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute commend'
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s
        
        #if s.find('ST_t-') > 0 :
        #   print ' single top\n   => add lhesource=source'
        #        config.JobType.pyCfgParams = ["isData=False", "changeJECs=True","lhesource=source"]
        
    for ijob, job in enumerate(jobs) :

        ptbin = job.split('/')[1]
        cd = (job.split('/')[2]).split('-')[1:-1]
        cond = ""
        for f in cd: cond = cond +f+"-"
        cond = cond[:-1]
        print cond
    
        print "jobsplit11 ", job.split('/')[2]
        print "cfgparams befor", config.JobType.pyCfgParams
        channelData=getJECV(job.split('/')[2])
        config.JobType.pyCfgParams = ["isData=True", "changeJECs=True","channel="+channelData]
        print "cfgparams after", config.JobType.pyCfgParams

        #print cond.split('-')[1:3]
        #        for f in (job.split('/')[2].split['-'][2:4]):
        #           print f
        #B2GAnaFW_80x_V2p0
        #        ptbin
        config.General.requestName = '80xV2_' + ptbin + cond 
        config.Data.inputDataset = job
        config.Data.outputDatasetTag = '80xV2_' + ptbin + cond

        #        config.Data.inputFiles=

        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config
        try :
            from multiprocessing import Process
            print options.dryrun
            if options.dryrun == None:
                p = Process(target=submit, args=(config,))
                p.start()
                p.join()
            else: print 'dry run option set to True: run not submitted!'    
            #submit(config)
        except :
            print 'Not submitted.'
        





if __name__ == '__main__':
    main()            
