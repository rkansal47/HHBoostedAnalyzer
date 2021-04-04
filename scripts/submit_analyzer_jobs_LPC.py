#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict

queueType = "longlunch"
option = 1
#label = "option1"
label = "option1_ak15_ak8_cut"

analysis = "HHToBBWWNtupler"
outputfile = "HHToBBWWNtuple" + "_" + label

cmsswReleaseVersion = "CMSSW_10_6_19"

outputDirectoryBase = "/store/user/rkansal/analyzer/"+analysis+"/"+label+"/"

datasetList = OrderedDict()

#2017  hqu ntuples
datasetList['nano/v2_ak15/2017/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.list'] = [0, 2, "2017", "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8"]
datasetList['nano/v2_ak15/2017/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"]
datasetList['nano/v2_ak15/2017/TTToHadronic_TuneCP5_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "TTToHadronic_TuneCP5_13TeV-powheg-pythia8"]
datasetList['nano/v2_ak15/2017/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"]
# datasetList['nano/v2_ak15/2017/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8"]
# datasetList['nano/v2_ak15/2017/WWTo4Q_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "WWTo4Q_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8"]
# datasetList['nano/v2_ak15/2017/WW_TuneCP5_13TeV-pythia8.list'] = [0, 2, "2017", "WW_TuneCP5_13TeV-pythia8"]
# datasetList['nano/v2_ak15/2017/WZ_TuneCP5_13TeV-pythia8.list'] = [0, 2, "2017", "WZ_TuneCP5_13TeV-pythia8"]
# datasetList['nano/v2_ak15/2017/ZZ_TuneCP5_13TeV-pythia8.list'] = [0, 2, "2017", "ZZ_TuneCP5_13TeV-pythia8"]
# datasetList['nano/v2_ak15/2017/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"]
# datasetList['nano/v2_ak15/2017/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8"]
# datasetList['nano/v2_ak15/2017/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"]
# datasetList['nano/v2_ak15/2017/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8"]
#datasetList['nano/v2_ak15/2017/HHToBBVVToBBQQQQ_node_SM.list'] = [0, 2, "2017", "HHToBBVVToBBQQQQ_node_SM"]
#datasetList['nano/v2_ak15/2017/HHToVVVV_node_SM_Pt300.list'] = [0, 2, "2017", "HHToVVVV_node_SM_Pt300"]
datasetList['nano/v2_ak15/2017/HHToBBVVToBBQQQQ_cHHH1.list'] = [0, 2, "2017", "GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"]
datasetList['nano/v2_ak15/2017/GluGluToHHTo4V_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"]
#datasetList['nano/v2_ak15/2017/GluGluHToBB_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8.list'] = [0, 2, "2017", "GluGluHToBB_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8"]
#datasetList['nano/v2_ak15/2017/GluGluToHHTo4B_node_SM.list'] = [0, 2, "2017", "GluGluToHHTo4B_node_SM"]
# datasetList['nano/v2_ak15/2017/GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.list'] = [0, 2, "2017", "GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"]
datasetList['nano/v2_ak15/2017/JetHT_2017B.list'] = [1, 1, "2017", ""]
datasetList['nano/v2_ak15/2017/JetHT_2017C.list'] = [1, 1, "2017", ""]
datasetList['nano/v2_ak15/2017/JetHT_2017D.list'] = [1, 1, "2017", ""]
datasetList['nano/v2_ak15/2017/JetHT_2017E.list'] = [1, 1, "2017", ""]
datasetList['nano/v2_ak15/2017/JetHT_2017F.list'] = [1, 1, "2017", ""]

CMSSW_BASE_DIR = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE_DIR+"/src/HHBoostedAnalyzer/"

#create directory for condor jobs

for listfile in datasetList.keys():

    datasetName = listfile.replace(".list","")
    print "Preparing analyzer workflow for dataset :" + datasetName + "\n"
    if not os.path.exists(Analyzer_DIR+"/list/" + listfile):
        print "listfile: " + listfile + " does not exist. skipping."
        continue

    outputDirectory = outputDirectoryBase + datasetName + "/"
    tmpListFile = open(Analyzer_DIR + "/list/" + listfile,"r")

    year = datasetList[listfile][2]
    sampleName = datasetList[listfile][3]

    #####################################
    #Job Splitting
    #####################################
    isData = "no"
    if (datasetList[listfile][0] == 1):
        isData = "yes"
    filesPerJob = datasetList[listfile][1]
    tmpJobFileCount = 0
    nJobs = 1

    if os.path.exists(Analyzer_DIR+"/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt"):
        print "Warning: condor directory " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + " is not empty. Skipping."
        continue

    #create condor directories
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName )
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/log/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/out/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/err/")

    tmpOutputListFile = open( Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt","w")
    for line in tmpListFile:

        #open list file for new job
        if tmpJobFileCount >= filesPerJob:
            tmpOutputListFile.close()
            tmpJobFileCount = 0
            nJobs = nJobs + 1
            tmpOutputListFile = open( Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt","w")

        #write input file into job list file
        tmpOutputListFile.write(line)
        tmpJobFileCount += 1

    tmpOutputListFile.close()
    os.system("cd " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/; tar czf input_list.tgz input_list_*.txt")

    #####################################
    #Copy run script and executable
    #####################################
    os.system("cp " + Analyzer_DIR + "/scripts/run_job_LPC.sh " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("cp " + Analyzer_DIR + "Run" + analysis + " " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_" + year + ".root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/PileupWeights/")
    os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupWeights.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/PileupWeights/")


    #####################################
    #Create Condor JDL file
    #####################################
    tmpCondorJDLFile = open(Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/task.jdl","w+")
    tmpCondorJDLFileTemplate = """
Universe  = vanilla
Executable = ./run_job_LPC.sh
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("Arguments = " + analysis + " " + str(isData) + " " + str(option) + " " + "$(I) " + outputfile + " " + outputDirectory + " " + cmsswReleaseVersion + " " + year + " " + sampleName + "\n")

    tmpCondorJDLFileTemplate = """
Log = log/job.$(Cluster).$(Process).log
Output = out/job.$(Cluster).$(Process).out
Error = err/job.$(Cluster).$(Process).err
x509userproxy = $ENV(X509_USER_PROXY)
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("transfer_input_files = " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/run_job_LPC.sh, " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/input_list.tgz, " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/Run" + analysis + ", " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JetHTTriggerEfficiency_" + year + ".root, " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/PileupWeights/PileupWeights.root " + "\n")

    tmpCondorJDLFileTemplate = """
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Resources request
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("RequestMemory = 2000 \n")

    tmpCondorJDLFileTemplate = """

# Jobs selection
Queue I from (
"""

    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    for i in range(1,nJobs+1):
        tmpCondorJDLFile.write(str(i)+"\n")
    tmpCondorJDLFile.write(")\n")
    tmpCondorJDLFile.close()
