# Boosted DiHiggs analyzer

# Setup

```
# release does not matter..
cmsrel CMSSW_10_6_19 
cd CMSSW_9_4_2/src
cmsenv
git clone git@github.com:cmantill/HHBoostedAnalyzer.git
cd HHBoostedAnalyzer
make
```

# Submit jobs

Make
```
mkdir -p condor
voms-proxy-init --voms cms
```

```
# to produce condor scripts
python scripts/submit_analyzer_jobs_LPC.py
# to submit - go into the folder created and do:
for i in * ; do cd ${i}; condor_submit task.jdl; cd -; done
```

# To normalize (after hadding the output files)
List the dataset name and the filename in a .txt file such as: inputlist.txt:
```
HHToBBVVToBBQQQQ_node_SM /eos/uscms/store/user/cmantill/analyzer/HHToBBWWNtupler/option1/nano/v1/2017/HHToBBVVToBBQQQQ_node_SM.root
```
and then:
```
./NormalizeNtuple inputlist.txt 1
```