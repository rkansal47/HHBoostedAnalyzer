#include "HHToBBWWNtupler.h"
#include "JetTree.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void HHToBBWWNtupler::Analyze(bool isData, int Option, string outputfilename, string year, string pileupWeightName, bool addAK15)
{
 
    cout << "Initializing..." << endl;

    //----------------------------------------
    //Load auxiliary information
    //----------------------------------------  
    TH1F *pileupWeightHist = 0;
    
    if (!isData) {
      string CMSSWDir = std::getenv("CMSSW_BASE");

      string pileupWeightFilename = CMSSWDir + "/src/HHBoostedAnalyzer/data/PileupWeights/PileupWeights.root";
      TFile *pileupWeightFile = new TFile(pileupWeightFilename.c_str(),"READ");
      if (!pileupWeightFile) {
	cout << "Warning : pileupWeightFile " << pileupWeightFile << " could not be opened.\n";  
      } else {
	cout << "Opened pileupWeightFile " << pileupWeightFilename << "\n"; 
      }
      string pileupWeightHistname = "PUWeight_" + pileupWeightName + "_" + year;
      if (pileupWeightFile) {
	pileupWeightHist = (TH1F*)pileupWeightFile->Get(pileupWeightHistname.c_str());
      } 
      if (pileupWeightHist) {
	cout << "Found pileupWeightHist " << pileupWeightHistname << "in file " << pileupWeightFilename << "\n";
      } else {
	cout << "Warning :  could not find pileupWeightHist named " 
	     << pileupWeightHistname 
	     << " in file " << pileupWeightFilename << "\n";
      }
    }

    //----------------------------------------
    //Output file
    //----------------------------------------  
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "HHToBBWWNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");    
 
    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //output TTree
    TTree *outputTree = new TTree("tree", "");
 
    //------------------------
    //declare branch variables
    //------------------------  
    float weight = 0;
    float pileupWeight = 0;
    float totalWeight = 0;

    float genHiggs1Pt = -1;
    float genHiggs1Eta = -1;
    float genHiggs1Phi = -1;
    float genHiggs2Pt = -1;
    float genHiggs2Eta = -1;
    float genHiggs2Phi = -1;
    float genHH_pt = -99;
    float genHH_eta = -99;
    float genHH_phi = -99;
    float genHH_mass = -99;
    float genLeptonPt = -1;
    float genLeptonEta = -1;
    float genLeptonPhi = -1;
    int   genLeptonId = 0;
    int   genLeptonMotherId = 0;

    float genHiggs1W1Pt = -1;
    float genHiggs1W1Eta = -1;
    float genHiggs1W1Phi = -1;
    float genHiggs1W1M = -1;
    int genHiggs1W1Decay = -1;
    float genHiggs1W2Pt = -1;
    float genHiggs1W2Eta = -1;
    float genHiggs1W2Phi = -1;
    float genHiggs1W2M = -1;
    int genHiggs1W2Decay = -1;
    float genHiggs2W1Pt = -1;
    float genHiggs2W1Eta = -1;
    float genHiggs2W1Phi = -1;
    float genHiggs2W1M = -1;
    int genHiggs2W1Decay = -1;
    float genHiggs2W2Pt = -1;
    float genHiggs2W2Eta = -1;
    float genHiggs2W2Phi = -1;
    float genHiggs2W2M = -1;
    int genHiggs2W2Decay = -1;

    float genHiggs1b1Pt = -1;
    float genHiggs1b1Eta = -1;
    float genHiggs1b1Phi = -1;
    float genHiggs1b2Pt = -1;
    float genHiggs1b2Eta = -1;
    float genHiggs1b2Phi = -1;
    float genHiggs2b1Pt = -1;
    float genHiggs2b1Eta = -1;
    float genHiggs2b1Phi = -1;
    float genHiggs2b2Pt = -1;
    float genHiggs2b2Eta = -1;
    float genHiggs2b2Phi = -1;

    float genHiggs1W1dau1Pt = -1;
    float genHiggs1W1dau1Eta = -1;
    float genHiggs1W1dau1Phi = -1;
    int genHiggs1W1dau1Id = 0;
    float genHiggs1W1dau2Pt = -1;
    float genHiggs1W1dau2Eta = -1;
    float genHiggs1W1dau2Phi = -1;
    int genHiggs1W1dau2Id = 0;
    float genHiggs1W2dau1Pt = -1;
    float genHiggs1W2dau1Eta = -1;
    float genHiggs1W2dau1Phi = -1;
    int genHiggs1W2dau1Id = 0;
    float genHiggs1W2dau2Pt = -1;
    float genHiggs1W2dau2Eta = -1;
    float genHiggs1W2dau2Phi = -1;
    int genHiggs1W2dau2Id = 0;

    float genHiggs2W1dau1Pt = -1;
    float genHiggs2W1dau1Eta = -1;
    float genHiggs2W1dau1Phi = -1;
    int genHiggs2W1dau1Id = 0;
    float genHiggs2W1dau2Pt = -1;
    float genHiggs2W1dau2Eta = -1;
    float genHiggs2W1dau2Phi = -1;
    int genHiggs2W1dau2Id = 0;
    float genHiggs2W2dau1Pt = -1;
    float genHiggs2W2dau1Eta = -1;
    float genHiggs2W2dau1Phi = -1;
    int genHiggs2W2dau1Id = 0;
    float genHiggs2W2dau2Pt = -1;
    float genHiggs2W2dau2Eta = -1;
    float genHiggs2W2dau2Phi = -1;
    int genHiggs2W2dau2Id = 0;

    int NJets = 0;
    float MET = -1;
    float fatJet1Pt = -99;
    float fatJet1Eta = -99;
    float fatJet1Phi = -99;
    float fatJet1Mass = -99;
    float fatJet1MassSD = -99;
    float fatJet1DeepAK8_H = -99;
    float fatJet1DeepAK8MD_H4qvsQCD = -99;
    float fatJet1PNetXbb = -99;
    float fatJet1PNetXbb_alt = -99;
    float fatJet1PNetQCDb = -99;
    float fatJet1PNetQCDbb = -99;
    float fatJet1PNetQCDc = -99;
    float fatJet1PNetQCDcc = -99;
    float fatJet1PNetQCDothers = -99;
    int   fatJet1GenMatchIndex = -99;
    float fatJet1Tau4OverTau3 = -99;
    float fatJet1Tau4OverTau2 = -99;
    float fatJet1Tau4OverTau1 = -99;
    float fatJet1Tau2OverTau1 = -99;
    float fatJet1n2b1 = -99; 
    float fatJet1lsf3 = -99;
    bool fatJet1HasMuon = 0;
    bool fatJet1HasElectron = 0;
    bool fatJet1HasBJetCSVLoose = 0;
    bool fatJet1HasBJetCSVMedium = 0;
    bool fatJet1HasBJetCSVTight = 0;
    bool fatJet1OppositeHemisphereHasBJet = 0;

    float fatJet2Pt = -99;
    float fatJet2Eta = -99;
    float fatJet2Phi = -99;
    float fatJet2Mass = -99;
    float fatJet2MassSD = -99;
    float fatJet2DeepAK8_H = -99;
    float fatJet2DeepAK8MD_H4qvsQCD = -99;
    float fatJet2PNetXbb = -99;
    float fatJet2PNetXbb_alt = -99;
    float fatJet2PNetQCDb = -99;
    float fatJet2PNetQCDbb = -99;
    float fatJet2PNetQCDc = -99;
    float fatJet2PNetQCDcc = -99;
    float fatJet2PNetQCDothers = -99;
    int   fatJet2GenMatchIndex = -99;
    float fatJet2Tau4OverTau3 = -99;
    float fatJet2Tau4OverTau2 = -99;
    float fatJet2Tau4OverTau1 = -99;
    float fatJet2Tau2OverTau1 = -99;
    float fatJet2n2b1 = -99;
    float fatJet2lsf3 = -99;
    bool fatJet2HasMuon = 0;
    bool fatJet2HasElectron = 0;
    bool fatJet2HasBJetCSVLoose = 0;
    bool fatJet2HasBJetCSVMedium = 0;
    bool fatJet2HasBJetCSVTight = 0;

    float fatJet3Pt = -99;
    float fatJet3Eta = -99;
    float fatJet3Phi = -99;
    float fatJet3Mass = -99;
    float fatJet3MassSD = -99;
    float fatJet3DeepAK8_H = -99;
    float fatJet3DeepAK8MD_H4qvsQCD = -99;
    float fatJet3PNetXbb = -99;
    float fatJet3PNetXbb_alt = -99;
    float fatJet3PNetQCDb = -99;
    float fatJet3PNetQCDbb = -99;
    float fatJet3PNetQCDc = -99;
    float fatJet3PNetQCDcc = -99;
    float fatJet3PNetQCDothers = -99;
    float fatJet3Tau4OverTau3 = -99;
    float fatJet3Tau4OverTau2 = -99;
    float fatJet3Tau4OverTau1 = -99;
    float fatJet3Tau2OverTau1 = -99;
    bool fatJet3HasMuon = 0;
    bool fatJet3HasElectron = 0;
    bool fatJet3HasBJetCSVLoose = 0;
    bool fatJet3HasBJetCSVMedium = 0;
    bool fatJet3HasBJetCSVTight = 0;

    float ak15fatJet1Pt = -99;
    float ak15fatJet1Eta = -99;
    float ak15fatJet1Phi = -99;
    float ak15fatJet1Mass = -99;
    float ak15fatJet1MassSD = -99;
    float ak15fatJet1PNetMDXbb = -99;
    float ak15fatJet1PNetHqqqq = -99;
    float ak15fatJet1PNetQCDb = -99;
    float ak15fatJet1PNetQCDbb = -99;
    float ak15fatJet1PNetQCDc = -99;
    float ak15fatJet1PNetQCDcc = -99;
    float ak15fatJet1PNetQCDothers = -99;
    
    float ak15fatJet2Pt = -99;
    float ak15fatJet2Eta = -99;
    float ak15fatJet2Phi = -99;
    float ak15fatJet2Mass = -99;
    float ak15fatJet2MassSD = -99;
    float ak15fatJet2PNetMDXbb = -99;
    float ak15fatJet2PNetHqqqq = -99;
    float ak15fatJet2PNetQCDb = -99;
    float ak15fatJet2PNetQCDbb = -99;
    float ak15fatJet2PNetQCDc = -99;
    float ak15fatJet2PNetQCDcc = -99;
    float ak15fatJet2PNetQCDothers = -99;

    float hh_pt = -99;
    float hh_eta = -99;
    float hh_phi = -99;
    float hh_mass = -99;        
    float fatJet1PtOverMHH = -99;
    float fatJet1PtOverMSD = -99;
    float fatJet2PtOverMHH = -99;
    float fatJet2PtOverMSD = -99;
    float deltaEta_j1j2 = -99;
    float deltaPhi_j1j2 = -99;
    float deltaR_j1j2 = -99;    
    float ptj2_over_ptj1 = -99;
    float mj2_over_mj1 = -99;
    float lep1Pt = -99;
    float lep1Eta = -99;
    float lep1Phi = -99;
    int   lep1Id = 0;
    float lep2Pt = -99;
    float lep2Eta = -99;
    float lep2Phi = -99;
    int   lep2Id = 0;
    float pho1Pt = -99;
    float pho1Eta = -99;
    float pho1Phi = -99;
    float jet1Pt = -99;
    float jet1Eta = -99;
    float jet1Phi = -99;
    float jet1DeepJetBTag = -99;
    float jet2Pt = -99;
    float jet2Eta = -99;
    float jet2Phi = -99;
    float jet2DeepJetBTag = -99;
    float jet3Pt = -99;
    float jet3Eta = -99;
    float jet3Phi = -99;
    float jet3DeepJetBTag = -99;
    float jet4Pt = -99;
    float jet4Eta = -99;
    float jet4Phi = -99;
    float jet4DeepJetBTag = -99;
    int   nBTaggedJets = 0;

    //------------------------
    //set branches on big tree
    //------------------------
    outputTree->Branch("weight", &weight, "weight/F");
    outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
    outputTree->Branch("totalWeight", &totalWeight, "totalWeight/F");
    outputTree->Branch("run", &run, "run/i");
    outputTree->Branch("lumi", &luminosityBlock, "lumi/i");
    outputTree->Branch("event", &event, "event/l");
    outputTree->Branch("npu", &Pileup_nTrueInt, "npu/F");
    outputTree->Branch("rho", &fixedGridRhoFastjetAll, "rho/F");
    outputTree->Branch("MET", &MET, "MET/F");

    outputTree->Branch("jet1Pt", &jet1Pt, "jet1Pt/F");
    outputTree->Branch("jet1Eta", &jet1Eta, "jet1Eta/F");
    outputTree->Branch("jet1Phi", &jet1Phi, "jet1Phi/F");
    outputTree->Branch("jet1DeepJetBTag", &jet1DeepJetBTag, "jet1DeepJetBTag/F");      
    outputTree->Branch("jet2Pt", &jet2Pt, "jet2Pt/F");
    outputTree->Branch("jet2Eta", &jet2Eta, "jet2Eta/F");
    outputTree->Branch("jet2Phi", &jet2Phi, "jet2Phi/F");
    outputTree->Branch("jet2DeepJetBTag", &jet2DeepJetBTag, "jet2DeepJetBTag/F");      
    outputTree->Branch("jet3Pt", &jet3Pt, "jet3Pt/F");
    outputTree->Branch("jet3Eta", &jet3Eta, "jet3Eta/F");
    outputTree->Branch("jet3Phi", &jet3Phi, "jet3Phi/F");
    outputTree->Branch("jet3DeepJetBTag", &jet3DeepJetBTag, "jet3DeepJetBTag/F");      
    outputTree->Branch("jet4Pt", &jet4Pt, "jet4Pt/F");
    outputTree->Branch("jet4Eta", &jet4Eta, "jet4Eta/F");
    outputTree->Branch("jet4Phi", &jet4Phi, "jet4Phi/F");
    outputTree->Branch("jet4DeepJetBTag", &jet4DeepJetBTag, "jet4DeepJetBTag/F");      

    outputTree->Branch("genHiggs1Pt", &genHiggs1Pt, "genHiggs1Pt/F");
    outputTree->Branch("genHiggs1Eta", &genHiggs1Eta, "genHiggs1Eta/F");
    outputTree->Branch("genHiggs1Phi", &genHiggs1Phi, "genHiggs1Phi/F");
    outputTree->Branch("genHiggs2Pt", &genHiggs2Pt, "genHiggs2Pt/F");
    outputTree->Branch("genHiggs2Eta", &genHiggs2Eta, "genHiggs2Eta/F");
    outputTree->Branch("genHiggs2Phi", &genHiggs2Phi, "genHiggs2Phi/F");
    outputTree->Branch("genHH_pt",      &genHH_pt,     "genHH_pt/F");
    outputTree->Branch("genHH_eta",     &genHH_eta,    "genHH_eta/F");
    outputTree->Branch("genHH_phi",     &genHH_phi,    "genHH_phi/F");
    outputTree->Branch("genHH_mass",    &genHH_mass,   "genHH_mass/F");
    outputTree->Branch("genLeptonId", &genLeptonId, "genLeptonId/I");
    outputTree->Branch("genLeptonMotherId", &genLeptonMotherId, "genLeptonMotherId/I");
    outputTree->Branch("genLeptonPt", &genLeptonPt, "genLeptonPt/F");
    outputTree->Branch("genLeptonEta", &genLeptonEta, "genLeptonEta/F");
    outputTree->Branch("genLeptonPhi", &genLeptonPhi, "genLeptonPhi/F");
    
    outputTree->Branch("genHiggs1W1Pt", &genHiggs1W1Pt, "genHiggs1W1Pt/F");
    outputTree->Branch("genHiggs1W1Eta", &genHiggs1W1Eta, "genHiggs1W1Eta/F");
    outputTree->Branch("genHiggs1W1Phi", &genHiggs1W1Phi, "genHiggs1W1Phi/F");
    outputTree->Branch("genHiggs1W1M", &genHiggs1W1M, "genHiggs1W1M/F");
    outputTree->Branch("genHiggs1W1Decay", &genHiggs1W1Decay, "genHiggs1W1Decay/I");
    outputTree->Branch("genHiggs1W2Pt", &genHiggs1W2Pt, "genHiggs1W2Pt/F");
    outputTree->Branch("genHiggs1W2Eta", &genHiggs1W2Eta, "genHiggs1W2Eta/F");
    outputTree->Branch("genHiggs1W2Phi", &genHiggs1W2Phi, "genHiggs1W2Phi/F");
    outputTree->Branch("genHiggs1W2M", &genHiggs1W2M, "genHiggs1W2M/F");
    outputTree->Branch("genHiggs1W2Decay", &genHiggs1W2Decay, "genHiggs1W2Decay/I");

    outputTree->Branch("genHiggs2W1Pt", &genHiggs2W1Pt, "genHiggs2W1Pt/F");
    outputTree->Branch("genHiggs2W1Eta", &genHiggs2W1Eta, "genHiggs2W1Eta/F");
    outputTree->Branch("genHiggs2W1Phi", &genHiggs2W1Phi, "genHiggs2W1Phi/F");
    outputTree->Branch("genHiggs2W1M", &genHiggs2W1M, "genHiggs2W1M/F");
    outputTree->Branch("genHiggs2W1Decay", &genHiggs2W1Decay, "genHiggs2W1Decay/I");
    outputTree->Branch("genHiggs2W2Pt", &genHiggs2W2Pt, "genHiggs2W2Pt/F");
    outputTree->Branch("genHiggs2W2Eta", &genHiggs2W2Eta, "genHiggs2W2Eta/F");
    outputTree->Branch("genHiggs2W2Phi", &genHiggs2W2Phi, "genHiggs2W2Phi/F");
    outputTree->Branch("genHiggs2W2M", &genHiggs2W2M, "genHiggs2W2M/F");
    outputTree->Branch("genHiggs2W2Decay", &genHiggs2W2Decay, "genHiggs2W2Decay/I");

    outputTree->Branch("genHiggs1b1Pt", &genHiggs1b1Pt, "genHiggs1b1Pt/F");
    outputTree->Branch("genHiggs1b1Eta", &genHiggs1b1Eta, "genHiggs1b1Eta/F");
    outputTree->Branch("genHiggs1b1Phi", &genHiggs1b1Phi, "genHiggs1b1Phi/F");
    outputTree->Branch("genHiggs1b2Pt", &genHiggs1b2Pt, "genHiggs1b2Pt/F");
    outputTree->Branch("genHiggs1b2Eta", &genHiggs1b2Eta, "genHiggs1b2Eta/F");
    outputTree->Branch("genHiggs1b2Phi", &genHiggs1b2Phi, "genHiggs1b2Phi/F");
    outputTree->Branch("genHiggs2b1Pt", &genHiggs2b1Pt, "genHiggs2b1Pt/F");
    outputTree->Branch("genHiggs2b1Eta", &genHiggs2b1Eta, "genHiggs2b1Eta/F");
    outputTree->Branch("genHiggs2b1Phi", &genHiggs2b1Phi, "genHiggs2b1Phi/F");
    outputTree->Branch("genHiggs2b2Pt", &genHiggs2b2Pt, "genHiggs2b2Pt/F");
    outputTree->Branch("genHiggs2b2Eta", &genHiggs2b2Eta, "genHiggs2b2Eta/F");
    outputTree->Branch("genHiggs2b2Phi", &genHiggs2b2Phi, "genHiggs2b2Phi/F");

    outputTree->Branch("genHiggs1W1dau1Pt", &genHiggs1W1dau1Pt, "genHiggs1W1dau1Pt/F");
    outputTree->Branch("genHiggs1W1dau1Eta", &genHiggs1W1dau1Eta,"genHiggs1W1dau1Eta/F");
    outputTree->Branch("genHiggs1W1dau1Phi", &genHiggs1W1dau1Phi,"genHiggs1W1dau1Phi/F");
    outputTree->Branch("genHiggs1W1dau1Id", &genHiggs1W1dau1Id,"genHiggs1W1dau1Id/I");
    outputTree->Branch("genHiggs1W1dau2Pt", &genHiggs1W1dau2Pt,"genHiggs1W1dau2Pt/F");
    outputTree->Branch("genHiggs1W1dau2Eta", &genHiggs1W1dau2Eta,"genHiggs1W1dau2Eta/F");
    outputTree->Branch("genHiggs1W1dau2Phi", &genHiggs1W1dau2Phi,"genHiggs1W1dau2Phi/F");
    outputTree->Branch("genHiggs1W1dau2Id", &genHiggs1W1dau2Id,"genHiggs1W1dau2Id/I");
    outputTree->Branch("genHiggs1W2dau1Pt", &genHiggs1W2dau1Pt,"genHiggs1W2dau1Pt/F");
    outputTree->Branch("genHiggs1W2dau1Eta", &genHiggs1W2dau1Eta,"genHiggs1W2dau1Eta/F");
    outputTree->Branch("genHiggs1W2dau1Phi", &genHiggs1W2dau1Phi,"genHiggs1W2dau1Phi/F");
    outputTree->Branch("genHiggs1W2dau1Id", &genHiggs1W2dau1Id,"genHiggs1W2dau1Id/I");
    outputTree->Branch("genHiggs1W2dau2Pt", &genHiggs1W2dau2Pt,"genHiggs1W2dau2Pt/F");
    outputTree->Branch("genHiggs1W2dau2Eta", &genHiggs1W2dau2Eta,"genHiggs1W2dau2Eta/F");
    outputTree->Branch("genHiggs1W2dau2Phi", &genHiggs1W2dau2Phi,"genHiggs1W2dau2Phi/F");
    outputTree->Branch("genHiggs1W2dau2Id", &genHiggs1W2dau2Id,"genHiggs1W2dau2Id/I");

    outputTree->Branch("genHiggs2W1dau1Pt", &genHiggs2W1dau1Pt,"genHiggs2W1dau1Pt/F");
    outputTree->Branch("genHiggs2W1dau1Eta", &genHiggs2W1dau1Eta,"genHiggs2W1dau1Eta/F");
    outputTree->Branch("genHiggs2W1dau1Phi", &genHiggs2W1dau1Phi,"genHiggs2W1dau1Phi/F");
    outputTree->Branch("genHiggs2W1dau1Id", &genHiggs2W1dau1Id,"genHiggs2W1dau1Id/I");
    outputTree->Branch("genHiggs2W1dau2Pt", &genHiggs2W1dau2Pt,"genHiggs2W1dau2Pt/F");
    outputTree->Branch("genHiggs2W1dau2Eta", &genHiggs2W1dau2Eta,"genHiggs2W1dau2Eta/F");
    outputTree->Branch("genHiggs2W1dau2Phi", &genHiggs2W1dau2Phi,"genHiggs2W1dau2Phi/F");
    outputTree->Branch("genHiggs2W1dau2Id", &genHiggs2W1dau2Id,"genHiggs2W1dau2Id/I");
    outputTree->Branch("genHiggs2W2dau1Pt", &genHiggs2W2dau1Pt,"genHiggs2W2dau1Pt/F");
    outputTree->Branch("genHiggs2W2dau1Eta", &genHiggs2W2dau1Eta,"genHiggs2W2dau1Eta/F");
    outputTree->Branch("genHiggs2W2dau1Phi", &genHiggs2W2dau1Phi,"genHiggs2W2dau1Phi/F");
    outputTree->Branch("genHiggs2W2dau1Id", &genHiggs2W2dau1Id,"genHiggs2W2dau1Id/I");
    outputTree->Branch("genHiggs2W2dau2Pt", &genHiggs2W2dau2Pt,"genHiggs2W2dau2Pt/F");
    outputTree->Branch("genHiggs2W2dau2Eta", &genHiggs2W2dau2Eta,"genHiggs2W2dau2Eta/F");
    outputTree->Branch("genHiggs2W2dau2Phi", &genHiggs2W2dau2Phi,"genHiggs2W2dau2Phi/F");
    outputTree->Branch("genHiggs2W2dau2Id", &genHiggs2W2dau2Id,"genHiggs2W2dau2Id/I");

    outputTree->Branch("NJets", &NJets, "NJets/I");
    outputTree->Branch("fatJet1Pt", &fatJet1Pt, "fatJet1Pt/F");
    outputTree->Branch("fatJet1Eta", &fatJet1Eta, "fatJet1Eta/F");
    outputTree->Branch("fatJet1Phi", &fatJet1Phi, "fatJet1Phi/F");
    outputTree->Branch("fatJet1Mass", &fatJet1Mass, "fatJet1Mass/F");
    outputTree->Branch("fatJet1MassSD", &fatJet1MassSD, "fatJet1MassSD/F");
    outputTree->Branch("fatJet1DeepAK8_H", &fatJet1DeepAK8_H, "fatJet1DeepAK8_H/F");
    outputTree->Branch("fatJet1DeepAK8MD_H4qvsQCD", &fatJet1DeepAK8MD_H4qvsQCD, "fatJet1DeepAK8MD_H4qvsQCD/F");
    outputTree->Branch("fatJet1PNetXbb", &fatJet1PNetXbb, "fatJet1PNetXbb/F");
    outputTree->Branch("fatJet1PNetXbb_alt", &fatJet1PNetXbb_alt, "fatJet1PNetXbb_alt/F");
    outputTree->Branch("fatJet1PNetQCDb", &fatJet1PNetQCDb, "fatJet1PNetQCDb/F");
    outputTree->Branch("fatJet1PNetQCDbb", &fatJet1PNetQCDbb, "fatJet1PNetQCDbb/F");
    outputTree->Branch("fatJet1PNetQCDc", &fatJet1PNetQCDc, "fatJet1PNetQCDc/F");
    outputTree->Branch("fatJet1PNetQCDcc", &fatJet1PNetQCDcc, "fatJet1PNetQCDcc/F");
    outputTree->Branch("fatJet1PNetQCDothers", &fatJet1PNetQCDothers, "fatJet1PNetQCDothers/F");
    outputTree->Branch("fatJet1GenMatchIndex", &fatJet1GenMatchIndex, "fatJet1GenMatchIndex/I");
    outputTree->Branch("fatJet1Tau4OverTau3", &fatJet1Tau4OverTau3, "fatJet1Tau4OverTau3/F");
    outputTree->Branch("fatJet1Tau4OverTau2", &fatJet1Tau4OverTau2, "fatJet1Tau4OverTau2/F");
    outputTree->Branch("fatJet1Tau4OverTau1", &fatJet1Tau4OverTau1, "fatJet1Tau4OverTau1/F");
    outputTree->Branch("fatJet1Tau2OverTau1", &fatJet1Tau2OverTau1, "fatJet1Tau2OverTau1/F");
    outputTree->Branch("fatJet1n2b1", &fatJet1n2b1, "fatJet1n2b1/F");
    outputTree->Branch("fatJet1lsf3", &fatJet1lsf3, "fatJet1lsf3/F");
    outputTree->Branch("fatJet1HasMuon", &fatJet1HasMuon, "fatJet1HasMuon/O");
    outputTree->Branch("fatJet1HasElectron", &fatJet1HasElectron, "fatJet1HasElectron/O");
    outputTree->Branch("fatJet1HasBJetCSVLoose", &fatJet1HasBJetCSVLoose, "fatJet1HasBJetCSVLoose/O");
    outputTree->Branch("fatJet1HasBJetCSVMedium", &fatJet1HasBJetCSVMedium, "fatJet1HasBJetCSVMedium/O");
    outputTree->Branch("fatJet1HasBJetCSVTight", &fatJet1HasBJetCSVTight, "fatJet1HasBJetCSVTight/O");
    outputTree->Branch("fatJet1OppositeHemisphereHasBJet", &fatJet1OppositeHemisphereHasBJet, "fatJet1OppositeHemisphereHasBJet/O");

    outputTree->Branch("fatJet2Pt", &fatJet2Pt, "fatJet2Pt/F");
    outputTree->Branch("fatJet2Eta", &fatJet2Eta, "fatJet2Eta/F");
    outputTree->Branch("fatJet2Phi", &fatJet2Phi, "fatJet2Phi/F");
    outputTree->Branch("fatJet2Mass", &fatJet2Mass, "fatJet2Mass/F");
    outputTree->Branch("fatJet2MassSD", &fatJet2MassSD, "fatJet2MassSD/F");
    outputTree->Branch("fatJet2DeepAK8_H", &fatJet2DeepAK8_H, "fatJet2DeepAK8_H/F");
    outputTree->Branch("fatJet2DeepAK8MD_H4qvsQCD", &fatJet2DeepAK8MD_H4qvsQCD,"fatJet2DeepAK8MD_H4qvsQCD/F");
    outputTree->Branch("fatJet2PNetXbb", &fatJet2PNetXbb, "fatJet2PNetXbb/F");
    outputTree->Branch("fatJet2PNetXbb_alt", &fatJet2PNetXbb_alt, "fatJet2PNetXbb_alt/F");
    outputTree->Branch("fatJet2PNetQCDb", &fatJet2PNetQCDb, "fatJet2PNetQCDb/F");
    outputTree->Branch("fatJet2PNetQCDbb", &fatJet2PNetQCDbb, "fatJet2PNetQCDbb/F");
    outputTree->Branch("fatJet2PNetQCDc", &fatJet2PNetQCDc, "fatJet2PNetQCDc/F");
    outputTree->Branch("fatJet2PNetQCDcc", &fatJet2PNetQCDcc, "fatJet2PNetQCDcc/F");
    outputTree->Branch("fatJet2PNetQCDothers", &fatJet2PNetQCDothers, "fatJet2PNetQCDothers/F");
    outputTree->Branch("fatJet2GenMatchIndex", &fatJet2GenMatchIndex, "fatJet2GenMatchIndex/I");
    outputTree->Branch("fatJet2Tau4OverTau3", &fatJet2Tau4OverTau3, "fatJet2Tau4OverTau3/F");
    outputTree->Branch("fatJet2Tau4OverTau2", &fatJet2Tau4OverTau2, "fatJet2Tau4OverTau2/F");
    outputTree->Branch("fatJet2Tau4OverTau1", &fatJet2Tau4OverTau1, "fatJet2Tau4OverTau1/F");
    outputTree->Branch("fatJet2Tau2OverTau1", &fatJet2Tau2OverTau1, "fatJet2Tau2OverTau1/F");
    outputTree->Branch("fatJet2n2b1", &fatJet2n2b1, "fatJet2n2b1/F");
    outputTree->Branch("fatJet2lsf3", &fatJet2lsf3, "fatJet2lsf3/F");
    outputTree->Branch("fatJet2HasMuon", &fatJet2HasMuon, "fatJet2HasMuon/O");
    outputTree->Branch("fatJet2HasElectron", &fatJet2HasElectron, "fatJet2HasElectron/O");
    outputTree->Branch("fatJet2HasBJetCSVLoose", &fatJet2HasBJetCSVLoose, "fatJet2HasBJetCSVLoose/O");
    outputTree->Branch("fatJet2HasBJetCSVMedium", &fatJet2HasBJetCSVMedium, "fatJet2HasBJetCSVMedium/O");
    outputTree->Branch("fatJet2HasBJetCSVTight", &fatJet2HasBJetCSVTight, "fatJet2HasBJetCSVTight/O");
    outputTree->Branch("fatJet3Pt", &fatJet3Pt, "fatJet3Pt/F");
    outputTree->Branch("fatJet3Eta", &fatJet3Eta, "fatJet3Eta/F");
    outputTree->Branch("fatJet3Phi", &fatJet3Phi, "fatJet3Phi/F");
    outputTree->Branch("fatJet3Mass", &fatJet3Mass, "fatJet3Mass/F");
    outputTree->Branch("fatJet3MassSD", &fatJet3MassSD, "fatJet3MassSD/F");
    outputTree->Branch("fatJet3DeepAK8_H", &fatJet3DeepAK8_H, "fatJet3DeepAK8_H/F");
    outputTree->Branch("fatJet3DeepAK8MD_H4qvsQCD", &fatJet3DeepAK8MD_H4qvsQCD,"fatJet3DeepAK8MD_H4qvsQCD/F");
    outputTree->Branch("fatJet3PNetXbb", &fatJet3PNetXbb, "fatJet3PNetXbb/F");
    outputTree->Branch("fatJet3PNetXbb_alt", &fatJet3PNetXbb_alt, "fatJet3PNetXbb_alt/F");
    outputTree->Branch("fatJet3PNetQCDb", &fatJet3PNetQCDb, "fatJet3PNetQCDb/F");
    outputTree->Branch("fatJet3PNetQCDbb", &fatJet3PNetQCDbb, "fatJet3PNetQCDbb/F");
    outputTree->Branch("fatJet3PNetQCDc", &fatJet3PNetQCDc, "fatJet3PNetQCDc/F");
    outputTree->Branch("fatJet3PNetQCDcc", &fatJet3PNetQCDcc, "fatJet3PNetQCDcc/F");
    outputTree->Branch("fatJet3PNetQCDothers", &fatJet3PNetQCDothers, "fatJet3PNetQCDothers/F");
    outputTree->Branch("fatJet3Tau4OverTau3", &fatJet3Tau4OverTau3, "fatJet3Tau4OverTau3/F");
    outputTree->Branch("fatJet3Tau4OverTau2", &fatJet3Tau4OverTau2, "fatJet3Tau4OverTau2/F");
    outputTree->Branch("fatJet3Tau4OverTau1", &fatJet3Tau4OverTau1, "fatJet3Tau4OverTau1/F");
    outputTree->Branch("fatJet3Tau2OverTau1", &fatJet3Tau2OverTau1, "fatJet3Tau2OverTau1/F");
    outputTree->Branch("fatJet3HasMuon", &fatJet3HasMuon, "fatJet3HasMuon/O");
    outputTree->Branch("fatJet3HasElectron", &fatJet3HasElectron, "fatJet3HasElectron/O");
    outputTree->Branch("fatJet3HasBJetCSVLoose", &fatJet3HasBJetCSVLoose, "fatJet3HasBJetCSVLoose/O");
    outputTree->Branch("fatJet3HasBJetCSVMedium", &fatJet3HasBJetCSVMedium, "fatJet3HasBJetCSVMedium/O");
    outputTree->Branch("fatJet3HasBJetCSVTight", &fatJet3HasBJetCSVTight, "fatJet3HasBJetCSVTight/O");
    outputTree->Branch("hh_pt",      &hh_pt,     "hh_pt/F");
    outputTree->Branch("hh_eta",     &hh_eta,    "hh_eta/F");
    outputTree->Branch("hh_phi",     &hh_phi,    "hh_phi/F");
    outputTree->Branch("hh_mass",    &hh_mass,   "hh_mass/F");
    outputTree->Branch("fatJet1PtOverMHH",    &fatJet1PtOverMHH,   "fatJet1PtOverMHH/F");
    outputTree->Branch("fatJet1PtOverMSD",    &fatJet1PtOverMSD,   "fatJet1PtOverMSD/F");
    outputTree->Branch("fatJet2PtOverMHH",    &fatJet2PtOverMHH,   "fatJet2PtOverMHH/F");
    outputTree->Branch("fatJet2PtOverMSD",    &fatJet2PtOverMSD,   "fatJet2PtOverMSD/F");
    outputTree->Branch("deltaEta_j1j2",    &deltaEta_j1j2,   "deltaEta_j1j2/F");
    outputTree->Branch("deltaPhi_j1j2",    &deltaPhi_j1j2,   "deltaPhi_j1j2/F");
    outputTree->Branch("deltaR_j1j2",    &deltaR_j1j2,   "deltaR_j1j2/F");
    outputTree->Branch("ptj2_over_ptj1",    &ptj2_over_ptj1,   "ptj2_over_ptj1/F");
    outputTree->Branch("mj2_over_mj1",    &mj2_over_mj1,   "mj2_over_mj1/F");

    if(addAK15){
      outputTree->Branch("ak15fatJet1Pt", &ak15fatJet1Pt, "ak15fatJet1Pt/F");
      outputTree->Branch("ak15fatJet1Eta", &ak15fatJet1Eta, "ak15fatJet1Eta/F");
      outputTree->Branch("ak15fatJet1Phi", &ak15fatJet1Phi, "ak15fatJet1Phi/F");
      outputTree->Branch("ak15fatJet1Mass", &ak15fatJet1Mass, "ak15fatJet1Mass/F");
      outputTree->Branch("ak15fatJet1MassSD", &ak15fatJet1MassSD, "ak15fatJet1MassSD/F");
      outputTree->Branch("ak15fatJet1PNetMDXbb", &ak15fatJet1PNetMDXbb, "ak15fatJet1PNetMDXbb/F");
      outputTree->Branch("ak15fatJet1PNetHqqqq", &ak15fatJet1PNetHqqqq, "ak15fatJet1PNetHqqqq/F");
      outputTree->Branch("ak15fatJet1PNetQCDb", &ak15fatJet1PNetQCDb, "ak15fatJet1PNetQCDb/F");
      outputTree->Branch("ak15fatJet1PNetQCDbb", &ak15fatJet1PNetQCDbb, "ak15fatJet1PNetQCDbb/F");
      outputTree->Branch("ak15fatJet1PNetQCDc", &ak15fatJet1PNetQCDc, "ak15fatJet1PNetQCDc/F");
      outputTree->Branch("ak15fatJet1PNetQCDcc", &ak15fatJet1PNetQCDcc, "ak15fatJet1PNetQCDcc/F");
      outputTree->Branch("ak15fatJet1PNetQCDothers", &ak15fatJet1PNetQCDothers, "ak15fatJet1PNetQCDothers/F");

      outputTree->Branch("ak15fatJet2Pt", &ak15fatJet2Pt, "ak15fatJet2Pt/F");
      outputTree->Branch("ak15fatJet2Eta", &ak15fatJet2Eta, "ak15fatJet2Eta/F");
      outputTree->Branch("ak15fatJet2Phi", &ak15fatJet2Phi, "ak15fatJet2Phi/F");
      outputTree->Branch("ak15fatJet2Mass", &ak15fatJet2Mass, "ak15fatJet2Mass/F");
      outputTree->Branch("ak15fatJet2MassSD", &ak15fatJet2MassSD, "ak15fatJet2MassSD/F");
      outputTree->Branch("ak15fatJet2PNetMDXbb", &ak15fatJet2PNetMDXbb, "ak15fatJet2PNetMDXbb/F");
      outputTree->Branch("ak15fatJet2PNetHqqqq", &ak15fatJet2PNetHqqqq,"ak15fatJet2PNetHqqqq/F");
      outputTree->Branch("ak15fatJet2PNetQCDb", &ak15fatJet2PNetQCDb, "ak15fatJet2PNetQCDb/F");
      outputTree->Branch("ak15fatJet2PNetQCDbb", &ak15fatJet2PNetQCDbb, "ak15fatJet2PNetQCDbb/F");
      outputTree->Branch("ak15fatJet2PNetQCDc", &ak15fatJet2PNetQCDc, "ak15fatJet2PNetQCDc/F");
      outputTree->Branch("ak15fatJet2PNetQCDcc", &ak15fatJet2PNetQCDcc, "ak15fatJet2PNetQCDcc/F");
      outputTree->Branch("ak15fatJet2PNetQCDothers", &ak15fatJet2PNetQCDothers,"ak15fatJet2PNetQCDothers/F");
    }

    outputTree->Branch("lep1Pt", &lep1Pt, "lep1Pt/F");
    outputTree->Branch("lep1Eta", &lep1Eta, "lep1Eta/F");
    outputTree->Branch("lep1Phi", &lep1Phi, "lep1Phi/F");
    outputTree->Branch("lep1Id", &lep1Id, "lep1Id/I");
    outputTree->Branch("lep2Pt", &lep2Pt, "lep2Pt/F");
    outputTree->Branch("lep2Eta", &lep2Eta, "lep2Eta/F");
    outputTree->Branch("lep2Phi", &lep2Phi, "lep2Phi/F");
    outputTree->Branch("lep2Id", &lep2Id, "lep2Id/I");
    outputTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");      
    outputTree->Branch("pho1Pt", &pho1Pt, "pho1Pt/F");
    outputTree->Branch("pho1Eta", &pho1Eta, "pho1Eta/F");
    outputTree->Branch("pho1Phi", &pho1Phi, "pho1Phi/F");  
      
    outputTree->Branch("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, "HLT_Ele27_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele28_WPTight_Gsf", &HLT_Ele28_WPTight_Gsf, "HLT_Ele28_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf, "HLT_Ele30_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, "HLT_Ele32_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, "HLT_Ele35_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, "HLT_Ele38_WPTight_Gsf/O");
    outputTree->Branch("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, "HLT_Ele40_WPTight_Gsf/O");
    outputTree->Branch("HLT_IsoMu20", &HLT_IsoMu20, "HLT_IsoMu20/O");
    outputTree->Branch("HLT_IsoMu24", &HLT_IsoMu24, "HLT_IsoMu24/O");
    outputTree->Branch("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, "HLT_IsoMu24_eta2p1/O");
    outputTree->Branch("HLT_IsoMu27", &HLT_IsoMu27, "HLT_IsoMu27/O");
    outputTree->Branch("HLT_IsoMu30", &HLT_IsoMu30, "HLT_IsoMu30/O");
    outputTree->Branch("HLT_Mu50", &HLT_Mu50, "HLT_Mu50/O");
    outputTree->Branch("HLT_Mu55", &HLT_Mu55, "HLT_Mu55/O");
    outputTree->Branch("HLT_Photon175",                          &HLT_Photon175,                         "HLT_Photon175/O");
    
    outputTree->Branch("HLT_PFHT780",                                        &HLT_PFHT780,                                       "HLT_PFHT780/O");
    outputTree->Branch("HLT_PFHT890",                                        &HLT_PFHT890,                                       "HLT_PFHT890/O");
    outputTree->Branch("HLT_PFHT1050",                                        &HLT_PFHT1050,                                       "HLT_PFHT1050/O");
    outputTree->Branch("HLT_AK8PFJet360_TrimMass30",                          &HLT_AK8PFJet360_TrimMass30,                         "HLT_AK8PFJet360_TrimMass30/O");
    outputTree->Branch("HLT_AK8PFJet380_TrimMass30",                          &HLT_AK8PFJet380_TrimMass30,                         "HLT_AK8PFJet380_TrimMass30/O");
    outputTree->Branch("HLT_AK8PFJet400_TrimMass30",                          &HLT_AK8PFJet400_TrimMass30,                         "HLT_AK8PFJet400_TrimMass30/O");
    outputTree->Branch("HLT_AK8PFJet420_TrimMass30",                          &HLT_AK8PFJet420_TrimMass30,                         "HLT_AK8PFJet420_TrimMass30/O");
    outputTree->Branch("HLT_AK8PFHT750_TrimMass50",                           &HLT_AK8PFHT750_TrimMass50,                          "HLT_AK8PFHT750_TrimMass50/O");
    outputTree->Branch("HLT_AK8PFHT800_TrimMass50",                           &HLT_AK8PFHT800_TrimMass50,                          "HLT_AK8PFHT800_TrimMass50/O");
    outputTree->Branch("HLT_AK8PFHT850_TrimMass50",                           &HLT_AK8PFHT850_TrimMass50,                          "HLT_AK8PFHT850_TrimMass50/O");
    outputTree->Branch("HLT_AK8PFHT900_TrimMass50",                           &HLT_AK8PFHT900_TrimMass50,                          "HLT_AK8PFHT900_TrimMass50/O");
    outputTree->Branch("HLT_PFJet450",                                        &HLT_PFJet450,                                       "HLT_PFJet450/O");
    outputTree->Branch("HLT_PFJet500",                                        &HLT_PFJet500,                                       "HLT_PFJet500/O");
    outputTree->Branch("HLT_PFJet550",                                        &HLT_PFJet550,                                       "HLT_PFJet550/O");
    outputTree->Branch("HLT_AK8PFJet450",                                     &HLT_AK8PFJet450,                                    "HLT_AK8PFJet450/O");
    outputTree->Branch("HLT_AK8PFJet500",                                     &HLT_AK8PFJet500,                                    "HLT_AK8PFJet500/O");
    outputTree->Branch("HLT_AK8PFJet550",                                     &HLT_AK8PFJet550,                                    "HLT_AK8PFJet550/O");
    outputTree->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17",     &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17,    "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17/O");
    outputTree->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1",      &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1,     "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1/O");
    outputTree->Branch("HLT_AK8PFJet330_PFAK8BTagCSV_p17",                    &HLT_AK8PFJet330_PFAK8BTagCSV_p17,                   "HLT_AK8PFJet330_PFAK8BTagCSV_p17/O");
    outputTree->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02",  &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02/O");
    outputTree->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2",  &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2/O");
    outputTree->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4",  &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4/O"); 
    outputTree->Branch("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20",        &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20,       "HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20/O");
    outputTree->Branch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087",       &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087,      "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087/O");
    outputTree->Branch("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087",       &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087,      "HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087/O");
    outputTree->Branch("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20",     &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20,    "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20/O");
    outputTree->Branch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",        &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20,       "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20/O");
    outputTree->Branch("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20",        &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20,       "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20/O");	

    cout << "Run With Option = " << Option << "\n";
    
    if (Option == 1) cout << "Option = 1 : Select FatJets with pT > 250 GeV and MassSD>20\n";

    UInt_t NEventsFilled = 0;
 
    //begin loop
    if (fChain == 0) return;
    UInt_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    cout << "nentries = " << nentries << "\n";
    for (UInt_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //fill normalization histogram
      weight = genWeight / fabs(genWeight);
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + weight);

      //reset tree variables
      genHiggs1Pt = -99.0;
      genHiggs1Eta = -99.0;
      genHiggs1Phi = -99.0;
      genHiggs2Pt = -99.0;
      genHiggs2Eta = -99.0;
      genHiggs2Phi = -99.0;
      genHH_pt = -99;
      genHH_eta = -99;
      genHH_phi = -99;
      genHH_mass = -99;   
      genHiggs1W1Pt = -99.0;
      genHiggs1W1Eta = -99.0;
      genHiggs1W1Phi = -99.0;
      genHiggs1W1M = -99.0;
      genHiggs1W1Decay = -99;
      genHiggs1W2Pt = -99.0;
      genHiggs1W2Eta = -99.0;
      genHiggs1W2Phi = -99.0;
      genHiggs1W2M = -99.0;
      genHiggs1W2Decay = -99;
      genHiggs2W1Pt = -99.0;
      genHiggs2W1Eta = -99.0;
      genHiggs2W1Phi = -99.0;
      genHiggs2W1M = -99.0;
      genHiggs2W1Decay = -99;
      genHiggs2W2Pt = -99.0;
      genHiggs2W2Eta = -99.0;
      genHiggs2W2Phi = -99.0;
      genHiggs2W2M = -99.0;
      genHiggs2W2Decay = -99;
      genHiggs1b1Pt = -99.0;
      genHiggs1b1Eta = -99.0;
      genHiggs1b1Phi = -99.0;
      genHiggs1b2Pt = -99.0;
      genHiggs1b2Eta = -99.0;
      genHiggs1b2Phi = -99.0;
      genHiggs2b1Pt = -99.0;
      genHiggs2b1Eta = -99.0;
      genHiggs2b1Phi = -99.0;
      genHiggs2b2Pt = -99.0;
      genHiggs2b2Eta = -99.0;
      genHiggs2b2Phi = -99.0;
      genLeptonId = 0;
      genLeptonMotherId = 0;
      genLeptonPt = -99.0;
      genLeptonEta = -99.0;
      genLeptonPhi = -99.0;

      genHiggs1W1dau1Pt = -99.0;
      genHiggs1W1dau1Eta = -99.0;
      genHiggs1W1dau1Phi = -99.0;
      genHiggs1W1dau1Id = 0;
      genHiggs1W1dau2Pt = -99.0;
      genHiggs1W1dau2Eta = -99.0;
      genHiggs1W1dau2Phi = -99.0;
      genHiggs1W1dau2Id = 0;
      genHiggs1W2dau1Pt = -99.0;
      genHiggs1W2dau1Eta = -99.0;
      genHiggs1W2dau1Phi = -99.0;
      genHiggs1W2dau1Id = 0;
      genHiggs1W2dau2Pt = -99.0;
      genHiggs1W2dau2Eta = -99.0;
      genHiggs1W2dau2Phi = -99.0;
      genHiggs1W2dau2Id = 0;

      genHiggs2W1dau1Pt = -99.0;
      genHiggs2W1dau1Eta = -99.0;
      genHiggs2W1dau1Phi = -99.0;
      genHiggs2W1dau1Id = 0;
      genHiggs2W1dau2Pt = -99.0;
      genHiggs2W1dau2Eta = -99.0;
      genHiggs2W1dau2Phi = -99.0;
      genHiggs2W1dau2Id = 0;
      genHiggs2W2dau1Pt = -99.0;
      genHiggs2W2dau1Eta = -99.0;
      genHiggs2W2dau1Phi = -99.0;
      genHiggs2W2dau1Id = 0;
      genHiggs2W2dau2Pt = -99.0;
      genHiggs2W2dau2Eta = -99.0;
      genHiggs2W2dau2Phi = -99.0;
      genHiggs2W2dau2Id = 0;

      NJets = -1;
      MET = -99.0;

      fatJet1Pt = -99.0;
      fatJet1Eta = -99.0;
      fatJet1Phi = -99.0;
      fatJet1Mass = -99.0;
      fatJet1MassSD = -99.0;
      fatJet1DeepAK8_H = -99.0;
      fatJet1DeepAK8MD_H4qvsQCD = -99.0;
      fatJet1PNetXbb = -99;
      fatJet1PNetXbb_alt = -99;
      fatJet1PNetQCDb = -99;
      fatJet1PNetQCDbb = -99;
      fatJet1PNetQCDc = -99;
      fatJet1PNetQCDcc = -99;
      fatJet1PNetQCDothers = -99;
      fatJet1GenMatchIndex = -99.0;
      fatJet1Tau4OverTau3 = -99.0;
      fatJet1Tau4OverTau2 = -99.0;
      fatJet1Tau4OverTau1 = -99.0;
      fatJet1Tau2OverTau1 = -99.0;
      fatJet1n2b1 = -99.0; 
      fatJet1lsf3 = -99.0;
      fatJet1HasMuon = 0;
      fatJet1HasElectron = 0;
      fatJet1HasBJetCSVLoose = 0;
      fatJet1HasBJetCSVMedium = 0;
      fatJet1HasBJetCSVTight = 0;      
      fatJet1OppositeHemisphereHasBJet = 0;
      fatJet2Pt = -99.0;
      fatJet2Eta = -99.0;
      fatJet2Phi = -99.0;
      fatJet2Mass = -99.0;
      fatJet2MassSD = -99.0;
      fatJet2DeepAK8_H = -99.0;
      fatJet2DeepAK8MD_H4qvsQCD = -99.0;
      fatJet2PNetXbb = -99;
      fatJet2PNetXbb_alt = -99;
      fatJet2PNetQCDb = -99;
      fatJet2PNetQCDbb = -99;
      fatJet2PNetQCDc = -99;
      fatJet2PNetQCDcc = -99;
      fatJet2PNetQCDothers = -99;
      fatJet2GenMatchIndex = -99.0;
      fatJet2Tau4OverTau3 = -99.0;
      fatJet2Tau4OverTau2 = -99.0;
      fatJet2Tau4OverTau1 = -99.0;
      fatJet2Tau2OverTau1 = -99.0;
      fatJet2n2b1 = -99.0;
      fatJet2lsf3 = -99.0;
      fatJet2HasMuon = 0;
      fatJet2HasElectron = 0;
      fatJet2HasBJetCSVLoose = 0;
      fatJet2HasBJetCSVMedium = 0;
      fatJet2HasBJetCSVTight = 0;
      fatJet3Pt = -99.0;
      fatJet3Eta = -99.0;
      fatJet3Phi = -99.0;
      fatJet3Mass = -99.0;
      fatJet3MassSD = -99.0;
      fatJet3DeepAK8_H = -99.0;
      fatJet3DeepAK8MD_H4qvsQCD = -99.0;
      fatJet3PNetXbb = -99;      
      fatJet3PNetXbb_alt = -99;
      fatJet3PNetQCDb = -99;
      fatJet3PNetQCDbb = -99;
      fatJet3PNetQCDc = -99;
      fatJet3PNetQCDcc = -99;
      fatJet3PNetQCDothers = -99;
      fatJet3Tau4OverTau3 = -99.0;
      fatJet3Tau4OverTau2 = -99.0;
      fatJet3Tau4OverTau1 = -99.0;
      fatJet3Tau2OverTau1 = -99.0;
      fatJet3HasMuon = 0;
      fatJet3HasElectron = 0;
      fatJet3HasBJetCSVLoose = 0;
      fatJet3HasBJetCSVMedium = 0;
      fatJet3HasBJetCSVTight = 0;

      hh_pt = -99;
      hh_eta = -99;
      hh_phi = -99;
      hh_mass = -99;
      fatJet1PtOverMHH = -99;
      fatJet1PtOverMSD = -99;
      fatJet2PtOverMHH = -99;
      fatJet2PtOverMSD = -99;
      deltaEta_j1j2 = -99;
      deltaPhi_j1j2 = -99;
      deltaR_j1j2 = -99;    
      ptj2_over_ptj1 = -99;
      mj2_over_mj1 = -99;
      lep1Pt = -99;
      lep1Eta = -99;
      lep1Phi = -99;
      lep1Id = 0;
      lep2Pt = -99;
      lep2Eta = -99;
      lep2Phi = -99;
      lep2Id = 0;
      pho1Pt = -99;
      pho1Eta = -99;
      pho1Phi = -99;
      jet1Pt = -99;
      jet1Eta = -99;
      jet1Phi = -99;
      jet1DeepJetBTag = -99;
      jet2Pt = -99;
      jet2Eta = -99;
      jet2Phi = -99;
      jet2DeepJetBTag = -99;
      jet3Pt = -99;
      jet3Eta = -99;
      jet3Phi = -99;
      jet3DeepJetBTag = -99;
      jet4Pt = -99;
      jet4Eta = -99;
      jet4Phi = -99;
      jet4DeepJetBTag = -99;
      nBTaggedJets = 0;

      if(addAK15){
	ak15fatJet1Pt = -99.0;
	ak15fatJet1Eta = -99.0;
	ak15fatJet1Phi = -99.0;
	ak15fatJet1Mass = -99.0;
	ak15fatJet1MassSD = -99.0;
	ak15fatJet1PNetMDXbb = -99.0;
	ak15fatJet1PNetHqqqq = -99.0;
	ak15fatJet1PNetQCDb = -99.0;
	ak15fatJet1PNetQCDbb = -99.0;
	ak15fatJet1PNetQCDc = -99.0;
	ak15fatJet1PNetQCDcc = -99.0;
	ak15fatJet1PNetQCDothers = -99.0;

        ak15fatJet2Pt = -99.0;
        ak15fatJet2Eta = -99.0;
        ak15fatJet2Phi = -99.0;
        ak15fatJet2Mass = -99.0;
        ak15fatJet2MassSD = -99.0;
        ak15fatJet2PNetMDXbb = -99.0;
        ak15fatJet2PNetHqqqq = -99.0;
        ak15fatJet2PNetQCDb = -99.0;
        ak15fatJet2PNetQCDbb = -99.0;
        ak15fatJet2PNetQCDc = -99.0;
        ak15fatJet2PNetQCDcc = -99.0;
        ak15fatJet2PNetQCDothers = -99.0;
      }


      //------------------------------
      //----Event variables------------
      //------------------------------
      MET = MET_pt;


      // Gen variables
      std::vector< TLorentzVector > genHiggsVector;
      if (!isData) {

	// find gen-higgs
	int current_mIndex = -1;
	std::vector< int > genHiggsIndex;
	for(int i = 0; i < nGenPart; i++) {
	  if( (abs(GenPart_pdgId[i]) == 24 || abs(GenPart_pdgId[i]) == 5)  && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 25 && current_mIndex != GenPart_genPartIdxMother[i] ) {
	    //std::cout << GenPart_genPartIdxMother[i] << std::endl;
	    // std::cout << "mother: " << GenPart_pdgId[GenPart_genPartIdxMother[i]]
	    // << " PT: " << GenPart_pt[GenPart_genPartIdxMother[i]]
	    // << " eta: " << GenPart_eta[GenPart_genPartIdxMother[i]]
	    // << " phi: " << GenPart_phi[GenPart_genPartIdxMother[i]] << std::endl;
	    TLorentzVector h;
	    h.SetPtEtaPhiM( GenPart_pt[GenPart_genPartIdxMother[i]], GenPart_eta[GenPart_genPartIdxMother[i]], GenPart_phi[GenPart_genPartIdxMother[i]], GenPart_mass[GenPart_genPartIdxMother[i]] );
	    genHiggsVector.push_back(h);
	    genHiggsIndex.push_back(GenPart_genPartIdxMother[i]);
	    current_mIndex = GenPart_genPartIdxMother[i];
	  }

	  // save leading gen lepton from W decay or tau decay
          if ( (abs(GenPart_pdgId[i]) == 11 || abs(GenPart_pdgId[i]) == 13)
               && GenPart_pt[i] > 10
               && (abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24 || abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 15)
               && GenPart_pt[i] > genLeptonPt
               ) {
            genLeptonId = GenPart_pdgId[i];
            genLeptonMotherId = GenPart_pdgId[GenPart_genPartIdxMother[i]];
            genLeptonPt = GenPart_pt[i];
            genLeptonEta = GenPart_eta[i];
            genLeptonPhi = GenPart_phi[i];
          }
	}

	for(int j = 0; j < genHiggsIndex.size(); j++) {
	  std::vector< TLorentzVector > genWVector, genbVector;
	  std::vector< int > genWIndex;
	  for(int i = 0; i < nGenPart; i++) {
	    if( abs(GenPart_pdgId[i]) == 5 && GenPart_genPartIdxMother[i] == genHiggsIndex[j] ) {
	      TLorentzVector b;
	      b.SetPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	      genbVector.push_back(b);
	    }
	    if( abs(GenPart_pdgId[i]) == 24 && GenPart_genPartIdxMother[i] == genHiggsIndex[j] && GenPart_status[i]==22) {
	      TLorentzVector w;
	      w.SetPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	      genWVector.push_back(w);
	      genWIndex.push_back(i);
	    }
	  }
	  std::vector< TLorentzVector > genWdauVector;
	  std::vector< int > genWdecay, genWdauId;
	  for(int k = 0; k < genWIndex.size(); k++) {
	    int wdecay=0;
	    

	    bool found=true;
	    for(int i = 0; i < nGenPart; i++) {
	      if(abs(GenPart_pdgId[i]) == 24)
		continue;
	      // codes to label W decay
	      // 1: Wquarks, 2: Wenu, 3: Wmnu, 4: Wtaunu
	      if(GenPart_genPartIdxMother[i] == genWIndex[k]){

		// if a W daughter is a W keep looking further
		if(abs(GenPart_pdgId[i]) == 24) {
		  int windex_tmp = i;
		  for(int j = 0; j < nGenPart; j++) {
		    if(GenPart_genPartIdxMother[j] == windex_tmp){
		      if(abs(GenPart_pdgId[j]) <= 5) wdecay=1;
		      if(abs(GenPart_pdgId[j]) == 11 || abs(GenPart_pdgId[j]) == 12) wdecay=2;
		      if(abs(GenPart_pdgId[j]) == 13 || abs(GenPart_pdgId[j]) == 14) wdecay=3;
		      if(abs(GenPart_pdgId[j]) == 15 || abs(GenPart_pdgId[j]) == 16) wdecay=4;
		      if(abs(GenPart_pdgId[j]) == 22) wdecay=5;
		      if(found){
			genWdecay.push_back(wdecay);
			found=false;
		      }
		      TLorentzVector d;
		      d.SetPtEtaPhiM( GenPart_pt[j], GenPart_eta[j], GenPart_phi[j], GenPart_mass[j]);
		      genWdauVector.push_back(d);
		      genWdauId.push_back( abs(GenPart_pdgId[j]) );
		    }
		  }
		}
		else{
		  if(abs(GenPart_pdgId[i]) <= 5) wdecay=1;
		  if(abs(GenPart_pdgId[i]) == 11 || abs(GenPart_pdgId[i]) == 12) wdecay=2; 
		  if(abs(GenPart_pdgId[i]) == 13 || abs(GenPart_pdgId[i]) == 14) wdecay=3;
		  if(abs(GenPart_pdgId[i]) == 15 || abs(GenPart_pdgId[i]) == 16) wdecay=4;
		  if(abs(GenPart_pdgId[i]) == 22) wdecay=5;
		  if(found){
		    genWdecay.push_back(wdecay);
		    found=false;
		  }
		  TLorentzVector d;
		  d.SetPtEtaPhiM( GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
		  genWdauVector.push_back(d);
		  genWdauId.push_back( abs(GenPart_pdgId[i]) );
		}
		//break;
	      }
	    }
	  }
	  if(j==0){
	    genHiggs1Pt = genHiggsVector[j].Pt();
	    genHiggs1Eta = genHiggsVector[j].Eta();
	    genHiggs1Phi = genHiggsVector[j].Phi();
	    if(genbVector.size()>1){
	      genHiggs1b1Pt = genbVector[0].Pt();
              genHiggs1b1Eta = genbVector[0].Eta();
              genHiggs1b1Phi = genbVector[0].Phi();
	      genHiggs1b2Pt = genbVector[1].Pt();
              genHiggs1b2Eta = genbVector[1].Eta();
              genHiggs1b2Phi = genbVector[1].Phi();
	    }
	    if(genWVector.size()>1){
	      genHiggs1W1Pt = genWVector[0].Pt();
	      genHiggs1W1Eta = genWVector[0].Eta();
	      genHiggs1W1Phi = genWVector[0].Phi();
              genHiggs1W1M = genWVector[0].M();
              genHiggs1W1Decay = genWdecay[0];
	      genHiggs1W2Pt = genWVector[1].Pt();
	      genHiggs1W2Eta = genWVector[1].Eta();
	      genHiggs1W2Phi = genWVector[1].Phi();
              genHiggs1W2M = genWVector[1].M();
              genHiggs1W2Decay = genWdecay[1];
	      if(genWdauVector.size()>3){
		genHiggs1W1dau1Pt = genWdauVector[0].Pt();
		genHiggs1W1dau1Eta = genWdauVector[0].Eta();
		genHiggs1W1dau1Phi = genWdauVector[0].Phi();
		genHiggs1W1dau1Id = genWdauId[0];
		genHiggs1W1dau2Pt = genWdauVector[1].Pt();
		genHiggs1W1dau2Eta = genWdauVector[1].Eta();
		genHiggs1W1dau2Phi = genWdauVector[1].Phi();
		genHiggs1W1dau2Id = genWdauId[1];
		genHiggs1W2dau1Pt = genWdauVector[2].Pt();
		genHiggs1W2dau1Eta = genWdauVector[2].Eta();
		genHiggs1W2dau1Phi = genWdauVector[2].Phi();
		genHiggs1W2dau1Id = genWdauId[2];
		genHiggs1W2dau2Pt = genWdauVector[3].Pt();
		genHiggs1W2dau2Eta = genWdauVector[3].Eta();
		genHiggs1W2dau2Phi = genWdauVector[3].Phi();
                genHiggs1W2dau2Id = genWdauId[3];
	      }
	    }
	  }
          if(j==1){
            genHiggs2Pt = genHiggsVector[j].Pt();
            genHiggs2Eta = genHiggsVector[j].Eta();
            genHiggs2Phi = genHiggsVector[j].Phi();
	    if(genbVector.size()>1){
	      genHiggs2b1Pt = genbVector[0].Pt();
              genHiggs2b1Eta = genbVector[0].Eta();
              genHiggs2b1Phi = genbVector[0].Phi();
              genHiggs2b2Pt = genbVector[1].Pt();
              genHiggs2b2Eta = genbVector[1].Eta();
              genHiggs2b2Phi = genbVector[1].Phi();
            }
	    if(genWVector.size()>1){
              genHiggs2W1Pt = genWVector[0].Pt();
              genHiggs2W1Eta = genWVector[0].Eta();
              genHiggs2W1Phi = genWVector[0].Phi();
              genHiggs2W1M = genWVector[0].M();
              genHiggs2W1Decay = genWdecay[0];
              genHiggs2W2Pt = genWVector[1].Pt();
              genHiggs2W2Eta = genWVector[1].Eta();
              genHiggs2W2Phi = genWVector[1].Phi();
              genHiggs2W2M = genWVector[1].M();
              genHiggs2W2Decay = genWdecay[1];
              if(genWdauVector.size()>3){
                genHiggs2W1dau1Pt = genWdauVector[0].Pt();
                genHiggs2W1dau1Eta = genWdauVector[0].Eta();
                genHiggs2W1dau1Phi = genWdauVector[0].Phi();
		genHiggs2W1dau1Id = genWdauId[0];
                genHiggs2W1dau2Pt = genWdauVector[1].Pt();
                genHiggs2W1dau2Eta = genWdauVector[1].Eta();
                genHiggs2W1dau2Phi = genWdauVector[1].Phi();
                genHiggs2W1dau2Id = genWdauId[1];
                genHiggs2W2dau1Pt = genWdauVector[2].Pt();
                genHiggs2W2dau1Eta = genWdauVector[2].Eta();
                genHiggs2W2dau1Phi = genWdauVector[2].Phi();
                genHiggs2W2dau1Id = genWdauId[2];
                genHiggs2W2dau2Pt = genWdauVector[3].Pt();
                genHiggs2W2dau2Eta = genWdauVector[3].Eta();
                genHiggs2W2dau2Phi = genWdauVector[3].Phi();
                genHiggs2W2dau2Id = genWdauId[3];
              }
            }
          }
	}

	//gen level HH
	if(genHiggsVector.size() > 1) { 
	  genHH_pt = (genHiggsVector[0]+genHiggsVector[1]).Pt();
	  genHH_eta = (genHiggsVector[0]+genHiggsVector[1]).Eta();
	  genHH_phi = (genHiggsVector[0]+genHiggsVector[1]).Phi();
	  genHH_mass= (genHiggsVector[0]+genHiggsVector[1]).M();
	}
      
      }//end if !data

      //------------------------------
      //-------find fatJet------------
      //------------------------------
      vector<int> selectedFatJetIndices;

      for(unsigned int i = 0; i < nFatJet; i++ ) {       
	//Hbb fat jet pre-selection
	if (FatJet_pt[i] < 200) continue;
	selectedFatJetIndices.push_back(i);
      }

      // this should change at some point on a selection on the tagger
      int fatJet1Index = -1;
      int fatJet2Index = -1;
      int fatJet3Index = -1;
      for(unsigned int i = 0; i < selectedFatJetIndices.size(); i++ ) {
	if(i==0){
	  fatJet1Index = selectedFatJetIndices[i];
	}
	if(i==1){
	  fatJet2Index = selectedFatJetIndices[i];
	}
	if(i==2){
	  fatJet3Index = selectedFatJetIndices[i];
	}
      }

      //------------------------------------------------------
      //----------Fill higgs candidate information
      //------------------------------------------------------   
      TLorentzVector Higgs1Jet;
      Higgs1Jet.SetPtEtaPhiM(FatJet_pt[fatJet1Index],FatJet_eta[fatJet1Index],FatJet_phi[fatJet1Index],FatJet_msoftdrop[fatJet1Index]);
      float Higgs1MinDR = 999.;
      int Higgs1_match_idx = -1;
      for( int j = 0; j < genHiggsVector.size(); j++) {
	if(Higgs1Jet.DeltaR(genHiggsVector[j]) < Higgs1MinDR) {
	  Higgs1MinDR = Higgs1Jet.DeltaR(genHiggsVector[j]);
	  Higgs1_match_idx = j;
	}
      }
      fatJet1Pt = FatJet_pt[fatJet1Index];
      fatJet1Eta = FatJet_eta[fatJet1Index];
      fatJet1Phi = FatJet_phi[fatJet1Index];
      fatJet1Mass = FatJet_mass[fatJet1Index];
      fatJet1MassSD = FatJet_msoftdrop[fatJet1Index];
      fatJet1DeepAK8_H = FatJet_deepTag_H[fatJet1Index];
      fatJet1DeepAK8MD_H4qvsQCD = FatJet_deepTagMD_H4qvsQCD[fatJet1Index];
      fatJet1PNetXbb = FatJet_ParticleNetMD_probXbb[fatJet1Index]/(1.0 - FatJet_ParticleNetMD_probXcc[fatJet1Index] - FatJet_ParticleNetMD_probXqq[fatJet1Index]);
      fatJet1PNetXbb_alt = FatJet_particleNetMD_Xbb[fatJet1Index]/(1.0 - FatJet_particleNetMD_Xcc[fatJet1Index] - FatJet_particleNetMD_Xqq[fatJet1Index]);
      fatJet1PNetQCDb = FatJet_ParticleNetMD_probQCDb[fatJet1Index];
      fatJet1PNetQCDbb = FatJet_ParticleNetMD_probQCDbb[fatJet1Index];
      fatJet1PNetQCDc = FatJet_ParticleNetMD_probQCDc[fatJet1Index];
      fatJet1PNetQCDcc = FatJet_ParticleNetMD_probQCDcc[fatJet1Index];
      fatJet1PNetQCDothers = FatJet_ParticleNetMD_probQCDothers[fatJet1Index];
      if(Higgs1MinDR < 0.4) {
	fatJet1GenMatchIndex = Higgs1_match_idx;
      }
      fatJet1Tau4OverTau3 = FatJet_tau4[fatJet1Index] / FatJet_tau3[fatJet1Index];
      fatJet1Tau4OverTau2 = FatJet_tau4[fatJet1Index] / FatJet_tau2[fatJet1Index];
      fatJet1Tau4OverTau1 = FatJet_tau4[fatJet1Index] / FatJet_tau1[fatJet1Index];
      fatJet1Tau2OverTau1 = FatJet_tau2[fatJet1Index] / FatJet_tau1[fatJet1Index];
      fatJet1n2b1 = FatJet_n2b1[fatJet1Index];
      fatJet1lsf3 = FatJet_lsf3[fatJet1Index];
      //find any bjets in opposite hemisphere as fatJet1
      double btagMediumCut = -1;
      if (year == "2016") btagMediumCut = 0.3033;
      else if (year == "2017") btagMediumCut = 0.3033 ;
      else if (year == "2018") btagMediumCut = 0.2770;
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_pt[q] > 25 && Jet_btagDeepB[q] > btagMediumCut && 
	    deltaPhi(fatJet1Phi, Jet_phi[q]) > 2.5
	    ) {
	  fatJet1OppositeHemisphereHasBJet = true;
	  break;
	}
      }
	

      //find muon inside jet
      for(unsigned int q = 0; q < nMuon; q++ ) {       
	if (Muon_pt[q] > 30 && Muon_looseId[q] && 
	    deltaR(fatJet1Eta , fatJet1Phi, Muon_eta[q], Muon_phi[q]) < 1.0
	    ) {
	  fatJet1HasMuon = true;
	  break;
	}
      }
      //find electron inside jet
      for(unsigned int q = 0; q < nElectron; q++ ) {       
	if (Electron_pt[q] > 30 && Electron_mvaFall17V2noIso_WP90[q] && 
	    deltaR(fatJet1Eta , fatJet1Phi, Electron_eta[q], Electron_phi[q]) < 1.0
	    ) {
	  fatJet1HasElectron = true;
	  break;
	}
      }
      //find loose b-tagged jet inside jet
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_btagDeepB[q] > 0.0521 && 
	    deltaR(fatJet1Eta , fatJet1Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	  fatJet1HasBJetCSVLoose = true;
	  break;
	}
      }
     //find medium b-tagged jet inside jet
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_btagDeepB[q] > 0.3033 && 
	    deltaR(fatJet1Eta , fatJet1Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	  fatJet1HasBJetCSVMedium = true;
	  break;
	}
      }
      //find tight b-tagged jet inside jet
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_btagDeepB[q] > 0.7489 && 
	    deltaR(fatJet1Eta , fatJet1Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	  fatJet1HasBJetCSVTight = true;
	  break;
	}
      }


      TLorentzVector Higgs2Jet;
      Higgs2Jet.SetPtEtaPhiM(FatJet_pt[fatJet2Index],FatJet_eta[fatJet2Index],FatJet_phi[fatJet2Index],FatJet_msoftdrop[fatJet2Index]);
      float Higgs2MinDR = 999.;
      int Higgs2_match_idx = -1;
      for( int j = 0; j < genHiggsVector.size(); j++) {
	if(Higgs2Jet.DeltaR(genHiggsVector[j]) < Higgs2MinDR) {
	  Higgs2MinDR = Higgs2Jet.DeltaR(genHiggsVector[j]);
	  Higgs2_match_idx = j;
	}
      }
     
      fatJet2Pt = FatJet_pt[fatJet2Index];
      fatJet2Eta = FatJet_eta[fatJet2Index];
      fatJet2Phi = FatJet_phi[fatJet2Index];
      fatJet2Mass = FatJet_mass[fatJet2Index];
      fatJet2MassSD = FatJet_msoftdrop[fatJet2Index];
      fatJet2DeepAK8_H = FatJet_deepTag_H[fatJet2Index];
      fatJet2DeepAK8MD_H4qvsQCD = FatJet_deepTagMD_H4qvsQCD[fatJet2Index];
      fatJet2PNetXbb = FatJet_ParticleNetMD_probXbb[fatJet2Index]/(1.0 - FatJet_ParticleNetMD_probXcc[fatJet2Index] - FatJet_ParticleNetMD_probXqq[fatJet2Index]);
      fatJet2PNetXbb_alt = FatJet_particleNetMD_Xbb[fatJet2Index]/(1.0 - FatJet_particleNetMD_Xcc[fatJet2Index] - FatJet_particleNetMD_Xqq[fatJet2Index]);
      fatJet2PNetQCDb = FatJet_ParticleNetMD_probQCDb[fatJet2Index];
      fatJet2PNetQCDbb = FatJet_ParticleNetMD_probQCDbb[fatJet2Index];
      fatJet2PNetQCDc = FatJet_ParticleNetMD_probQCDc[fatJet2Index];
      fatJet2PNetQCDcc = FatJet_ParticleNetMD_probQCDcc[fatJet2Index];
      fatJet2PNetQCDothers = FatJet_ParticleNetMD_probQCDothers[fatJet2Index];
      if(Higgs2MinDR < 0.4) {
	fatJet2GenMatchIndex = Higgs2_match_idx;
      }
      fatJet2Tau4OverTau3 = FatJet_tau4[fatJet2Index] / FatJet_tau3[fatJet2Index];
      fatJet2Tau4OverTau2 = FatJet_tau4[fatJet2Index] / FatJet_tau2[fatJet2Index];
      fatJet2Tau4OverTau1 = FatJet_tau4[fatJet2Index] / FatJet_tau1[fatJet2Index];
      fatJet2Tau2OverTau1 = FatJet_tau2[fatJet2Index] / FatJet_tau1[fatJet2Index];
      fatJet2n2b1 = FatJet_n2b1[fatJet2Index];
      fatJet2lsf3 = FatJet_lsf3[fatJet2Index];      
      //find muon inside jet
      for(unsigned int q = 0; q < nMuon; q++ ) {       
	if (Muon_pt[q] > 30 && Muon_looseId[q] && 
	    deltaR(fatJet2Eta , fatJet2Phi, Muon_eta[q], Muon_phi[q]) < 1.0
	    ) {
	  fatJet2HasMuon = true;
	  break;
	}
      }
      //find electron inside jet
      for(unsigned int q = 0; q < nElectron; q++ ) {       
	if (Electron_pt[q] > 30 && Electron_mvaFall17V2noIso_WP90[q] && 
	    deltaR(fatJet2Eta , fatJet2Phi, Electron_eta[q], Electron_phi[q]) < 1.0
	    ) {
	  fatJet2HasElectron = true;
	  break;
	}
      }
      //find loose b-tagged jet inside jet
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_btagDeepB[q] > 0.0521 && 
	    deltaR(fatJet2Eta , fatJet2Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	  fatJet2HasBJetCSVLoose = true;
	  break;
	}
      }
      //find medium b-tagged jet inside jet
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_btagDeepB[q] > 0.3033 && 
	    deltaR(fatJet2Eta , fatJet2Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	  fatJet2HasBJetCSVMedium = true;
	  break;
	}
      }
      //find tight b-tagged jet inside jet
      for(unsigned int q = 0; q < nJet; q++ ) {       
	if (Jet_btagDeepB[q] > 0.7489 && 
	    deltaR(fatJet2Eta , fatJet2Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	  fatJet2HasBJetCSVTight = true;
	  break;
	}
      }


      //------------------------------------------------------
      //----------Fill Jet 3 information
      //------------------------------------------------------
      if(fatJet3Index!=-1) {
	fatJet3Pt = FatJet_pt[fatJet3Index];
	fatJet3Eta = FatJet_eta[fatJet3Index];
	fatJet3Phi = FatJet_phi[fatJet3Index];
	fatJet3Mass = FatJet_mass[fatJet3Index];
	fatJet3MassSD = FatJet_msoftdrop[fatJet3Index];
	fatJet3DeepAK8_H = FatJet_deepTag_H[fatJet3Index];
	fatJet3DeepAK8MD_H4qvsQCD = FatJet_deepTagMD_H4qvsQCD[fatJet3Index];
	fatJet3PNetXbb = FatJet_ParticleNetMD_probXbb[fatJet3Index]/(1.0 - FatJet_ParticleNetMD_probXcc[fatJet3Index] - FatJet_ParticleNetMD_probXqq[fatJet3Index]);
	fatJet3PNetXbb_alt = FatJet_particleNetMD_Xbb[fatJet3Index]/(1.0 - FatJet_particleNetMD_Xcc[fatJet3Index] - FatJet_particleNetMD_Xqq[fatJet3Index]);
	fatJet3PNetQCDb = FatJet_ParticleNetMD_probQCDb[fatJet3Index];
	fatJet3PNetQCDbb = FatJet_ParticleNetMD_probQCDbb[fatJet3Index];
	fatJet3PNetQCDc = FatJet_ParticleNetMD_probQCDc[fatJet3Index];
	fatJet3PNetQCDcc = FatJet_ParticleNetMD_probQCDcc[fatJet3Index];
	fatJet3PNetQCDothers = FatJet_ParticleNetMD_probQCDothers[fatJet3Index];
	fatJet3Tau4OverTau3 = FatJet_tau4[fatJet3Index] / FatJet_tau3[fatJet3Index];
	fatJet3Tau4OverTau2 = FatJet_tau4[fatJet3Index] / FatJet_tau2[fatJet3Index];
	fatJet3Tau4OverTau1 = FatJet_tau4[fatJet3Index] / FatJet_tau1[fatJet3Index];
	fatJet3Tau2OverTau1 = FatJet_tau2[fatJet3Index] / FatJet_tau1[fatJet3Index];	

	//find muon inside jet
	for(unsigned int q = 0; q < nMuon; q++ ) {       
	  if (Muon_pt[q] > 30 && Muon_looseId[q] && 
	      deltaR(fatJet3Eta , fatJet3Phi, Muon_eta[q], Muon_phi[q]) < 1.0
	      ) {
	    fatJet3HasMuon = true;
	    break;
	  }
	}
	//find electron inside jet
	for(unsigned int q = 0; q < nElectron; q++ ) {       
	  if (Electron_pt[q] > 30 && Electron_mvaFall17V2noIso_WP90[q] && 
	      deltaR(fatJet3Eta , fatJet3Phi, Electron_eta[q], Electron_phi[q]) < 1.0
	    ) {
	    fatJet3HasElectron = true;
	    break;
	  }
	}
	//find loose b-tagged jet inside jet
	for(unsigned int q = 0; q < nJet; q++ ) {       
	  if (Jet_btagDeepB[q] > 0.0521 && 
	      deltaR(fatJet3Eta , fatJet3Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	    ) {
	    fatJet3HasBJetCSVLoose = true;
	    break;
	  }
	}
	//find medium b-tagged jet inside jet
	for(unsigned int q = 0; q < nJet; q++ ) {       
	  if (Jet_btagDeepB[q] > 0.3033 && 
	    deltaR(fatJet3Eta , fatJet3Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	      ) {
	    fatJet3HasBJetCSVMedium = true;
	    break;
	  }
	}
	//find tight b-tagged jet inside jet
	for(unsigned int q = 0; q < nJet; q++ ) {       
	  if (Jet_btagDeepB[q] > 0.7489 && 
	      deltaR(fatJet3Eta , fatJet3Phi, Jet_eta[q], Jet_phi[q]) < 1.0
	      ) {
	    fatJet3HasBJetCSVTight = true;
	  break;
	  }
	}
      }

      //------------------------------------------------------
      //----------Fill hh candidate information
      //------------------------------------------------------
      hh_pt = (Higgs1Jet+Higgs2Jet).Pt();
      hh_eta = (Higgs1Jet+Higgs2Jet).Eta();
      hh_phi = (Higgs1Jet+Higgs2Jet).Phi();
      hh_mass = (Higgs1Jet+Higgs2Jet).M();      
    
      fatJet1PtOverMHH = fatJet1Pt / hh_mass;
      fatJet1PtOverMSD = fatJet1Pt / fatJet1MassSD;
      fatJet2PtOverMHH = fatJet2Pt / hh_mass;
      fatJet2PtOverMSD = fatJet2Pt / fatJet1MassSD;
      deltaEta_j1j2 = fabs(fatJet1Eta - fatJet2Eta);
      deltaPhi_j1j2 = deltaPhi(fatJet1Phi, fatJet2Phi);
      deltaR_j1j2 = deltaR(fatJet1Eta, fatJet1Phi, fatJet2Eta, fatJet2Phi);
      ptj2_over_ptj1 = fatJet2Pt / fatJet1Pt;
      mj2_over_mj1 = fatJet2MassSD / fatJet1MassSD;             
      
      // AK15 info
      vector<int> ak15selectedFatJetIndices;
      for(unsigned int i = 0; i < nFatJetAK15; i++ ) {
        if (FatJetAK15_pt[i] < 200) continue;
        ak15selectedFatJetIndices.push_back(i);
      }
      int ak15fatJet1Index = -1;
      int ak15fatJet2Index = -1;
      for(unsigned int i = 0; i < ak15selectedFatJetIndices.size(); i++ ) {
        if(i==0){
          ak15fatJet1Index = ak15selectedFatJetIndices[i];
        }
        if(i==1){
          ak15fatJet2Index = ak15selectedFatJetIndices[i];
        }
      }

      if(addAK15){
	if(ak15fatJet1Index>-1){
	  ak15fatJet1Pt = FatJetAK15_pt[ak15fatJet1Index];
	  ak15fatJet1Eta = FatJetAK15_eta[ak15fatJet1Index];
	  ak15fatJet1Phi = FatJetAK15_phi[ak15fatJet1Index];
	  ak15fatJet1Mass = FatJetAK15_mass[ak15fatJet1Index];
	  ak15fatJet1MassSD = FatJetAK15_msoftdrop[ak15fatJet1Index];
	  ak15fatJet1PNetMDXbb = FatJetAK15_ParticleNetMD_probXbb[ak15fatJet1Index]/(1-FatJetAK15_ParticleNetMD_probXcc[ak15fatJet1Index]-FatJetAK15_ParticleNetMD_probXqq[ak15fatJet1Index]);
	  ak15fatJet1PNetHqqqq = FatJetAK15_ParticleNet_probHqqqq[ak15fatJet1Index];
	  ak15fatJet1PNetQCDb = FatJetAK15_ParticleNet_probQCDb[ak15fatJet1Index];
	  ak15fatJet1PNetQCDbb = FatJetAK15_ParticleNet_probQCDbb[ak15fatJet1Index];
	  ak15fatJet1PNetQCDc = FatJetAK15_ParticleNet_probQCDc[ak15fatJet1Index];
	  ak15fatJet1PNetQCDcc = FatJetAK15_ParticleNet_probQCDcc[ak15fatJet1Index];
	  ak15fatJet1PNetQCDothers = FatJetAK15_ParticleNet_probQCDothers[ak15fatJet1Index];
	}
	if(ak15fatJet2Index>-1){
	  ak15fatJet2Pt = FatJetAK15_pt[ak15fatJet2Index];
          ak15fatJet2Eta = FatJetAK15_eta[ak15fatJet2Index];
          ak15fatJet2Phi = FatJetAK15_phi[ak15fatJet2Index];
          ak15fatJet2Mass = FatJetAK15_mass[ak15fatJet2Index];
          ak15fatJet2MassSD = FatJetAK15_msoftdrop[ak15fatJet2Index];
          ak15fatJet2PNetMDXbb = FatJetAK15_ParticleNetMD_probXbb[ak15fatJet2Index]/(1-FatJetAK15_ParticleNetMD_probXcc[ak15fatJet2Index]-FatJetAK15_ParticleNetMD_probXqq[ak15fatJet2Index]);
          ak15fatJet2PNetHqqqq = FatJetAK15_ParticleNet_probHqqqq[ak15fatJet2Index];
          ak15fatJet2PNetQCDb = FatJetAK15_ParticleNet_probQCDb[ak15fatJet2Index];
          ak15fatJet2PNetQCDbb = FatJetAK15_ParticleNet_probQCDbb[ak15fatJet2Index];
          ak15fatJet2PNetQCDc = FatJetAK15_ParticleNet_probQCDc[ak15fatJet2Index];
          ak15fatJet2PNetQCDcc = FatJetAK15_ParticleNet_probQCDcc[ak15fatJet2Index];
          ak15fatJet2PNetQCDothers = FatJetAK15_ParticleNet_probQCDothers[ak15fatJet2Index];
	}
      }
      //------------------------------------------------------
      //----------Find Leptons
      //------------------------------------------------------     
      for(unsigned int i = 0; i < nMuon; i++ ) {       

	if (Muon_pt[i] < 30) continue;
	if (fabs(Muon_eta[i]) > 2.4) continue;
	if (Muon_miniPFRelIso_all[i] > 0.2) continue;
	if (!Muon_tightId) continue;
	
	if (lep1Id == 0) {
	  lep1Pt = Muon_pt[i];
	  lep1Eta = Muon_eta[i];
	  lep1Phi = Muon_phi[i];
	  lep1Id = Muon_charge[i] * (13);
	} else if (Muon_pt[i] > lep1Pt) {
	  lep2Pt = lep1Pt;
	  lep2Eta = lep1Eta;
	  lep2Phi = lep1Phi;
	  lep2Id = lep1Id;
	  lep1Pt = Muon_pt[i];
	  lep1Eta = Muon_eta[i];
	  lep1Phi = Muon_phi[i];
	  lep1Id = Muon_charge[i] * (13);
	} else if (lep2Id == 0 || Muon_pt[i] > lep2Pt) {
	  lep2Pt = Muon_pt[i];
	  lep2Eta = Muon_eta[i];
	  lep2Phi = Muon_phi[i];
	  lep2Id = Muon_charge[i] * (13);
	} 
      } //loop over muons

      for(unsigned int i = 0; i < nElectron; i++ ) {       

	if (Electron_pt[i] < 35) continue;
	if (fabs(Electron_eta[i]) > 2.5) continue;
	if (Electron_miniPFRelIso_all[i] > 0.2) continue;
	if (!Electron_cutBased[i]) continue;
	
	if (lep1Id == 0) {
	  lep1Pt = Electron_pt[i];
	  lep1Eta = Electron_eta[i];
	  lep1Phi = Electron_phi[i];
	  lep1Id = Electron_charge[i] * (11);
	} else if (Electron_pt[i] > lep1Pt) {
	  lep2Pt = lep1Pt;
	  lep2Eta = lep1Eta;
	  lep2Phi = lep1Phi;
	  lep2Id = lep1Id;
	  lep1Pt = Electron_pt[i];
	  lep1Eta = Electron_eta[i];
	  lep1Phi = Electron_phi[i];
	  lep1Id = Electron_charge[i] * (11);
	} else if (lep2Id == 0 || Electron_pt[i] > lep2Pt) {
	  lep2Pt = Electron_pt[i];
	  lep2Eta = Electron_eta[i];
	  lep2Phi = Electron_phi[i];
	  lep2Id = Electron_charge[i] * (11);
	} 
      } //loop over electrons


      //------------------------------------------------------
      //----------Find Photons
      //------------------------------------------------------     
      for(unsigned int i = 0; i < nPhoton; i++ ) {       

	if (Photon_pt[i] < 30) continue;
	if (fabs(Photon_eta[i]) > 2.5) continue;
	if (fabs(Photon_eta[i]) < 1.5) {
	  if (Photon_mvaID[i] > 0.42) continue;
	} else {
	  if (Photon_mvaID[i] > 0.14) continue;
	}
	if (!Photon_electronVeto[i]) continue;
	
	pho1Pt = Photon_pt[i];
	pho1Eta = Photon_eta[i];
	pho1Phi = Photon_phi[i];

      }

      //*******************************
      //Count additional AK4 jets 
      //*******************************
      for(int i = 0; i < nJet; i++) {
	if (Jet_pt[i] > 30 && fabs(Jet_eta[i]) < 2.5
	    && deltaR(Jet_eta[i] , Jet_phi[i], fatJet1Eta, fatJet1Phi) > 0.8
	    && deltaR(Jet_eta[i] , Jet_phi[i], fatJet2Eta, fatJet2Phi) > 0.8
	    ) {
	  NJets++;
	}

	bool passBTag = false;
	if (year == "2018") {
	  if (Jet_pt[i] > 40 && fabs(Jet_eta[i]) < 2.5
	      && Jet_btagDeepFlavB[i] > 0.2770
	      && Jet_puId[i] >= 2	
	      && Jet_jetId[i] >= 4
	      ) {
	    passBTag = true;
	  }
	} else if (year == "2017" ) {
	  if (Jet_pt[i] > 40 && fabs(Jet_eta[i]) < 2.5
	      && Jet_btagDeepFlavB[i] > 0.3033
	      && Jet_puId[i] >= 2
	      && Jet_jetId[i] >= 4
	      ) {
	    passBTag = true;
	  }
	} else if (year == "2016") {
	  if (Jet_pt[i] > 30 && fabs(Jet_eta[i]) < 2.4
	      && Jet_btagDeepFlavB[i] > 0.3033
	      && Jet_puId[i] >= 2
	      && Jet_jetId[i] >= 4
	      ) {
	    passBTag = true;
	  }
	}
	
	if (passBTag) {
	  nBTaggedJets++;

	  if (nBTaggedJets==1) {
	    jet1Pt = Jet_pt[i];
	    jet1Eta = Jet_eta[i];
	    jet1Phi = Jet_phi[i];
	  }
	  if (nBTaggedJets==2) {
	    jet2Pt = Jet_pt[i];
	    jet2Eta = Jet_eta[i];
	    jet2Phi = Jet_phi[i];
	  }
	  if (nBTaggedJets==3) {
	    jet3Pt = Jet_pt[i];
	    jet3Eta = Jet_eta[i];
	    jet3Phi = Jet_phi[i];
	  }
	  if (nBTaggedJets==4) {
	    jet4Pt = Jet_pt[i];
	    jet4Eta = Jet_eta[i];
	    jet4Phi = Jet_phi[i];
	  }


	}

      } //loop over AK4 jets
      
       
      //****************************************************
      //Fill Event - skim for events with two jets found
      //****************************************************
      if (
	  Option == 0 || 
	  (Option == 1 && fatJet1Pt > 250 && fatJet2Pt > 250 && fatJet1MassSD > 20 
	   && fatJet2MassSD > 20) 
	  ) {
	 

	//****************************************************
	//Compute pileupWeight
	//****************************************************      
	if (pileupWeightHist) {
	  pileupWeight = pileupWeightHist->GetBinContent( pileupWeightHist->GetXaxis()->FindFixBin(Pileup_nTrueInt));
	}

	//****************************************************
	//Compute totalWeight
	//****************************************************      
	totalWeight = weight * pileupWeight;
	
	
	NEventsFilled++;            
	outputTree->Fill();      
      }
    }//end of event loop

    cout << "Filled Total of " << NEventsFilled << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}
