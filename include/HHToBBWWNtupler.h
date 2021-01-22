#ifndef DEF_HHToBBWWNtupler
#define DEF_HHToBBWWNtupler

#include "EventAnalyzer.h"
#include "TH2F.h"

class HHToBBWWNtupler: public EventAnalyzer {
    public: 
        HHToBBWWNtupler(TTree *tree=0): EventAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label, string pileupWeightName);

};

#endif
