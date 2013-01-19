#ifndef HWWANALYSIS_H
#define HWWANALYSIS_H 1

#include "packages/TCounterUI/TCounterUI.h"
#include "packages/CMSAnalysisSelectorMiniTrees/CMSAnalysisSelectorMiniTrees.h"


#include "TH1F.h"
#include "TH2F.h"
#include <vector>
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "GlobalVariables.h"

class TLorentzVector;

class HWWAnalysis: public CMSAnalysisSelectorMiniTrees {
 public:
  HWWAnalysis(TTree *tree=0);
  virtual ~HWWAnalysis() {}
  
  class lepton{
  public:
    lepton(){};
    lepton(TLorentzVector vec, int ch, int ty, int ind){
      p = vec;
      charge = ch;
      type = ty;
      index = ind;
    };
    TLorentzVector p;
    int charge;
    int type; // -1(unknown), 0(mu), 1(ele)
    int index;
  };
  enum WorkingPoint {
        VETO,
        LOOSE,
        MEDIUM,
        TIGHT
  };
    
 enum PFisolationType {
        neutralHadron,
        chargedHadron,
        gamma,
        //TIGHT
 };
 protected:

    TH1F *passCuts;

    TH1F* hkfchi2[9];
    TH1F* hkfhits[9];
    TH1F* hgsfchi2[9];
    TH1F* hdetacalo[9];
    TH1F* hsee[9];
    TH1F* hspp[9];
    TH1F* hetawidth[9];
    TH1F* hphiwidth[9];
    TH1F* he1x5e5x5[9];
    TH1F* hr9[9];
    TH1F* hEoP[9];
    TH1F* hIoEmIeP[9];
    TH1F* hEoPout[9];
    TH1F* hPreShowerOverRaw[9];
    TH1F* hd0[9];
    TH1F* hIP3D[9];
    TH1F* hdZ[9];
    TH1F* hConbIsoHWW[9];
    TH1F* hMVAid[9];
    TH1F* hMvaiso[9];
    TH1F* hRadiaIso[9];
    TH1F* hRadiaIsoVeto[9];
    TH1F* hRadiaIsoVetoMore[9];
    TH1F* hchargedHadronIso[9];
    TH1F* hneutralHadronIso[9];
    TH1F* hphotonIso[9];
    TH1F* hchargedHadronIso04[9];
    TH1F* hneutralHadronIso04[9];
    TH1F* hphotonIso04[9];
    TH1F* hpassConversionVeto[9];
    TH1F* hsigmaIetaIeta[9];
    TH1F* hdeltaPhiIn[9];
    TH1F* hdeltaEtaIn[9];
    TH1F* hisEcalDriven[9];
    TH1F* hnLost[9];
    TH1F* hnHits[9];
    TH1F* hfBrem[9];
    TH1F* hPt[9];
    TH1F* hEta[9];
    
    TH1F* hMass;
    
    TH2F *noBias[2];
    TH2F *eventHLT[2];
    TH2F *eventFindApair[2];
    TH2F *eventFindAtag[2];
    TH2F *eventTagMatchedWithHlt[2];
    TH2F *eventInMassWindows[2];
    TH2F *eventProbeOnly0[2];
    TH2F *eventProbeOnly1[2];
    TH2F *eventProbeOnly2[2];
    TH2F *eventBeforeSaving[2];

    
    TTree *theTree; 
    TTree *testTree;
    TTree *isoInfo;
    
    
  virtual void Initialise();
  virtual void InsideLoop();
  virtual void Summary();  
  virtual int  isGenMatched(int, int, int);
  virtual void compWeights(int);
  virtual void fillTheHisto(int);
  virtual void fillTheTestTree(int, int);
  virtual bool passTightIdNoIso(int, float, float, float, float, float, float, float, int, int);
  virtual bool passTightIdNoIsoModified(int, float, float, float, float, float, float, float, int, int);
  virtual float giveThePOGiso(float, float, float, float, float, float);
  virtual float giveThePOGisoPU2012(float, float, float, float, float, float);
  virtual float giveMyPOGiso(float, TLorentzVector *,float, bool, float, float);
  virtual float calcPFIso(TLorentzVector*, float, PFisolationType, bool);
  virtual float calcPFRadIso(TLorentzVector *,float , float , PFisolationType , bool , bool , bool );
  virtual float giveRadIso(TLorentzVector *, float, bool , bool , bool );
  virtual float calcPFRadIsoFonc(TLorentzVector *,float, float, PFisolationType, bool, bool, bool, int);
  virtual float giveRadIsoFonc(TLorentzVector *, float, bool, bool, bool, int);

  virtual float PFisolationWithDeltaBeta(int);
  virtual void FillThePFtree(TLorentzVector *, float , bool , bool);
  virtual float congPhi(float );


    
    bool isMC;
    TString signal;
    float mass;
    float pt;
    float eta;
    float SCeta;
    float absEta;
    float absSCeta;
    float weight;
    float weight_runA;
    float weight_runB;
    float weight_runC;
    float weight_runD;
    int isFSR;
    int tag_isFSR;
    float tag_pt;
    float tag_absSCEta;
    float tp_deltaR;
    int nVtx;
    int eventMatched;
    int eventMatched2;
    int eventMatched3;
    int passTight;
    int passLoose; 
    int passFO;
    int passBDT; 
    int passISO;
    int passFO_BDT;
    int passFO_ISO;
    int passFO_BDT_ISO;
    int isSameSign;
    
    //for testTree
    int isTriggering;
    int theType;
    int isPair;
    float MVAidTrig;
    float MVAiso;
    float radIso;
    float radIsoVeto;
    float radIsoVetoMore;
    float combIsoHWW;
    int   passMVA;
    int  PDGid;
    int passTightIdNoIsoPart;
    int passTightLooserIdNoIsoPart;
    float POGisolation;
    float POGisolation2012A;
    float myPOGisolation;
    float POGnoRho;
    float PFisoDeltaBeta;
    float radIsoStandard;
    float radIsoNoInner;
    float radisoNoThreshold;
    float radisoNoThresholdNoInner;
    
    float radisoFonc0;
    float radisoFonc1;
    float radisoFonc2;
    float radisoFonc3;
    
    float chRadisoFonc0;
    float chRadisoFonc1;
    float chRadisoFonc2;
    float chRadisoFonc3;
    
    float nhRadisoFonc0;
    float nhRadisoFonc1;
    float nhRadisoFonc2;
    float nhRadisoFonc3;
    
    float gRadisoFonc0;
    float gRadisoFonc1;
    float gRadisoFonc2;
    float gRadisoFonc3;
    
    
    
    float radIsoBeta;

    int myTightID;
    int myTightIDiso;
    int nbRecoVertex;
    
    
    float chIso0;
    float chIso1;
    float chIso2;
    float chIso3;
    float chIso4;
    
    float nhIso0;
    float nhIso1;
    float nhIso2;
    float nhIso3;
    float nhIso4;
    
    float gIso0;
    float gIso1;
    float gIso2;
    float gIso3;
    float gIso4;

    float gComb0;
    float gComb1;
    float gComb2;
    float gComb3;
    float gComb4;


    //iso Info
    float PFdeltaEta;
    float PFdeltaPhi;
    float PFet;
    int   PFtype;
    int   PFisPU;
    float mainEta;
    float mainSCEta;
    float mainPhi;
    float mainEnergy;
    float mainPt;
    int isSignalBg;
    float PFdeltaR;
                                     
    
  //virtual int isPassingPOGid();
 

  ClassDef(HWWAnalysis,0);
};
#endif

