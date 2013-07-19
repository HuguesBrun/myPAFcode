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
    
    
  enum RunRange {RunAB, RunC, RunD};
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
    TH1F* hMassUncut[9];
    
    
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
    TTree *treeClusterShape;
    
    
  virtual void Initialise();
  virtual void InsideLoop();
  virtual void Summary();  
  virtual int  isGenMatched(int, int, int);
  virtual void compWeights(int);
  virtual void fillTheHisto(int);
  virtual void fillTheHistoMass(int , float );
  virtual void fillTheTestTree(int, int);
  virtual bool passTightIdNoIso(int, float, float, float, float, float, float, float, int, int);
  virtual bool passTightIdNoIsoModified(int, float, float, float, float, float, float, float, int, int);
  virtual float giveThePOGiso(float, float, float, float, float, float);
  virtual float giveThePOGisoPU2012(float, float, float, float, float, float);
  virtual float giveMyPOGiso(float, TLorentzVector *,float, bool, float, float);
    virtual float calc03Iso(float,  float, float,int);
  virtual float calcPFIso(TLorentzVector*, float, PFisolationType, bool);
  virtual float calcPFRadIso(TLorentzVector *,float , float , PFisolationType , bool , bool , bool );
  virtual float giveRadIso(TLorentzVector *, float, bool , bool , bool );
  virtual float calcPFRadIsoFonc(TLorentzVector *,float, float, PFisolationType, bool, bool, bool, int);
  virtual float giveRadIsoFonc(TLorentzVector *, float, bool, bool, bool, int);

  virtual float PFisolationWithDeltaBeta(int);
  virtual void FillThePFtree(TLorentzVector *, float , bool , bool);
  virtual float congPhi(float );
  virtual bool  passPreCuts(float, int, float, float, float,float);
  virtual bool  passIPcuts(float, float);
    virtual bool  FOnoDeta(float , int , float ,float ,float ,float , float ,bool , int , float , float , float );
    virtual bool  FO_full(float, int, float,float, float,float,float, float,bool, int, float, float, float);

    virtual  bool passMissItCons(bool , int );
    virtual  bool FOnoIso(float , int , float , float , float ,float ,float , float ,bool , int );
    virtual void fillTheCLTree(int );
    virtual float correctForNoise(float , bool , HWWAnalysis::RunRange , bool );
    virtual float correctForHLTDefinition(float iso, bool isBarrel, HWWAnalysis::RunRange runRange);


 // virtual bool isAGoodEvent(int, int);


    RunRange runRange;

    
    bool isMC;
    bool doMuons;
    bool runPF; 
    TString signal;
    int runNumber;
    int lumiSec;
    int eventNumber;
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
    float weight_runAsingle;
    float weight_runDsingle;
    float weight_runAB_Rdep;
    float weight_runC_Rdep;
    float weight_runD_Rdep;
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
    int tag_passISO;
    int pair_passISO;
    int passFO_ISO;
    int passFO_BDT_ISO;
    int passBDT_ISO;
    int passAllNoIsoDet;
    int passNM1IP;
    int passNM1presel;
    int passNM1convs;
    int trigSingle;
    int trigDoubleLeg0;
    int trigDoubleLeg1;
    int passPreselec;
    int passIP;
    int passConvs;
    int passFOnoIso;
    int isSameSign;
    int topSelection;
    int passFO_plus;
    int passFO_noDeta;
    int passFO_full;
    int passFO_full_HLT;
    int passFO_full_HLT_noise;
    int passFO_full_noise;
    int topFO;

    
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
                                     
    // treeClusterShape
    //general variables

    
     // MVA tracking 
    float CL_fBrem;
    float CL_hkfchi2;
    float CL_hkfhits;
    // Geometrical matiching
    float CL_deltaPhiIn;
    float CL_deltaEtaIn;
    float CL_detacalo;
    // ECAL shower variables
    float CL_see;
    float CL_spp;
    float CL_etawidth;
    float CL_phiwidth;
    float CL_e1x5e5x5;
    float CL_R9;
    
    // energy matching
    float CL_HoE;
    float CL_EoP;
    float CL_IoEmIoP;
    float CL_eleEoPout;
    float CL_preshowerOverRaw;

    // IP
    float CL_d0;
    float CL_ip3d;
    
    // more vars
    float CL_sigmaIetaIeta;
    float CL_passConversionVeto;
    float CL_isEcalDriven;
    float CL_hnHits;
    float CL_dZ;
    float CL_MVA;
    float CL_CombIsoHWW; 
    
    // the detIsos
    float CL_isoECAL;
    float CL_isoECAL_HLT;
    float CL_isoECAL_noise;
    float CL_isoECAL_HLT_noise;
    float CL_isoHCAL;
    float CL_isoTracker;

    // the detIsos
    float CL_isoECALRelat;
    float CL_isoECALRelatModif;
    float CL_isoHCALRelat;
    float CL_isoTrackerRelat;
    
    // the non trig MVA
    float CL_nonTrigMVA;
    
    // passPreselection
    float CL_passPreselection;
    
    // PF isolation
    float CL_PFchargedIso;
    float CL_POGcombIso;
    
    // relatIsoPF
    float CL_relatPFiso03;
    


  ClassDef(HWWAnalysis,0);
};
#endif

