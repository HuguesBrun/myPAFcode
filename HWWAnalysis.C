///////////////////////////////////////////////////////////////
//   Analyzer for the Fake Rate Studies within HWW group               
//
//         AUTHOR:    Santiago Folgueras             
//////////////////////////////////////////////////////////////
//
#include "HWWAnalysis.h" 
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "GlobalVariables.h"
#include "MuonEffectiveArea.h"
#include "MuonEffectiveArea.h"
#if !defined(__CINT__)
ClassImp(HWWAnalysis);
#endif

#define DEBUG
using namespace std;

HWWAnalysis::HWWAnalysis(TTree* tree):
  CMSAnalysisSelectorMiniTrees(tree)  {} 


void HWWAnalysis::Initialise() {
	cout << "initialisation ...  " << endl;
    
    GetInputParameters()->TheNamedBool("isMC", isMC);
    signal = GetInputParameters()->TheNamedString("Signal");
    cout << "le MC est a " << isMC << endl;
    
    passCuts = CreateH1F("cut_flow", "", 8,-0.5,7.5);
    hMass    = CreateH1F("diElecMass", "", 100,50,150);

    
    for (int index=0; index<9; index++){
        hkfchi2[index] = CreateH1F(Form("hkfchi2_%i",index),"",60,0,6);
        hkfhits[index] = CreateH1F(Form("hkfhits_%i",index),"",19,-1,18);
        hgsfchi2[index] = CreateH1F(Form("hgsfchi2_%i",index),"",80,0,80);
        hdetacalo[index] = CreateH1F(Form("hdetacalo_%i",index),"",80,-0.3,0.3 );
        hsee[index] = CreateH1F(Form("hsee_%i",index),"",100,0,0.04);
        hspp[index] = CreateH1F(Form("hspp_%i",index),"",100,0,0.1);
        hetawidth[index] = CreateH1F(Form("hetawidth_%i",index),"",50,0,0.06 );
        hphiwidth[index] = CreateH1F(Form("hphiwidth_%i",index),"",50,0,0.15);
        he1x5e5x5[index] = CreateH1F(Form("he1x5e5x5_%i",index),"",120,-0.1,1);
        hr9[index] = CreateH1F(Form("hr9_%i",index),"",110,0,1.1);
        hEoP[index] = CreateH1F(Form("hEoP_%i",index),"",100,0,10);
        hIoEmIeP[index] = CreateH1F(Form("hIoEmIeP_%i",index),"",60,-0.01,0.05);
        hEoPout[index] = CreateH1F(Form("hEoPout_%i",index),"",90,0,30);
        hPreShowerOverRaw[index] = CreateH1F(Form("hPreShowerOverRaw_%i",index),"",100,0,0.2);
        hd0[index] = CreateH1F(Form("hd0_%i",index),"",100,-0.1,0.1);
        hIP3D[index] = CreateH1F(Form("hIP3D_%i",index),"",100,-0.1,0.1);
        hdZ[index] = CreateH1F(Form("hdZ_%i",index),"",100,-0.5,0.5);
        hConbIsoHWW[index] = CreateH1F(Form("hConbIsoHWW_%i",index),"",60,0,30);
        hMVAid[index] = CreateH1F(Form("hMVAid_%i",index),"",44,-1.1,1.1 );
        hMvaiso[index] = CreateH1F(Form("hMvaiso_%i",index),"",44,-1.1,1.1);
        hRadiaIso[index] = CreateH1F(Form("hRadiaIso_%i",index),"",50,0,2);
        hRadiaIsoVeto[index] = CreateH1F(Form("hRadiaIsoVeto_%i",index),"",50,0,2);
        hRadiaIsoVetoMore[index] = CreateH1F(Form("hRadiaIsoVetoMore_%i",index),"",50,0,2);
        hchargedHadronIso[index] = CreateH1F(Form("hchargedHadronIso_%i",index),"",120,0,8);
        hneutralHadronIso[index] = CreateH1F(Form("hneutralHadronIso_%i",index),"",90,0,8);
        hphotonIso[index] = CreateH1F(Form("hphotonIso_%i",index),"",100,0,10 );
        hchargedHadronIso04[index] = CreateH1F(Form("hchargedHadronIso04_%i",index),"",120,0,15);
        hneutralHadronIso04[index] = CreateH1F(Form("hneutralHadronIso04_%i",index),"",80,0,10);
        hphotonIso04[index] = CreateH1F(Form("hphotonIso04_%i",index),"",120,0,20);
        hpassConversionVeto[index] = CreateH1F(Form("hpassConversionVeto_%i",index),"",2,0,2);
        hsigmaIetaIeta[index] = CreateH1F(Form("hsigmaIetaIeta_%i",index),"",100,0,0.04);
        hdeltaPhiIn[index] = CreateH1F(Form("hdeltaPhiIn_%i",index),"",80,-0.2,0.2);
        hdeltaEtaIn[index] = CreateH1F(Form("hdeltaEtaIn_%i",index),"",80,-0.05,0.05);
        hisEcalDriven[index] = CreateH1F(Form("hisEcalDriven_%i",index),"",2,0,2);
        hnLost[index] = CreateH1F(Form("hnLost_%i",index),"",7,0,7);
        hnHits[index] = CreateH1F(Form("hnHits_%i",index),"",7,0,7);
        hfBrem[index] = CreateH1F(Form("hfBrem_%i",index),"",60,-1,5);
        hPt[index] = CreateH1F(Form("hPt_%i",index),"",200,0,200);
        hEta[index] = CreateH1F(Form("hEta_%i",index),"",100,-3.5,3.5);
        hMassUncut[index] =  CreateH1F(Form("hMassUncut_%i",index),"",100,40,120);
    }
    
    double binX[6] = {0, 0.8, 1.4442, 1.556, 2, 2.5};
    double binY[7] = {10,15,20,30,40,50,150};
    for (int index=0 ; index < 2 ; index++){
        noBias[index] = CreateH2F(Form("noBias_%i",index),"",5,binX,6,binY);
        eventHLT[index] = CreateH2F(Form("eventHLT_%i",index),"",5,binX,6,binY);
        eventFindApair[index] = CreateH2F(Form("eventFindApair_%i",index),"",5,binX,6,binY);
        eventFindAtag[index] = CreateH2F(Form("eventFindAtag_%i",index),"",5,binX,6,binY);
        eventTagMatchedWithHlt[index] = CreateH2F(Form("eventTagMatchedWithHlt_%i",index),"",5,binX,6,binY);
        eventInMassWindows[index] = CreateH2F(Form("eventInMassWindows_%i",index),"",5,binX,6,binY);
        eventProbeOnly0[index] = CreateH2F(Form("eventProbeOnly0_%i",index),"",5,binX,6,binY);
        eventProbeOnly1[index] = CreateH2F(Form("eventProbeOnly1_%i",index),"",5,binX,6,binY);
        eventProbeOnly2[index] = CreateH2F(Form("eventProbeOnly2_%i",index),"",5,binX,6,binY);
        eventBeforeSaving[index] = CreateH2F(Form("eventBeforeSaving_%i",index),"",5,binX,6,binY);
     
    }
    
    
    theTree = CreateTree("fitter_tree","the tree for Tag and Probe");
    theTree->Branch("mass",          &mass,          "mass/F");
    theTree->Branch("runNumber",          &runNumber,          "runNumber/I");
    theTree->Branch("lumiSec",          &lumiSec,          "lumiSec/I");
    theTree->Branch("pt",          &pt,          "pt/F");
    theTree->Branch("absEta",          &absEta,          "absEta/F");
    theTree->Branch("SCeta",          &SCeta,          "SCeta/F");
    theTree->Branch("tag_pt",          &tag_pt,          "tag_pt/F");
    theTree->Branch("tag_absSCEta",          &tag_absSCEta,          "tag_absSCEta/F");
    theTree->Branch("tp_deltaR",          &tp_deltaR,          "tp_deltaR/F");
    theTree->Branch("absSCeta",          &absSCeta,          "absSCeta/F");
    theTree->Branch("isFSR",          &isFSR,          "isFSR/I");
    theTree->Branch("eventMatched",   &eventMatched, "eventMatched/I");
    theTree->Branch("eventMatched2",   &eventMatched2, "eventMatched2/I");
    theTree->Branch("eventMatched3",   &eventMatched3, "eventMatched3/I");
    theTree->Branch("tag_isFSR",          &tag_isFSR,          "tag_isFSR/I");
    theTree->Branch("nVtx",          &nVtx,          "nVtx/I");
    theTree->Branch("weight",&weight,"weight/F");
    theTree->Branch("weight_runA",&weight_runA,"weight_runA/F");
    theTree->Branch("weight_runB",&weight_runB,"weight_runB/F");
    theTree->Branch("weight_runC",&weight_runC,"weight_runC/F");
    theTree->Branch("weight_runD",&weight_runD,"weight_runD/F");
    theTree->Branch("weight_runAsingle",&weight_runAsingle,"weight_runAsingle/F");
    theTree->Branch("weight_runDsingle",&weight_runDsingle,"weight_runDsingle/F");
    theTree->Branch("passTight",          &passTight,          "passTight/I");
    theTree->Branch("passLoose",          &passLoose,          "passLoose/I");
    theTree->Branch("passFO",          &passFO,          "passFO/I");
    theTree->Branch("passBDT",          &passBDT,          "passBDT/I");
    theTree->Branch("passISO",          &passISO,          "passISO/I");
    theTree->Branch("passFO_BDT",          &passFO_BDT,          "passFO_BDT/I");
    theTree->Branch("passFO_ISO",          &passFO_ISO,          "passFO_ISO/I");
    theTree->Branch("passFO_BDT_ISO",          &passFO_BDT_ISO,          "passFO_BDT_ISO/I");
    theTree->Branch("passBDT_ISO",          &passBDT_ISO,          "passBDT_ISO/I");
    theTree->Branch("passNM1IP",          &passNM1IP,          "passNM1IP/I");
    theTree->Branch("passNM1presel",          &passNM1presel,          "passNM1presel/I");
    theTree->Branch("passNM1convs",          &passNM1convs,          "passNM1convs/I");
    theTree->Branch("passAllNoIsoDet",          &passAllNoIsoDet,          "passAllNoIsoDet/I");
    theTree->Branch("trigSingle",          &trigSingle,          "trigSingle/I");
    theTree->Branch("trigDoubleLeg0",          &trigDoubleLeg0,          "trigDoubleLeg0/I");
    theTree->Branch("trigDoubleLeg1",          &trigDoubleLeg1,          "trigDoubleLeg1/I");
    theTree->Branch("passPreselec",          &passPreselec,          "passPreselec/I");
    theTree->Branch("passIP",          &passIP,          "passIP/I");
    theTree->Branch("passFOnoIso",          &passFOnoIso,          "passFOnoIso/I");
    theTree->Branch("passConvs",          &passConvs,          "passConvs/I");
    theTree->Branch("isSameSign",          &isSameSign,          "isSameSign/I");
    theTree->Branch("topSelection",          &topSelection,          "topSelection/I");
    theTree->Branch("topFO",          &topFO,          "topFO/I");

    testTree = CreateTree("signalTree","treeWithTheSigEvents");
    testTree->Branch("pt",          &pt,          "pt/F");
  //  testTree->Branch("absEta",          &absEta,          "absEta/F");
    testTree->Branch("SCeta",          &SCeta,          "SCeta/F");
    testTree->Branch("absSCeta",          &absSCeta,          "absSCeta/F");
    testTree->Branch("passTight",          &passTight,          "passTight/I");
  //  testTree->Branch("passLoose",          &passLoose,          "passLoose/I");
    testTree->Branch("passFO",          &passFO,          "passFO/I");
    testTree->Branch("passBDT",          &passBDT,          "passBDT/I");
   // testTree->Branch("passFO_BDT",          &passFO_BDT,          "passFO_BDT/I");
    //testTree->Branch("isTriggering",          &isTriggering,          "isTriggering/I");
    testTree->Branch("theType",          &theType,          "theType/I");
    testTree->Branch("isPair",          &isPair,          "isPair/I");
   // testTree->Branch("MVAidTrig",          &MVAidTrig,          "MVAidTrig/F");
    testTree->Branch("MVAiso",          &MVAiso,          "MVAiso/F");
  //  testTree->Branch("radIso",          &radIso,          "radIso/F");
  //  testTree->Branch("radIsoVeto",          &radIsoVeto,          "radIsoVeto/F");
  //  testTree->Branch("radIsoVetoMore",          &radIsoVetoMore,          "radIsoVetoMore/F");
    testTree->Branch("combIsoHWW",          &combIsoHWW,          "combIsoHWW/F");
  /// testTree->Branch("passMVA",          &passMVA,          "passMVA/I");
    testTree->Branch("PDGid",          &PDGid,          "PDGid/I");
    testTree->Branch("passTightIdNoIsoPart", &passTightIdNoIsoPart, "passTightIdNoIsoPart/I");
  //  testTree->Branch("passTightLooserIdNoIsoPart", &passTightLooserIdNoIsoPart, "passTightLooserIdNoIsoPart/I");
    testTree->Branch("POGisolation", &POGisolation, "POGisolation/F");
//    testTree->Branch("POGisolation2012A", &POGisolation2012A, "POGisolation2012A/F");
 //   testTree->Branch("myPOGisolation", &myPOGisolation, "myPOGisolation/F");
  ///  testTree->Branch("radIsoStandard", &radIsoStandard, "radIsoStandard/F");
  //  testTree->Branch("radIsoNoInner", &radIsoNoInner, "radIsoNoInner/F");
   // testTree->Branch("radisoNoThreshold", &radisoNoThreshold, "radisoNoThreshold/F");
   // testTree->Branch("radisoNoThresholdNoInner", &radisoNoThresholdNoInner, "radisoNoThresholdNoInner/F");
  //  testTree->Branch("radIsoBeta", &radIsoBeta, "radIsoBeta/F");
  //  testTree->Branch("POGnoRho", &POGnoRho, "POGnoRho/F");
  //  testTree->Branch("PFisoDeltaBeta", &PFisoDeltaBeta, "PFisoDeltaBeta/F");
  //  testTree->Branch("myTightID", &myTightID, "myTightID/I");
  //  testTree->Branch("myTightIDiso", &myTightIDiso, "myTightIDiso/I");
    testTree->Branch("nbRecoVertex", &nbRecoVertex, "nbRecoVertex/I");
    
  //  testTree->Branch("passFO_BDT_ISO",          &passFO_BDT_ISO,          "passFO_BDT_ISO/I");

    
    testTree->Branch("radisoFonc0", &radisoFonc0, "radisoFonc0/F");
    testTree->Branch("radisoFonc1", &radisoFonc1, "radisoFonc1/F");
    testTree->Branch("radisoFonc2", &radisoFonc2, "radisoFonc2/F");
    testTree->Branch("radisoFonc3", &radisoFonc3, "radisoFonc3/F");

    testTree->Branch("chRadisoFonc0", &chRadisoFonc0, "chRadisoFonc0/F");
    testTree->Branch("chRadisoFonc1", &chRadisoFonc1, "chRadisoFonc1/F");
    testTree->Branch("chRadisoFonc2", &chRadisoFonc2, "chRadisoFonc2/F");
    testTree->Branch("chRadisoFonc3", &chRadisoFonc3, "chRadisoFonc3/F");
    
    testTree->Branch("nhRadisoFonc0", &nhRadisoFonc0, "nhRadisoFonc0/F");
    testTree->Branch("nhRadisoFonc1", &nhRadisoFonc1, "nhRadisoFonc1/F");
    testTree->Branch("nhRadisoFonc2", &nhRadisoFonc2, "nhRadisoFonc2/F");
    testTree->Branch("nhRadisoFonc3", &nhRadisoFonc3, "nhRadisoFonc3/F");
   
    testTree->Branch("gRadisoFonc0", &gRadisoFonc0, "gRadisoFonc0/F");
    testTree->Branch("gRadisoFonc1", &gRadisoFonc1, "gRadisoFonc1/F");
    testTree->Branch("gRadisoFonc2", &gRadisoFonc2, "gRadisoFonc2/F");
    testTree->Branch("gRadisoFonc3", &gRadisoFonc3, "gRadisoFonc3/F");
   
    testTree->Branch("chIso0", &chIso0, "chIso0/F");
    testTree->Branch("chIso1", &chIso1, "chIso1/F");
    testTree->Branch("chIso2", &chIso2, "chIso2/F");
    testTree->Branch("chIso3", &chIso3, "chIso3/F");
    testTree->Branch("chIso4", &chIso4, "chIso4/F");
	
    testTree->Branch("nhIso0", &nhIso0, "nhIso0/F");
    testTree->Branch("nhIso1", &nhIso1, "nhIso1/F");
    testTree->Branch("nhIso2", &nhIso2, "nhIso2/F");
    testTree->Branch("nhIso3", &nhIso3, "nhIso3/F");
    testTree->Branch("nhIso4", &nhIso4, "nhIso4/F");

    testTree->Branch("gIso0", &gIso0, "gIso0/F");
    testTree->Branch("gIso1", &gIso1, "gIso1/F");
    testTree->Branch("gIso2", &gIso2, "gIso2/F");
    testTree->Branch("gIso3", &gIso3, "gIso3/F");
    testTree->Branch("gIso4", &gIso4, "gIso4/F");

    testTree->Branch("gComb0", &gComb0, "gComb0/F");
    testTree->Branch("gComb1", &gComb1, "gComb1/F");
    testTree->Branch("gComb2", &gComb2, "gComb2/F");
    testTree->Branch("gComb3", &gComb3, "gComb3/F");
    testTree->Branch("gComb4", &gComb4, "gComb4/F");

    
    isoInfo  = CreateTree("isoInfo","iso info");
    isoInfo->Branch("PFdeltaEta", &PFdeltaEta, "PFdeltaEta/F");
    isoInfo->Branch("PFdeltaPhi", &PFdeltaPhi, "PFdeltaPhi/F");
    isoInfo->Branch("PFet", &PFet, "PFet/F");
    isoInfo->Branch("PFtype", &PFtype, "PFtype/I");
    isoInfo->Branch("PFisPU", &PFisPU, "PFisPU/I");
    isoInfo->Branch("mainEta", &mainEta, "mainEta/F");
    isoInfo->Branch("mainSCEta", &mainSCEta, "mainSCEta/F");
    isoInfo->Branch("mainPhi", &mainPhi, "mainPhi/F");
    isoInfo->Branch("mainEnergy", &mainEnergy, "mainEnergy/F");
    isoInfo->Branch("mainPt", &mainPt, "mainPt/F");
    isoInfo->Branch("isSignalBg", &isSignalBg, "isSignalBg/I");
    isoInfo->Branch("PFdeltaR", &PFdeltaR, "PFdeltaR/F");




}

void HWWAnalysis::InsideLoop() {
    
  // if ((T_Event_RunNumber>207099)&&(T_Event_RunNumber<208686)) return;
   // if (!((T_Event_RunNumber==191226)&&(T_Event_LuminosityBlock>878)&&(T_Event_LuminosityBlock<1003))) return;
   // if (!((T_Event_RunNumber==206745)&&(T_Event_LuminosityBlock>765)&&(T_Event_LuminosityBlock<874))) return;
//if (!((T_Event_RunNumber==203912)&&(T_Event_LuminosityBlock>710)&&(T_Event_LuminosityBlock<818))) return;
  // print the event number
   // cout << "event=" << T_Event_EventNumber << endl;
  // The InsideLoop() function is called for each entry in the tree to be 
  // processed  
/*    if (signal == "electrons_runA"){
        if ((T_Event_RunNumber==190949)||(T_Event_RunNumber==191090)||(T_Event_RunNumber==191367)||(T_Event_RunNumber==193112)||(T_Event_RunNumber==193116)) return;
    }*/
   // if (!(isAGoodEvent(T_Event_RunNumber, T_Event_LuminosityBlock))) return;

    runNumber = T_Event_RunNumber;
    lumiSec = T_Event_LuminosityBlock;
    TLorentzVector* elecTag;
    TLorentzVector*  elecProb;
    TLorentzVector sumElec;
	passCuts->Fill("no cut",1);
//	cout << "nb of elec in the event = " << T_Elec_Pt->size() << endl;
    int nbElec = T_Elec_Pt->size();

    if (nbElec == 0) return;
    else if((nbElec == 1)&&(T_METPF_ET<20)&&(T_Elec_Pt->at(0)>10)){
      // fillTheTestTree(0,0);
        return;
    }
    else if ((nbElec == 2)&&(T_METPF_ET<20)){
        elecTag = new TLorentzVector(T_Elec_Px->at(0),T_Elec_Py->at(0),T_Elec_Pz->at(0),T_Elec_Energy->at(0));
        elecProb = new TLorentzVector(T_Elec_Px->at(1),T_Elec_Py->at(1),T_Elec_Pz->at(1),T_Elec_Energy->at(1));
        sumElec = *elecTag + *elecProb;
        mass = sumElec.M();
        if ((mass<60)&&(mass>120)){
            if ( T_Elec_Pt->at(0) >10){
               // fillTheTestTree(0,0);
            }
            if ( T_Elec_Pt->at(1) >10){
             //  fillTheTestTree(0,0);
            }
        }
    }
    for (int j = 0 ; j < nbElec; j++){
        if ((isMC) && (fabs(T_Gen_Elec_PDGid->at(j))==11)) {noBias[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) noBias[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); }
    }
    if (T_Event_HLT_Ele27_WP80 == 0) return;
   // if ((T_Event_HLT_Ele17_Ele8_M50_TnP == 0)&&(T_Event_HLT_Ele20_SC4_M50_TnP==0)) return;
    //passCuts->Fill("passTrig",1);
   /* cout << "Event Rho " << T_Event_Rho << endl;
    cout << "Event Run Number = " << T_Event_RunNumber << endl;
    cout << "nb of electrons " << nbElec << endl;
    cout << "pass trigger " << T_Event_HLT_Ele27_WP80 << endl;
    cout << "test size = " << T_Elec_HLT_Elec27_WP80->size() << endl;*/
    if (isMC) compWeights(T_Event_nTruePU);
	//TLorentzVector *elec1 = new TLorentzVector(1,1,1,1);
    for (int j = 0 ; j < nbElec; j++){
        if ((isMC) && (fabs(T_Gen_Elec_PDGid->at(j))==11)) {eventHLT[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventHLT[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1);}
        if ((isMC) && (fabs(T_Gen_Elec_PDGid->at(j))==11)) {if (nbElec>1) { eventFindApair[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventFindApair[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); }}
        bool findAGoodTag = false;
        bool findATagMatched = false;
        bool findEventInMassWindows = false;
        for (int i = 0 ; i < nbElec ; i++){
            if (i==j) continue;
	/*  if (T_Event_EventNumber==996760) {
           cout << "nb " << nbElec << endl;

            elecTag = new TLorentzVector(T_Elec_Px->at(i),T_Elec_Py->at(i),T_Elec_Pz->at(i),T_Elec_Energy->at(i));
            elecProb = new TLorentzVector(T_Elec_Px->at(j),T_Elec_Py->at(j),T_Elec_Pz->at(j),T_Elec_Energy->at(j));
            sumElec = *elecTag + *elecProb;
            mass = sumElec.M();	
          cout << "trig = " << T_Event_HLT_Ele17_Ele8_M50_TnP << " " << T_Event_HLT_Ele20_SC4_M50_TnP << endl;
          cout << "pres =" << (T_Elec_HLT_Ele17_Ele8_TnP_Ele17Leg->at(i)) << " " << (T_Elec_HLT_Ele20_SC4_TnP_Ele20Leg->at(i)) << endl;
          cout << "pres2 =" << (T_Elec_HLT_Ele17_Ele8_TnP_Ele8Leg->at(j)) << " " << (T_Elec_HLT_Ele20_SC4_TnP_SC4Leg->at(j)) << endl;
            cout <<  T_Event_RunNumber << " " << T_Event_EventNumber << " " << mass << " " << T_Elec_Pt->at(i) << " " << T_Elec_Charge->at(i) << " " << T_Elec_Pt->at(j) << " " << T_Elec_Charge->at(j) << endl;
           cout << "elec1 = eta " << T_Elec_SC_Eta->at(i) << "ID " << T_Elec_passTight->at(i) << endl;
           cout << "elec2 = eta " << T_Elec_SC_Eta->at(j) << "ID " << T_Elec_passTight->at(j) << endl;
	   }
	    continue;*/
            passCuts->Fill("nb pairs",1);
           // if (!((T_Elec_Pt->at(i)>20)&&(T_Elec_passTight->at(i)==1)&&(fabs(T_Elec_SC_Eta->at(i))<2.5)&&(T_Elec_HLT_Elec27_WP80->at(i)==1))) continue;
            if (!(T_Elec_Pt->at(i)>25)) continue;
            passCuts->Fill("tag>20GeV",1);
            if (!(fabs(T_Elec_SC_Eta->at(i))<2.5)) continue;
            passCuts->Fill("tagEta<2.5",1);
          //  float theTagIso = giveThePOGiso(T_Elec_Pt->at(i), T_Elec_neutralHadronIso->at(i), T_Elec_chargedHadronIso->at(i), T_Elec_photonIso->at(i), SCeta, T_Event_RhoIso);
           // float ooemoopTag = 1.0/T_Elec_EcalEnergy->at(i) - T_Elec_eSuperClusterOverP->at(i)/T_Elec_EcalEnergy->at(i);
           // bool tagPassLoose = passTightIdNoIsoModified(T_Elec_isEB->at(i), T_Elec_deltaEtaIn->at(i), T_Elec_deltaPhiIn->at(i), T_Elec_sigmaIetaIeta->at(i), T_Elec_HtoE->at(i), ooemoopTag , T_Elec_d0->at(i), T_Elec_dZ->at(i), T_Elec_nHits->at(i), T_Elec_passConversionVeto->at(i));
           // float cutIsoValue = 0.1;
          //  if ((T_Elec_Pt->at(i)<20)&&(T_Elec_isEB==0)){
         //       cutIsoValue = 0.07;
         //   }
          //  bool caPasse = tagPassLoose&&(theTagIso<cutIsoValue);
            //cout << "on passe tight ? " << T_Elec_passTight->at(i) << endl;
            //cout << "on passe mytight ? " << caPasse << endl;
            //if (!(tagPassLoose&&(theTagIso<cutIsoValue))) continue;
            if (!(T_Elec_passTight->at(i)==1)) continue;
            findAGoodTag = true;
            //if (!(T_Elec_passMedium->at(i)==1)) continue;
            passCuts->Fill("tagPassTight",1);
            if (!(T_Elec_HLT_Elec27_WP80->at(i)==1)) continue;
            findATagMatched = true;
            passCuts->Fill("tagMatchWithTrigger",1);
            elecTag = new TLorentzVector(T_Elec_Px->at(i),T_Elec_Py->at(i),T_Elec_Pz->at(i),T_Elec_Energy->at(i));
            elecProb = new TLorentzVector(T_Elec_Px->at(j),T_Elec_Py->at(j),T_Elec_Pz->at(j),T_Elec_Energy->at(j));
            sumElec = *elecTag + *elecProb;
            mass = sumElec.M();
           if ((mass > 60)&&(mass<120)){
                findEventInMassWindows = true;
                passCuts->Fill("inMassWindow",1);
                isFSR = 0;
                tag_isFSR = 0;
               eventMatched = 0;
               eventMatched2 = 0;
               eventMatched3 = 0;
                if (isMC){
                    int statusOftag = isGenMatched(T_Gen_Elec_PDGid->at(i), T_Gen_Elec_MotherID->at(i), T_Gen_Elec_status->at(i));
                    int statusOfprobe = isGenMatched(T_Gen_Elec_PDGid->at(j), T_Gen_Elec_MotherID->at(j), T_Gen_Elec_status->at(j));
                    if (((statusOftag>0)&&(statusOfprobe>0))) eventMatched=1;
                    if (((statusOftag>1)&&(statusOfprobe>1))) eventMatched2=1;
                    if (((statusOftag>2)&&(statusOfprobe>2))) eventMatched3=1;
                    if (statusOftag==2) tag_isFSR = 1;
                    if (statusOfprobe==2) isFSR = 1;
               //     cout << "le tag : PDGID=" << T_Gen_Elec_PDGid->at(i) << ", status=" << T_Gen_Elec_status->at(i) << ", mother=" << T_Gen_Elec_MotherID->at(i) << " deltaR=" << T_Gen_Elec_deltaR->at(i) << endl;
                 //   cout << "le probe : PDGID=" << T_Gen_Elec_PDGid->at(j) << ", status=" << T_Gen_Elec_status->at(j) << ", mother=" << T_Gen_Elec_MotherID->at(j) << " deltaR=" << T_Gen_Elec_deltaR->at(j) << endl;
                
             //   cout << "we will store this event ! " << endl;
                    //if (statusOfprobe>0) cout << "le j de l'interieur " << endl;
                    if (statusOfprobe>0) { eventProbeOnly0[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventProbeOnly0[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1);}
                    if (statusOfprobe>1) { eventProbeOnly1[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventProbeOnly1[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1);}
                    if (statusOfprobe>2) { eventProbeOnly2[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventProbeOnly2[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1);}
               }
               if (eventMatched==1) { eventBeforeSaving[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventBeforeSaving[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1);}
                pt = T_Elec_Pt->at(j);
                eta = T_Elec_Eta->at(j);
                SCeta = T_Elec_SC_Eta->at(j);
                absEta = fabs(eta);
                absSCeta = fabs(SCeta);
                tag_absSCEta = fabs(T_Elec_SC_Eta->at(i));
                tag_pt = T_Elec_Pt->at(i);
                TLorentzVector *theTag = new TLorentzVector(T_Elec_Px->at(j),T_Elec_Py->at(j),T_Elec_Pz->at(j),T_Elec_Energy->at(j));
                TLorentzVector *theProbe = new TLorentzVector(T_Elec_Px->at(i),T_Elec_Py->at(i),T_Elec_Pz->at(i),T_Elec_Energy->at(i));
                tp_deltaR = sqrt(pow(eta-T_Elec_Eta->at(j),2)+ pow(acos(cos(theProbe->Phi()-theTag->Phi())),2)) ;
                nVtx = T_Vertex_z->size();
                passLoose = T_Elec_passLoose->at(j);
                passTight = T_Elec_passTight->at(j);
                passFO = T_Elec_isFO->at(j);
                passBDT = T_Elec_passMVA->at(j);
                passISO = ((T_Elec_CombIsoHWW->at(j)/pt)<0.15);
                passFO_BDT = passFO && passBDT;
                passFO_ISO = passFO && passISO;
                passFO_BDT_ISO = passFO_BDT && passISO;
               passBDT_ISO = passBDT && passISO;
                passPreselec = passPreCuts(T_Elec_Pt->at(j), T_Elec_isEB->at(j), T_Elec_sigmaIetaIeta->at(j), T_Elec_deltaEtaIn->at(j), T_Elec_deltaPhiIn->at(j) ,T_Elec_HtoE->at(j));
                passIP = passIPcuts( T_Elec_d0->at(j), T_Elec_dZ->at(j));
               passConvs = passMissItCons(T_Elec_passConversionVeto->at(j), T_Elec_nHits->at(j));
               passFOnoIso = FOnoIso(T_Elec_Pt->at(j), T_Elec_isEB->at(j), T_Elec_sigmaIetaIeta->at(j), T_Elec_deltaEtaIn->at(j), T_Elec_deltaPhiIn->at(j) ,T_Elec_HtoE->at(j),T_Elec_d0->at(j), T_Elec_dZ->at(j), T_Elec_passConversionVeto->at(j), T_Elec_nHits->at(j));
               passAllNoIsoDet = passFOnoIso && passISO && passBDT;
               passNM1IP = passPreselec && passConvs && passBDT_ISO;
               passNM1convs = passPreselec && passIP && passBDT_ISO;
               passNM1presel = passConvs && passIP && passBDT_ISO;
               trigSingle = T_Elec_HLT_Elec27_WP80->at(j);
               trigDoubleLeg0 = T_Elec_HLT_Ele17_Ele8_Ele8Leg->at(j);
               trigDoubleLeg1 = T_Elec_HLT_Ele17_Ele8_Ele17Leg->at(j);
               float the03RelIsolation = calc03Iso(pt, SCeta , T_Event_RhoIso, j);
               topSelection = (T_Elec_MVAid_trig->at(j)>0.5)&&passConvs&&(the03RelIsolation<0.15)&&(T_Elec_d0->at(j)<0.04);
               //cout << "pass? " << (T_Elec_MVAid_trig->at(j)>0.5) << " " << passConvs << " " << (the03RelIsolation<0.15) << " " << (T_Elec_d0->at(j)<0.04) << endl;
               topFO = (T_Elec_MVAid_trig->at(j)>-0.1)&&passConvs&&(the03RelIsolation<1.0)&&(T_Elec_d0->at(j)<0.1);

                isSameSign = (T_Elec_Charge->at(i)*T_Elec_Charge->at(j)==1 ? 1 : 0);
               
               float pass=0;
               //if (passTight==1 && T_Elec_passTight->at(j)==1) pass=1;
	           //if ((T_Elec_Charge->at(i)*T_Elec_Charge->at(j))==1) continue;
               //if (pt < 10 ) continue;
               //if ((T_Elec_HLT_Ele17_Ele8_TnP_Ele17Leg->at(i)==0)&&(T_Elec_HLT_Ele20_SC4_TnP_Ele20Leg->at(i)==0)) continue;
                // cout <<  T_Event_RunNumber << " " << T_Event_EventNumber << " " << mass << " " << T_Elec_Pt->at(i) << " " << T_Elec_Charge->at(i) << " " << pt << " " << T_Elec_Charge->at(j) << " " << pass << endl;
               theTree->Fill();
               hMass->Fill(mass);
               fillTheHistoMass(j, mass);
               if ((mass<100)&&(mass>80)&&(T_Elec_Pt->at(j)>10)) {
                   fillTheHisto(j);
                 //  fillTheTestTree(j,1);
               }

            }
        }
        if ((isMC) && (fabs(T_Gen_Elec_PDGid->at(j))==11)) { if (findAGoodTag) { eventFindAtag[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventFindAtag[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); }}
        if ((isMC) && (fabs(T_Gen_Elec_PDGid->at(j))==11)) { if (findATagMatched) { eventTagMatchedWithHlt[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventTagMatchedWithHlt[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); }}
        if ((isMC) && (fabs(T_Gen_Elec_PDGid->at(j))==11)) { if (findEventInMassWindows) { eventInMassWindows[0]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); if(T_Elec_passTight->at(j)==1) eventInMassWindows[1]->Fill(fabs(T_Elec_SC_Eta->at(j)),T_Elec_Pt->at(j),1); }}
      //  if (findEventInMassWindows) cout << "on a un tag j=" << j << endl;
        
    }
    
    
 /*   pt = T_Elec_Pt->at(0);
    eta = T_Elec_Eta->at(0);
    theTree->Fill();*/


}  

int HWWAnalysis::isGenMatched(int part_PDGID, int mother_PDGID, int status){
    //cout << "coucou mother id " << mother_PDGID << " status = " << status << endl;
    if ((fabs(part_PDGID)==11)&&(fabs(mother_PDGID==23))&&(status>0)) return 3; //before FSR case
    if ((fabs(part_PDGID)==11)&&(fabs(mother_PDGID==22))&&(status>0)) return 2; //before FSR case
    if ((fabs(part_PDGID)==11)&&(fabs(mother_PDGID==11))&&(status<3)) return 1; //after FSR case
    return 0;
}

void HWWAnalysis::compWeights(int nbVtx){
    if (nbVtx>=50) nbVtx = 49;
    if (nbVtx==0) nbVtx = 1;
    float runAandBweights[50] = {0.366495, 0.728276, 3.19901, 204.661, 18.9332, 5.71718, 22.7843, 44.2697, 69.0199, 92.2806, 92.6908, 102.041, 82.3959, 72.988, 65.5631, 52.8622, 36.2083, 23.0384, 14.9948, 10.596, 8.10421, 6.56447, 5.48442, 4.65863, 3.98905, 3.37014, 2.83076, 2.33678, 1.9154, 1.54498, 1.2326, 0.985676, 0.788231, 0.625784, 0.494525, 0.39218, 0.312244, 0.245845, 0.197023, 0.157868, 0.126589, 0.104034, 0.0864918, 0.0725297, 0.0643046, 0.0561723, 0.0520573, 0.050254, 0.0476669, 1.59261};
    float runAweights[50] = {0.318681, 0.339003, 0.354653, 0.403267, 0.510374, 0.656627, 0.825227, 1.0201, 1.21247, 1.3808, 1.52815, 1.62851, 1.69215, 1.71533, 1.71016, 1.67747, 1.61846, 1.54947, 1.46749, 1.37338, 1.27144, 1.16833, 1.05714, 0.947377, 0.838663, 0.733132, 0.634314, 0.54097, 0.455521, 0.379485, 0.311971, 0.254901, 0.205755, 0.163343, 0.129027, 0.100699, 0.0779285, 0.0599513, 0.04564, 0.0343054, 0.0256398, 0.0191075, 0.0140853, 0.0102895, 0.00746788, 0.00539845, 0.00388476, 0.00275792, 0.00196963, 0.00139967};
    float runBweights[50] = {0.755029, 1.04812, 1.19818, 1.26512, 1.35359, 1.42638, 1.47226, 1.52087, 1.54476, 1.53856, 1.52282, 1.48216, 1.43442, 1.37939, 1.32735, 1.27745, 1.22843, 1.18992, 1.15664, 1.12601, 1.0979, 1.07451, 1.04573, 1.01647, 0.982746, 0.943448, 0.900252, 0.849388, 0.792951, 0.733358, 0.669741, 0.60798, 0.545072, 0.480276, 0.420655, 0.363559, 0.311096, 0.264178, 0.221564, 0.183077, 0.150064, 0.122329, 0.0983644, 0.0781441, 0.0614781, 0.0480073, 0.0371809, 0.0282992, 0.0215801, 0.0163058};
    float runCweights[50] = {0.486018, 0.56173, 0.6633, 0.753877, 0.874244, 0.995233, 1.10175, 1.21003, 1.29483, 1.34688, 1.38135, 1.38372, 1.37062, 1.34307, 1.3124, 1.27903, 1.24252, 1.21311, 1.18587, 1.15839, 1.13088, 1.1061, 1.07431, 1.04137, 1.00408, 0.962099, 0.917782, 0.867696, 0.81409, 0.759285, 0.701982, 0.647782, 0.592907, 0.535719, 0.483308, 0.432177, 0.384313, 0.34062, 0.299427, 0.260391, 0.225524, 0.195003, 0.166934, 0.141691, 0.119506, 0.100373, 0.0838754, 0.069088, 0.0571805, 0.0470227};
    float runDweights[50] = {0.189659, 0.271387, 0.34027, 0.404113, 0.48973, 0.583095, 0.674846, 0.773579, 0.861875, 0.93084, 0.988508, 1.02284, 1.04455, 1.05387, 1.05958, 1.06244, 1.06242, 1.06872, 1.07772, 1.08753, 1.09836, 1.11294, 1.12126, 1.12865, 1.13108, 1.12727, 1.11908, 1.10141, 1.07595, 1.04489, 1.00573, 0.965985, 0.91994, 0.864468, 0.810683, 0.753109, 0.695325, 0.639457, 0.582911, 0.525348, 0.471278, 0.421855, 0.373684, 0.328072, 0.28612, 0.248437, 0.214596, 0.182713, 0.156328, 0.132928};
    float runAsingleweights[50] = {0.0556421, 0.131913, 0.245752, 0.403764, 0.641765, 0.958889, 1.33926, 1.78689, 2.23844, 2.62731, 2.9303, 3.07544, 3.07365, 2.92522, 2.67165, 2.34206, 1.97036, 1.60524, 1.26312, 0.959508, 0.704922, 0.503016, 0.346199, 0.231395, 0.149965, 0.0943136, 0.057754, 0.034332, 0.0198652, 0.011222, 0.00617867, 0.00334219, 0.00176676, 0.000909211, 0.000461103, 0.000228939, 0.00011173, 5.37535e-05, 2.53839e-05, 1.1742e-05, 5.35913e-06, 2.4203e-06, 1.07307e-06, 4.67931e-07, 2.01209e-07, 8.55278e-08, 3.59175e-08, 1.47683e-08, 6.06217e-09, 2.4572e-09};
    float runDsingleweights[50] = {0.0734205, 0.166531, 0.298051, 0.472309, 0.726821, 1.05524, 1.43715, 1.87606, 2.30684, 2.66601, 2.93663, 3.05277, 3.0305, 2.87258, 2.61995, 2.29942, 1.94154, 1.59134, 1.2627, 0.969442, 0.721415, 0.522547, 0.365825, 0.24922, 0.164952, 0.106148, 0.0666358, 0.0406822, 0.0242188, 0.0141006, 0.00801517, 0.00448351, 0.00245493, 0.00131066, 0.000690659, 0.00035685, 0.000181502, 9.11369e-05, 4.49824e-05, 2.17787e-05, 1.0418e-05, 4.9379e-06, 2.30068e-06, 1.05567e-06, 4.78258e-07, 2.14453e-07, 9.51201e-08, 4.1358e-08, 1.79734e-08, 7.72183e-09};
    float allHCPdataset[50] = {0.430197, 0.55548, 0.652449, 0.725975, 0.825586, 0.926318, 1.01602, 1.11065, 1.18723, 1.23701, 1.27325, 1.28168, 1.27676, 1.25874, 1.23778, 1.21416, 1.18747, 1.16768, 1.15034, 1.13333, 1.11695, 1.10396, 1.08458, 1.06439, 1.03981, 1.01008, 0.977264, 0.937321, 0.89226, 0.844345, 0.791948, 0.74129, 0.688104, 0.63042, 0.576587, 0.522621, 0.471031, 0.42311, 0.376965, 0.332277, 0.291743, 0.25579, 0.222104, 0.19129, 0.163786, 0.139728, 0.118671, 0.0994171, 0.0837516, 0.070164};
    weight =  allHCPdataset[nbVtx-1];
    weight_runA = runAweights[nbVtx-1];
    weight_runB = runBweights[nbVtx-1];
    weight_runC = runCweights[nbVtx-1];
    weight_runD = runDweights[nbVtx-1];
    weight_runAsingle = runAsingleweights[nbVtx-1];
    weight_runDsingle = runDsingleweights[nbVtx-1];
    
}
void HWWAnalysis::fillTheHisto(int theElec){
    for (int index=0 ; index < 9 ; index++){
        switch (index) {
            case 1:
                if (!(fabs(T_Elec_SC_Eta->at(theElec))<1.479)) continue; 
                break;
            case 2:
                if (!(fabs(T_Elec_SC_Eta->at(theElec))>1.479)) continue; 
                break;
            case 3:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))<0.8)&&(T_Elec_Pt->at(theElec))<20)) continue; 
                break;    
            case 4:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))<0.8)&&(T_Elec_Pt->at(theElec))>20)) continue; 
                break;    
            case 5:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>0.8)&&(fabs(T_Elec_SC_Eta->at(theElec))<1.479)&&(T_Elec_Pt->at(theElec)<20))) continue; 
                break;    
            case 6:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>0.8)&&(fabs(T_Elec_SC_Eta->at(theElec))<1.479)&&(T_Elec_Pt->at(theElec)>20))) continue;
                break;    
            case 7:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>1.479)&&(fabs(T_Elec_SC_Eta->at(theElec))<2.5)&&(T_Elec_Pt->at(theElec)<20))) continue; 
                break;    
            case 8:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>1.479)&&(fabs(T_Elec_SC_Eta->at(theElec))<2.5)&&(T_Elec_Pt->at(theElec)>20))) continue; 
                break;    
            default:
                break;
        }
        hkfchi2[index]->Fill(T_Elec_kfchi2->at(theElec));
        hkfhits[index]->Fill(T_Elec_kfhits->at(theElec));
        hgsfchi2[index]->Fill(T_Elec_gsfchi2->at(theElec));
        hdetacalo[index]->Fill(T_Elec_detacalo->at(theElec));
        hsee[index]->Fill(T_Elec_see->at(theElec));
        hspp[index]->Fill(T_Elec_spp->at(theElec));
        hetawidth[index]->Fill(T_Elec_etawidth->at(theElec));
        hphiwidth[index]->Fill(T_Elec_phiwidth->at(theElec));
        he1x5e5x5[index]->Fill(T_Elec_e1x5e5x5->at(theElec));
        hr9[index]->Fill(T_Elec_R9->at(theElec));
        hEoP[index]->Fill(T_Elec_EoP->at(theElec));
        hIoEmIeP[index]->Fill(T_Elec_IoEmIoP->at(theElec));
        hEoPout[index]->Fill(T_Elec_eleEoPout->at(theElec));
        hPreShowerOverRaw[index]->Fill(T_Elec_PreShowerOverRaw->at(theElec));
        hd0[index]->Fill(T_Elec_d0->at(theElec));
        hIP3D[index]->Fill(T_Elec_IP3D->at(theElec));
        hdZ[index]->Fill(T_Elec_dZ->at(theElec));
        hConbIsoHWW[index]->Fill(T_Elec_CombIsoHWW->at(theElec));
        hMVAid[index]->Fill(T_Elec_MVAid_trig->at(theElec));
        hMvaiso[index]->Fill(T_Elec_Mvaiso->at(theElec));
        hRadiaIso[index]->Fill(T_Elec_RadialIso->at(theElec));
        hRadiaIsoVeto[index]->Fill(T_Elec_RadialIsoVeto->at(theElec));
        hRadiaIsoVetoMore[index]->Fill(T_Elec_RadialIsoVetoMore->at(theElec));
        hchargedHadronIso[index]->Fill(T_Elec_chargedHadronIso->at(theElec));
        hneutralHadronIso[index]->Fill(T_Elec_neutralHadronIso->at(theElec));
        hphotonIso[index]->Fill(T_Elec_photonIso->at(theElec));
        hchargedHadronIso04[index]->Fill(T_Elec_chargedHadronIso04->at(theElec));
        hneutralHadronIso04[index]->Fill(T_Elec_neutralHadronIso04->at(theElec));
        hphotonIso04[index]->Fill(T_Elec_photonIso04->at(theElec));
        hpassConversionVeto[index]->Fill(T_Elec_passConversionVeto->at(theElec));
        hsigmaIetaIeta[index]->Fill(T_Elec_sigmaIetaIeta->at(theElec));
        hdeltaPhiIn[index]->Fill(T_Elec_deltaPhiIn->at(theElec));
        hdeltaEtaIn[index]->Fill(T_Elec_deltaEtaIn->at(theElec));
        hisEcalDriven[index]->Fill(T_Elec_isEcalDriven->at(theElec));
        hnLost[index]->Fill(T_Elec_nLost->at(theElec));
        hnHits[index]->Fill(T_Elec_nHits->at(theElec));
        hfBrem[index]->Fill(T_Elec_fBrem->at(theElec));
        hPt[index]->Fill(T_Elec_Pt->at(theElec));
        hEta[index]->Fill(T_Elec_Eta->at(theElec));
    }
}
void HWWAnalysis::fillTheHistoMass(int theElec, float theMass){
    for (int index=0 ; index < 9 ; index++){
        switch (index) {
            case 1:
                if (!(fabs(T_Elec_SC_Eta->at(theElec))<1.479)) continue;
                break;
            case 2:
                if (!(fabs(T_Elec_SC_Eta->at(theElec))>1.479)) continue;
                break;
            case 3:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))<0.8)&&(T_Elec_Pt->at(theElec))<20)) continue;
                break;
            case 4:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))<0.8)&&(T_Elec_Pt->at(theElec))>20)) continue;
                break;
            case 5:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>0.8)&&(fabs(T_Elec_SC_Eta->at(theElec))<1.479)&&(T_Elec_Pt->at(theElec)<20))) continue;
                break;
            case 6:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>0.8)&&(fabs(T_Elec_SC_Eta->at(theElec))<1.479)&&(T_Elec_Pt->at(theElec)>20))) continue;
                break;
            case 7:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>1.479)&&(fabs(T_Elec_SC_Eta->at(theElec))<2.5)&&(T_Elec_Pt->at(theElec)<20))) continue;
                break;
            case 8:
                if (!((fabs(T_Elec_SC_Eta->at(theElec))>1.479)&&(fabs(T_Elec_SC_Eta->at(theElec))<2.5)&&(T_Elec_Pt->at(theElec)>20))) continue;
                break;
            default:
                break;
        }
        hMassUncut[index]->Fill(theMass);
    }
}

void HWWAnalysis::fillTheTestTree(int ite, int localType){
    pt = T_Elec_Pt->at(ite);
    eta = T_Elec_Eta->at(ite);
    SCeta = T_Elec_SC_Eta->at(ite);
    absEta = fabs(eta);
    absSCeta = fabs(SCeta);
    passLoose = T_Elec_passLoose->at(ite);
    passTight = T_Elec_passTight->at(ite);
    passFO = T_Elec_isFO->at(ite);
    passBDT = T_Elec_passMVA->at(ite);
    passISO = ((T_Elec_CombIsoHWW->at(ite)/pt)<0.15);
    passFO_BDT = passFO && passBDT;
    passFO_BDT_ISO = passFO_BDT && passISO;

    isTriggering = T_Elec_isTrig->at(ite);
    theType=localType;
    MVAidTrig= T_Elec_MVAid_trig->at(ite);
    MVAiso = T_Elec_Mvaiso->at(ite);
    radIso = T_Elec_RadialIso->at(ite);
    radIsoVeto = T_Elec_RadialIsoVeto->at(ite);
    radIsoVetoMore = T_Elec_RadialIsoVetoMore->at(ite);
    combIsoHWW = T_Elec_CombIsoHWW->at(ite);
    passMVA = T_Elec_passMVA->at(ite);
    float ooemoop = 1.0/T_Elec_EcalEnergy->at(ite) - T_Elec_eSuperClusterOverP->at(ite)/T_Elec_EcalEnergy->at(ite);
    passTightIdNoIsoPart = passTightIdNoIso(T_Elec_isEB->at(ite), T_Elec_deltaEtaIn->at(ite), T_Elec_deltaPhiIn->at(ite), T_Elec_sigmaIetaIeta->at(ite), T_Elec_HtoE->at(ite), ooemoop , T_Elec_d0->at(ite), T_Elec_dZ->at(ite), T_Elec_nHits->at(ite), T_Elec_passConversionVeto->at(ite));
    passTightLooserIdNoIsoPart = passTightIdNoIsoModified(T_Elec_isEB->at(ite), T_Elec_deltaEtaIn->at(ite), T_Elec_deltaPhiIn->at(ite), T_Elec_sigmaIetaIeta->at(ite), T_Elec_HtoE->at(ite), ooemoop , T_Elec_d0->at(ite), T_Elec_dZ->at(ite), T_Elec_nHits->at(ite), T_Elec_passConversionVeto->at(ite));

	POGisolation = giveThePOGiso(pt, T_Elec_neutralHadronIso->at(ite), T_Elec_chargedHadronIso->at(ite), T_Elec_photonIso->at(ite), SCeta, T_Event_RhoIso);
    POGisolation2012A = giveThePOGisoPU2012(pt, T_Elec_neutralHadronIso->at(ite), T_Elec_chargedHadronIso->at(ite), T_Elec_photonIso->at(ite), SCeta, T_Event_RhoIso);
    POGnoRho = giveThePOGiso(pt, T_Elec_neutralHadronIso->at(ite), T_Elec_chargedHadronIso->at(ite), T_Elec_photonIso->at(ite), SCeta, 0);
	myTightID = 0;
	if (pt>20) myTightID=(POGisolation<0.15);
	else {
		if (T_Elec_isEB->at(ite)) myTightID=(POGisolation<0.15);
		else myTightID=(POGisolation<0.10);
	}
	myTightIDiso = ((myTightID&&passTightIdNoIsoPart)? 1 :0);
	nbRecoVertex = T_Vertex_z->size();
    if (isMC) PDGid = T_Gen_Elec_PDGid->at(ite);
    TLorentzVector *monElec = new TLorentzVector(T_Elec_Px->at(ite),T_Elec_Py->at(ite),T_Elec_Pz->at(ite),T_Elec_Energy->at(ite));
   /* if ((T_Elec_isEB->at(ite))&&(T_Elec_isFO->at(ite))){
        cout << "event = " << T_Event_EventNumber << endl;
        cout << "on a un elec dans EB" << endl;
        float neutralIso = calcPFIso(monElec, 0.3, neutralHadron,T_Elec_isEB->at(ite));
        if (T_Elec_neutralHadronIso->at(ite) != neutralIso) cout << "NEUTRAL ISO = " << "isoPOG="<< T_Elec_neutralHadronIso->at(ite) << "  myIso=" << neutralIso << endl;
        float gammaIso = calcPFIso(monElec,0.3, gamma,T_Elec_isEB->at(ite));
        if (T_Elec_photonIso->at(ite)!=gammaIso) cout << "GAMMA ISO = " << "isoPOG="<< T_Elec_photonIso->at(ite) << "  myIso=" << gammaIso << endl;
        float chargedIso = calcPFIso(monElec,0.3,chargedHadron,T_Elec_isEB->at(ite));
        if ((T_Elec_chargedHadronIso->at(ite)-chargedIso)>0.0001) cout << "CHARGED ISO = " << "isoPOG="<< T_Elec_chargedHadronIso->at(ite) << "  myIso=" << chargedIso << endl;
    }*/
    myPOGisolation = giveMyPOGiso(T_Elec_Pt->at(ite), monElec, 0.3, T_Elec_isEB->at(ite), SCeta, T_Event_RhoIso);
    PFisoDeltaBeta = PFisolationWithDeltaBeta(ite);
    radIsoStandard = giveRadIso(monElec, 0.01, T_Elec_isEB->at(ite),  true,  false);
    radIsoNoInner = giveRadIso(monElec, 0.0, T_Elec_isEB->at(ite),  true,  false);
    radisoNoThreshold = giveRadIso(monElec, 0.01, T_Elec_isEB->at(ite),  false,  false);
    radisoNoThresholdNoInner = giveRadIso(monElec, 0.0, T_Elec_isEB->at(ite),  false,  false);
    radIsoBeta = giveRadIso(monElec, 0.01, T_Elec_isEB->at(ite),  true,  true);
    //FillThePFtree( monElec,  absSCeta,  localType,  T_Elec_isEB->at(ite));
    
    radisoFonc0 = giveRadIsoFonc(monElec, 0.0, T_Elec_isEB->at(ite),  false,  false, 0);
    radisoFonc1 = giveRadIsoFonc(monElec, 0.0, T_Elec_isEB->at(ite),  false,  false, 1);
    radisoFonc2 = giveRadIsoFonc(monElec, 0.0, T_Elec_isEB->at(ite),  false,  false, 2);
    radisoFonc3 = giveRadIsoFonc(monElec, 0.0, T_Elec_isEB->at(ite),  false,  false, 3);

    
    chRadisoFonc0 = ((calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 0)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 0)/pt : 0.3);
    chRadisoFonc1 = ((calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 1)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 1)/pt : 0.3);
    chRadisoFonc2 = ((calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 2)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 2)/pt : 0.3);
    chRadisoFonc3 = ((calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 3)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, chargedHadron, T_Elec_isEB->at(ite), false, false, 3)/pt : 0.3);

    nhRadisoFonc0 = ((calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 0)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 0)/pt : 0.3);
    nhRadisoFonc1 = ((calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 1)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 1)/pt : 0.3);
    nhRadisoFonc2 = ((calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 2)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 2)/pt : 0.3);
    nhRadisoFonc3 = ((calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 3)/pt) < 0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, neutralHadron, T_Elec_isEB->at(ite), false, false, 3)/pt : 0.3);

    
    gRadisoFonc0 = ((calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 0)/pt)<0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 0)/pt : 0.3);
    gRadisoFonc1 = ((calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 1)/pt)<0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 1)/pt : 0.3);
    gRadisoFonc2 = ((calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 2)/pt)<0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 2)/pt : 0.3);
    gRadisoFonc3 = ((calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 3)/pt)<0.3 ? calcPFRadIsoFonc(monElec, 0.3, 0, gamma, T_Elec_isEB->at(ite), false, false, 3)/pt : 0.3);

    
  ///cou  cout << "radisoFonc0=" << radisoFonc0 << " radisoFonc1=" << radisoFonc1 << " radisoFonc2=" << radisoFonc2 << " radisoFonc3=" << radisoFonc3 << endl;
    
    isPair = T_Event_EventNumber%2;
    
    chIso0 = T_Elec_ChargedIso_DR0p1To0p1_DR0p0To0p1->at(ite);
    chIso1 = T_Elec_ChargedIso_DR0p1To0p2->at(ite);
    chIso2 = T_Elec_ChargedIso_DR0p2To0p3->at(ite);
    chIso3 = T_Elec_ChargedIso_DR0p3To0p4->at(ite);
    chIso4 = T_Elec_ChargedIso_DR0p4To0p5->at(ite);
    
    nhIso0 = T_Elec_NeutralHadronIso_DR0p0To0p1->at(ite);
    nhIso1 = T_Elec_NeutralHadronIso_DR0p1To0p2->at(ite);
    nhIso2 = T_Elec_NeutralHadronIso_DR0p2To0p3->at(ite);
    nhIso3 = T_Elec_NeutralHadronIso_DR0p3To0p4->at(ite);
    nhIso4 = T_Elec_NeutralHadronIso_DR0p4To0p5->at(ite);
    
    gIso0 = T_Elec_GammaIso_DR0p0To0p1->at(ite);
    gIso1 = T_Elec_GammaIso_DR0p1To0p2->at(ite);
    gIso2 = T_Elec_GammaIso_DR0p2To0p3->at(ite);
    gIso3 = T_Elec_GammaIso_DR0p3To0p4->at(ite);
    gIso4 = T_Elec_GammaIso_DR0p4To0p5->at(ite);
    
    
    
    if (passFO) testTree->Fill();
}

bool HWWAnalysis::passTightIdNoIso(int isEB, float dEtaIn, float dPhiIn, float sigmaIEtaIEta, float hoe, float ooemoop, float d0vtx, float dzvtx, int mHit, int convVeto){
    float cut_dEtaIn;
    float cut_dPhiIn;
    float cut_sigmaIEtaIEta;
    float cut_hoe;
    float cut_ooemoop;
    float cut_d0vtx;
    float cut_dzvtx;
    int   cut_vtxFit;
    int   cut_mHits;
    
    if (isEB){
        cut_dEtaIn        = 0.004; 
        cut_dPhiIn        = 0.030; 
        cut_sigmaIEtaIEta = 0.010;
        cut_hoe           = 0.120;
        cut_ooemoop       = 0.050;
        cut_d0vtx         = 0.020;
        cut_dzvtx         = 0.100;
        cut_vtxFit        = 1;
        cut_mHits         = 0;
    }
    else {
        cut_dEtaIn        = 0.005;
        cut_dPhiIn        = 0.020;
        cut_sigmaIEtaIEta = 0.030;
        cut_hoe           = 0.100;
        cut_ooemoop       = 0.050;
        cut_d0vtx         = 0.020;
        cut_dzvtx         = 0.100;
        cut_vtxFit        = 1;
        cut_mHits         = 0;
    }
    if (fabs(dEtaIn)>cut_dEtaIn) return false;
    if (fabs(dPhiIn)>cut_dPhiIn) return false;
  if (sigmaIEtaIEta> cut_sigmaIEtaIEta) return false;
    if (hoe>cut_hoe) return false;
    if (fabs(ooemoop)> cut_ooemoop) return false;
    if (fabs(d0vtx) > cut_d0vtx) return false;
    if (fabs(dzvtx) > cut_dzvtx) return false;
    if (mHit > cut_mHits) return false;
    if (cut_vtxFit&&(convVeto)) return false;
    
    return true;
}

bool HWWAnalysis::passTightIdNoIsoModified(int isEB, float dEtaIn, float dPhiIn, float sigmaIEtaIEta, float hoe, float ooemoop, float d0vtx, float dzvtx, int mHit, int convVeto){
    float cut_dEtaIn;
    float cut_dPhiIn;
    float cut_sigmaIEtaIEta;
    float cut_hoe;
    float cut_ooemoop;
    float cut_d0vtx;
    float cut_dzvtx;
    int   cut_vtxFit;
    int   cut_mHits;
    
    if (isEB){
        cut_dEtaIn        = 0.004; 
        cut_dPhiIn        = 0.030; 
        cut_sigmaIEtaIEta = 0.010;
        cut_hoe           = 0.120;
        cut_ooemoop       = 0.050;
        cut_d0vtx         = 0.020;
        cut_dzvtx         = 0.100;
        cut_vtxFit        = 1;
        cut_mHits         = 0;
    }
    else {
        cut_dEtaIn        = 0.005;
        cut_dPhiIn        = 0.020;
        cut_sigmaIEtaIEta = 0.030;
        cut_hoe           = 0.100;
        cut_ooemoop       = 0.050;
        cut_d0vtx         = 0.020;
        cut_dzvtx         = 0.100;
        cut_vtxFit        = 1;
        cut_mHits         = 0;
    }
    if (fabs(dEtaIn)>cut_dEtaIn) return false;
    if (fabs(dPhiIn)>cut_dPhiIn) return false;
  if (sigmaIEtaIEta> cut_sigmaIEtaIEta) return false;
    if (hoe>cut_hoe) return false;
    if (fabs(ooemoop)> cut_ooemoop) return false;
    //if (fabs(d0vtx) > cut_d0vtx) return false;
    //if (fabs(dzvtx) > cut_dzvtx) return false;
    if (mHit > cut_mHits) return false;
    //if (cut_vtxFit&&(convVeto)) return false;
    
    return true;
}




float HWWAnalysis::giveThePOGiso(float thePt, float isoNeutralHadrons, float isoCharged, float isoEm, float theSCEta, float rho){
	float EffectiveArea = 0;
	if (fabs(theSCEta) >= 0.0 && fabs(theSCEta) < 1.0 ) EffectiveArea = 0.100;
	if (fabs(theSCEta) >= 1.0 && fabs(theSCEta) < 1.479 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 1.479 && fabs(theSCEta) < 2.0 ) EffectiveArea = 0.085;
	if (fabs(theSCEta) >= 2.0 && fabs(theSCEta) < 2.2 ) EffectiveArea = 0.110;
	if (fabs(theSCEta) >= 2.2 && fabs(theSCEta) < 2.3 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 2.3 && fabs(theSCEta) < 2.4 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 2.4) EffectiveArea = 0.130;

	// apply to neutrals
    double rhoPrime = (rho>0 ? rho : 0);
    double iso_n = std::max(isoNeutralHadrons + isoEm - rhoPrime * EffectiveArea
							, 0.0);
	
    // compute final isolation
    double iso = (iso_n + isoCharged) / thePt;
	
	return iso; 

}
float HWWAnalysis::giveThePOGisoPU2012(float thePt, float isoNeutralHadrons, float isoCharged, float isoEm, float theSCEta, float rho){
	float EffectiveArea = 0;
	if (fabs(theSCEta) >= 0.0 && fabs(theSCEta) < 1.0 ) EffectiveArea = (0.122+0.013);
	if (fabs(theSCEta) >= 1.0 && fabs(theSCEta) < 1.479 ) EffectiveArea = (0.147+0.147);
	if (fabs(theSCEta) >= 1.479 && fabs(theSCEta) < 2.0 ) EffectiveArea = (0.055+0.010);
	if (fabs(theSCEta) >= 2.0 && fabs(theSCEta) < 2.2 ) EffectiveArea = (0.106+0.010);
	if (fabs(theSCEta) >= 2.2 && fabs(theSCEta) < 2.3 ) EffectiveArea = (0.138+ 0.024);
	if (fabs(theSCEta) >= 2.3 && fabs(theSCEta) < 2.4 ) EffectiveArea = (0.221+0.020);
	if (fabs(theSCEta) >= 2.4) EffectiveArea = 0.019;
    
	// apply to neutrals
    double rhoPrime = (rho>0 ? rho : 0);
    double iso_n = std::max(isoNeutralHadrons + isoEm - rhoPrime * EffectiveArea
							, 0.0);
	
    // compute final isolation
    double iso = (iso_n + isoCharged) / thePt;
	
	return iso;
    
}
float HWWAnalysis::giveMyPOGiso(float thePt, TLorentzVector *theElec,float deltaRmax, bool barrelEndcap, float theSCEta, float rho){
	float EffectiveArea = 0;
	if (fabs(theSCEta) >= 0.0 && fabs(theSCEta) < 1.0 ) EffectiveArea = 0.100;
	if (fabs(theSCEta) >= 1.0 && fabs(theSCEta) < 1.479 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 1.479 && fabs(theSCEta) < 2.0 ) EffectiveArea = 0.085;
	if (fabs(theSCEta) >= 2.0 && fabs(theSCEta) < 2.2 ) EffectiveArea = 0.110;
	if (fabs(theSCEta) >= 2.2 && fabs(theSCEta) < 2.3 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 2.3 && fabs(theSCEta) < 2.4 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 2.4) EffectiveArea = 0.211;
    
    // obtain the isolations
    float isoNeutralHadrons = calcPFIso(theElec,deltaRmax,neutralHadron,barrelEndcap);
    float isoEm = calcPFIso(theElec,deltaRmax,gamma,barrelEndcap);
    float isoCharged = calcPFIso(theElec,deltaRmax,chargedHadron,barrelEndcap);
    
	// apply to neutrals
    double rhoPrime = (rho>0 ? rho : 0);
    double iso_n = std::max(isoNeutralHadrons + isoEm - rhoPrime * EffectiveArea
							, 0.0);
	
    // compute final isolation
    double iso = (iso_n + isoCharged) / thePt;
	
	return iso;
    
}


float HWWAnalysis::calc03Iso(float thePt,  float theSCEta, float rho,int j){
	float EffectiveArea = 0;
	if (fabs(theSCEta) >= 0.0 && fabs(theSCEta) < 1.0 ) EffectiveArea = 0.100;
	if (fabs(theSCEta) >= 1.0 && fabs(theSCEta) < 1.479 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 1.479 && fabs(theSCEta) < 2.0 ) EffectiveArea = 0.085;
	if (fabs(theSCEta) >= 2.0 && fabs(theSCEta) < 2.2 ) EffectiveArea = 0.110;
	if (fabs(theSCEta) >= 2.2 && fabs(theSCEta) < 2.3 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 2.3 && fabs(theSCEta) < 2.4 ) EffectiveArea = 0.120;
	if (fabs(theSCEta) >= 2.4) EffectiveArea = 0.211;
    
    // obtain the isolations
    float isoNeutralHadrons =T_Elec_neutralHadronIso->at(j);
    float isoEm = T_Elec_photonIso->at(j);
    float isoCharged = T_Elec_chargedHadronIso->at(j);
    //cout << "Neutral=" << isoNeutralHadrons << " EM=" << isoEm << " Charged=" << isoCharged << endl;
	// apply to neutrals
    double rhoPrime = (rho>0 ? rho : 0);
    double iso_n = std::max(isoNeutralHadrons + isoEm - rhoPrime * EffectiveArea
							, 0.0);
	
    // compute final isolation
    double iso = (iso_n + isoCharged) / thePt;
	
	return iso;
    
}

float HWWAnalysis::calcPFIso(TLorentzVector *theElec,float deltaRmax, PFisolationType isoType, bool barrelEndcap){
 //   cout << "eta = " << eta << "phi = " << phi << endl;
    int nbOfPF = T_PF_Et->size();
    //cout << "nbOf PF particles=" << nbOfPF << endl;
    //cout << "eta " << theElec->Eta() << endl;
    float ptInCone = 0;
    int thePDGID = 0;
    if (isoType==neutralHadron) thePDGID=130;
    else if (isoType==gamma) thePDGID=22;
    else thePDGID = 211;
    for (int i = 0 ; i < nbOfPF ; i++){
        float deltaR = sqrt(pow(T_PF_Eta->at(i)-theElec->Eta(),2)+ pow(acos(cos(T_PF_Phi->at(i)-theElec->Phi())),2)) ;
//        if (isoType==neutralHadron) cout << "deltaR =" << deltaR << " pt=" << T_PF_Pt->at(i) << " PartId=" << T_PF_pdgID->at(i) <<  "issharingTrack=" << T_PF_hasTrack->at(i) << " isPU=" << T_PF_isPU->at(i)<< endl;
        
        if ((deltaR<deltaRmax)&&(fabs(T_PF_pdgID->at(i))==thePDGID)){
            if ((isoType==gamma)&&(!(barrelEndcap))&&(deltaR<0.08)) continue;
            if ((isoType==chargedHadron)&&(!(barrelEndcap))&&(deltaR<0.015)) continue;
            if ((!(isoType==chargedHadron))||(T_PF_isPU->at(i)==0)){//if charged particle, only from vertex
                // cout << "deltaR =" << deltaR << " pt=" << T_PF_Pt->at(i) << " PartId=" << T_PF_pdgID->at(i) <<  "issharingTrack=" << T_PF_hasTrack->at(i) << endl;
                ptInCone+= T_PF_Pt->at(i);
            }
        }
    }
    return ptInCone;
}

float HWWAnalysis::calcPFRadIso(TLorentzVector *theElec,float deltaRmax, float deltaRmin, PFisolationType isoType, bool barrelEndcap, bool cutForNonCharged, bool PUenergy){
    //   cout << "eta = " << eta << "phi = " << phi << endl;
    int nbOfPF = T_PF_Et->size();
    //cout << "nbOf PF particles=" << nbOfPF << endl;
    //cout << "eta " << theElec->Eta() << endl;
    float ptInCone = 0;
    int thePDGID = 0;
    if (isoType==neutralHadron) thePDGID=130;
    else if (isoType==gamma) thePDGID=22;
    else thePDGID = 211;
    for (int i = 0 ; i < nbOfPF ; i++){
        float deltaR = sqrt(pow(T_PF_Eta->at(i)-theElec->Eta(),2)+ pow(acos(cos(T_PF_Phi->at(i)-theElec->Phi())),2)) ;
        //        if (isoType==neutralHadron) cout << "deltaR =" << deltaR << " pt=" << T_PF_Pt->at(i) << " PartId=" << T_PF_pdgID->at(i) <<  "issharingTrack=" << T_PF_hasTrack->at(i) << " isPU=" << T_PF_isPU->at(i)<< endl;
        
        if ((deltaR<deltaRmax)&&(deltaR>deltaRmin)&&(fabs(T_PF_pdgID->at(i))==thePDGID)){
            if ((isoType==gamma)&&(!(barrelEndcap))&&(deltaR<0.08)) continue;
            if ((isoType==chargedHadron)&&(!(barrelEndcap))&&(deltaR<0.015)) continue;
            if (cutForNonCharged&&((!(isoType==chargedHadron))&&(T_PF_Pt->at(i)<1.0))) continue;
            if ((!(isoType==chargedHadron))||(((!(PUenergy))&&(T_PF_isPU->at(i)==0))||((PUenergy)&&(T_PF_isPU->at(i)==1)))){//if charged particle, only from vertex
                // cout << "deltaR =" << deltaR << " pt=" << T_PF_Pt->at(i) << " PartId=" << T_PF_pdgID->at(i) <<  "issharingTrack=" << T_PF_hasTrack->at(i) << endl;
                ptInCone+= T_PF_Pt->at(i)*(1 - 3*deltaR);
            }
        }
    }
    return ptInCone;
}



float HWWAnalysis::giveRadIso(TLorentzVector *theElec, float innerCone, bool barrelEndcap, bool cutForNonCharged, bool doDeltaBeta){
    double theNeutralIso = calcPFRadIso(theElec, 0.3, innerCone, neutralHadron, barrelEndcap, cutForNonCharged, false);
    double theGammaIso = calcPFRadIso(theElec, 0.3, innerCone, gamma, barrelEndcap, cutForNonCharged, false);
    double theChargedIso = calcPFRadIso(theElec, 0.3, innerCone, chargedHadron, barrelEndcap, cutForNonCharged, false);
    double theChargedIsoPU = calcPFRadIso(theElec, 0.3, innerCone, chargedHadron, barrelEndcap, cutForNonCharged, true);
    
  /*  cout << "neutral iso = " << theNeutralIso;
    cout << " charged iso = " << theChargedIso;
    cout << " gamma iso = " << theGammaIso;
    cout << " charged iso PU = " << theChargedIsoPU << endl;*/
    
    float allNeutralIso = std::max(theNeutralIso+theGammaIso,0.0);
    if (doDeltaBeta) allNeutralIso = std::max(theNeutralIso+theGammaIso-0.5*theChargedIsoPU,0.0);
    
    float allRelatIso = 1.0*(allNeutralIso + theChargedIso)/theElec->Pt();
    
    return allRelatIso;
    

}


float HWWAnalysis::calcPFRadIsoFonc(TLorentzVector *theElec,float deltaRmax, float deltaRmin, PFisolationType isoType, bool barrelEndcap, bool cutForNonCharged, bool PUenergy, int theFonct){
    //   cout << "eta = " << eta << "phi = " << phi << endl;
    int nbOfPF = T_PF_Et->size();
    //cout << "nbOf PF particles=" << nbOfPF << endl;
    //cout << "eta " << theElec->Eta() << endl;
    float ptInCone = 0;
    int thePDGID = 0;
    if (isoType==neutralHadron) thePDGID=130;
    else if (isoType==gamma) thePDGID=22;
    else thePDGID = 211;
    for (int i = 0 ; i < nbOfPF ; i++){
        float deltaR = sqrt(pow(T_PF_Eta->at(i)-theElec->Eta(),2)+ pow(acos(cos(T_PF_Phi->at(i)-theElec->Phi())),2)) ;
        //        if (isoType==neutralHadron) cout << "deltaR =" << deltaR << " pt=" << T_PF_Pt->at(i) << " PartId=" << T_PF_pdgID->at(i) <<  "issharingTrack=" << T_PF_hasTrack->at(i) << " isPU=" << T_PF_isPU->at(i)<< endl;
        
        if ((deltaR<deltaRmax)&&(deltaR>deltaRmin)&&(fabs(T_PF_pdgID->at(i))==thePDGID)){
            if ((isoType==gamma)&&(!(barrelEndcap))&&(deltaR<0.08)) continue;
            if ((isoType==chargedHadron)&&(!(barrelEndcap))&&(deltaR<0.015)) continue;
            if (cutForNonCharged&&((!(isoType==chargedHadron))&&(T_PF_Pt->at(i)<1.0))) continue;
            if ((!(isoType==chargedHadron))||(((!(PUenergy))&&(T_PF_isPU->at(i)==0))||((PUenergy)&&(T_PF_isPU->at(i)==1)))){//if charged particle, only from vertex
                // cout << "deltaR =" << deltaR << " pt=" << T_PF_Pt->at(i) << " PartId=" << T_PF_pdgID->at(i) <<  "issharingTrack=" << T_PF_hasTrack->at(i) << endl;
                deltaR = deltaR/0.3;
                if (theFonct==0){
                    ptInCone+= T_PF_Pt->at(i)*deltaR*deltaR*deltaR;
                }
                else if (theFonct==1){
                    ptInCone+= T_PF_Pt->at(i)*3*deltaR*deltaR*(1-deltaR);
                }
                else if (theFonct==2){
                    ptInCone+= T_PF_Pt->at(i)*3*deltaR*(1-deltaR)*(1-deltaR);
                }
                else {
                    ptInCone+= T_PF_Pt->at(i)*(1-deltaR)*(1-deltaR)*(1-deltaR);
                }
            }
        }
    }
  //  cout << "ptInCone=" << ptInCone << endl;
    return ptInCone;
    
}



float HWWAnalysis::giveRadIsoFonc(TLorentzVector *theElec, float innerCone, bool barrelEndcap, bool cutForNonCharged, bool doDeltaBeta, int theFonct){
    double theNeutralIso = calcPFRadIsoFonc(theElec, 0.3, innerCone, neutralHadron, barrelEndcap, cutForNonCharged, false, theFonct);
    double theGammaIso = calcPFRadIsoFonc(theElec, 0.3, innerCone, gamma, barrelEndcap, cutForNonCharged, false, theFonct);
    double theChargedIso = calcPFRadIsoFonc(theElec, 0.3, innerCone, chargedHadron, barrelEndcap, cutForNonCharged, false, theFonct);
    double theChargedIsoPU = calcPFRadIsoFonc(theElec, 0.3, innerCone, chargedHadron, barrelEndcap, cutForNonCharged, true, theFonct);
    
    /*  cout << "neutral iso = " << theNeutralIso;
     cout << " charged iso = " << theChargedIso;
     cout << " gamma iso = " << theGammaIso;
     cout << " charged iso PU = " << theChargedIsoPU << endl;*/
    
    float allNeutralIso = std::max(theNeutralIso+theGammaIso,0.0);
    if (doDeltaBeta) allNeutralIso = std::max(theNeutralIso+theGammaIso-0.5*theChargedIsoPU,0.0);
    
    float allRelatIso = 1.0*(allNeutralIso + theChargedIso)/theElec->Pt();
    
  //  cout << "allRelatIso=" << allRelatIso << endl;
    if (allRelatIso>0.3) return 0.3;
    else return allRelatIso;
    
}

float HWWAnalysis::PFisolationWithDeltaBeta(int theIte){
    
    float isoWithDeltaBeta = T_Elec_chargedHadronIso->at(theIte) + std::max( T_Elec_neutralHadronIso->at(theIte) + T_Elec_photonIso->at(theIte) - 1.0*T_Elec_puChargedHadronIso->at(theIte)/2, 0.0);
    float relatIso = isoWithDeltaBeta / T_Elec_Pt->at(theIte);
    
    return relatIso;
}

void HWWAnalysis::FillThePFtree(TLorentzVector *theElec, float elecSCeta, bool sigBg, bool barrelEndcap){
    int nbOfPF = T_PF_Et->size();

    for (int i = 0 ; i < nbOfPF ; i++){
        float deltaR = sqrt(pow(T_PF_Eta->at(i)-theElec->Eta(),2)+ pow(acos(cos(T_PF_Phi->at(i)-theElec->Phi())),2)) ;        
        if ((deltaR<0.4)){
            if ((fabs(T_PF_pdgID->at(i))==22)&&(!(barrelEndcap))&&(deltaR<0.08)) continue;
            if ((fabs(T_PF_pdgID->at(i))==211)&&(!(barrelEndcap))&&(deltaR<0.015)) continue;
            mainEta = theElec->Eta();
            mainPhi = theElec->Phi();
            mainEnergy = theElec->E();
            mainPt = theElec->Pt();
            mainSCEta = elecSCeta;
            PFdeltaEta = T_PF_Eta->at(i) - mainEta;
            PFdeltaPhi = congPhi(T_PF_Phi->at(i) - mainPhi);
            PFet = T_PF_Pt->at(i);
            PFtype = T_PF_pdgID->at(i);
            PFisPU  = T_PF_isPU->at(i);
            isSignalBg = sigBg;
            PFdeltaR = deltaR;
            isoInfo->Fill();
        }
    }
    
}


float HWWAnalysis::congPhi(float thePhi){
    if (fabs(thePhi)<TMath::Pi()) return thePhi;
    else if (thePhi>TMath::Pi()) return thePhi-2*TMath::Pi();
    else if (thePhi<(-1*TMath::Pi())) return thePhi+2*TMath::Pi();
    return 0;
}

bool HWWAnalysis::passPreCuts(float myPt, int isEB, float theSigIeta, float deltaEtaSC, float deltaPhiSC,float HOE)
{

    if (isEB) {
        if (theSigIeta                               > 0.01)  return false;
        if (fabs(deltaEtaSC)        > 0.007) return false;
        if (fabs(deltaPhiSC)        > 0.15)  return false;
        if (HOE                              > 0.12)  return false;
    } else {
        if (theSigIeta                               > 0.03)  return false;
        if (fabs(deltaEtaSC)        > 0.009) return false;
        if (fabs(deltaPhiSC)        > 0.10)  return false;
        if (HOE                              > 0.10)  return false;
    }
    
    return true;
    
}
bool HWWAnalysis::FOnoIso(float myPt, int isEB, float theSigIeta, float deltaEtaSC, float deltaPhiSC,float HOE,float d0, float dz,bool conversion, int nMissHit)
{
    
    if (isEB) {
        if (theSigIeta                               > 0.01)  return false;
        if (fabs(deltaEtaSC)        > 0.007) return false;
        if (fabs(deltaPhiSC)        > 0.15)  return false;
        if (HOE                              > 0.12)  return false;
    } else {
        if (theSigIeta                               > 0.03)  return false;
        if (fabs(deltaEtaSC)        > 0.009) return false;
        if (fabs(deltaPhiSC)        > 0.10)  return false;
        if (HOE                              > 0.10)  return false;
    }
    if (d0 > 0.02)                                  return false;
    if (dz > 0.1)                                   return false;
    if (nMissHit > 0)      return false;
    if (conversion)           return false;
    return true;
    
}


bool HWWAnalysis::passIPcuts(float d0, float dz)
{
    if (d0 > 0.02)                                  return false;
    if (dz > 0.1)                                   return false;
    return true;
}

bool HWWAnalysis::passMissItCons(bool conversion, int nMissHit)
{
    if (nMissHit > 0)      return false;
    if (conversion)           return false;
    return true;
}

void HWWAnalysis::Summary() {
  cout << " ---------------------------------------------------" << endl;
  InitialiseParameters();
    cout << "====================================================================================" << endl;
  cout << endl;
  
}




