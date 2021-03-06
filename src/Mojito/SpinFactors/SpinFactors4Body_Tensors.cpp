// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// created: 9 Oct 2013
#include "Mint/SpinFactors.h"
#include "Mint/SpinFactors4Body_Tensors.h"
#include "Mint/DecayTree.h"
#include "Mint/Utils.h"
#include "Mint/CLHEPSystemOfUnits.h"

#include "Mint/ZTspin1.h"
#include "Mint/ZTspin2.h"
#include "Mint/LeviCivita.h"
#include "Mint/SpinSumT.h"


using namespace std;
using namespace MINT;

DecayTree* SF_DtoVT_VtoP0P1_TtoP2P3_P::_exampleDecay=0;
DecayTree* SF_DtoVT_VtoP0P1_TtoP2P3_D::_exampleDecayD=0;

DecayTree* SF_DtoT1T2_T1toP0P1_T2toP2P3_S::_exampleDecay=0;
DecayTree* SF_DtoT1T2_T1toP0P1_T2toP2P3_P::_exampleDecay=0;
DecayTree* SF_DtoT1T2_T1toP0P1_T2toP2P3_D::_exampleDecay=0;

DecayTree* SF_DtoTP0_TtoVP1_VtoP2P3::_exampleDecay=0;
DecayTree* SF_DtoT1P0_T1toT2P1_T2toP2P3::_exampleDecay=0;

DecayTree* SF_DtoTS_TtoP0P1_StoP2P3::_exampleDecay=0;

DecayTree* SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::_exampleDecay=0;
DecayTree* SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::_exampleDecay=0;
DecayTree* SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::_exampleDecay=0;


DecayTree* SF_DtoAP0_AtoTP1_TtoP2P3::_exampleDecay=0;

//============================================================

// ----------- D->VT, various L states --------
const DecayTree& SF_DtoVT_VtoP0P1_TtoP2P3_P::getExampleDecay(){
  if(0==_exampleDecay){
    _exampleDecay = new DecayTree(421);
    // remark: addDgtr always returns a pointer to the 
    // last daughter that was added, thus allowing these
    // chains:
    _exampleDecay->addDgtr( 225)->addDgtr( 221, -211);
    _exampleDecay->addDgtr( 333)->addDgtr(-321,  321);
  }
  return *_exampleDecay;
}
const DecayTree& SF_DtoVT_VtoP0P1_TtoP2P3_P::exampleDecay(){
  return getExampleDecay();
}

const DecayTree& SF_DtoVT_VtoP0P1_TtoP2P3_D::getExampleDecay(){
  if(0==_exampleDecayD){
    _exampleDecayD = new DecayTree(421);
    // remark: addDgtr always returns a pointer to the 
    // last daughter that was added, thus allowing these
    // chains:
    _exampleDecayD->addDgtr( 225)->addDgtr( 221, -211);
    _exampleDecayD->addDgtr( 333)->addDgtr(-321,  321);

    _exampleDecayD->getVal().setL(2);
  }
  return *_exampleDecayD;
}
const DecayTree& SF_DtoVT_VtoP0P1_TtoP2P3_D::exampleDecay(){
  return getExampleDecay();
}

// ------- D-> TT, varous L states ----------
const DecayTree& SF_DtoT1T2_T1toP0P1_T2toP2P3_S::getExampleDecay(){
  // D->f2 f2 (with L=0)
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(225)->addDgtr(211, -211);
  _exampleDecay->addDgtr(225)->addDgtr(211, -211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoT1T2_T1toP0P1_T2toP2P3_S::exampleDecay(){
  return getExampleDecay();
}

const DecayTree& SF_DtoT1T2_T1toP0P1_T2toP2P3_P::getExampleDecay(){
  // D->f2 f2 (with L=1)
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(225)->addDgtr(211, -211);
  _exampleDecay->addDgtr(225)->addDgtr(211, -211);
  _exampleDecay->getVal().setL(1);
  return *_exampleDecay;
}
const DecayTree& SF_DtoT1T2_T1toP0P1_T2toP2P3_P::exampleDecay(){
  return getExampleDecay();
}

const DecayTree& SF_DtoT1T2_T1toP0P1_T2toP2P3_D::getExampleDecay(){
  // D->f2 f2 (with L=2)
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(225)->addDgtr(211, -211);
  _exampleDecay->addDgtr(225)->addDgtr(211, -211);
  _exampleDecay->getVal().setL(2);
  return *_exampleDecay;
}
const DecayTree& SF_DtoT1T2_T1toP0P1_T2toP2P3_D::exampleDecay(){
  return getExampleDecay();
}

// ------- D->TP0, various subesquent decays of the T ----

const DecayTree& SF_DtoTP0_TtoVP1_VtoP2P3::getExampleDecay(){
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(-211, 215)->addDgtr(211, 113)->addDgtr(-211, 211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoT1P0_T1toT2P1_T2toP2P3::exampleDecay(){
  return getExampleDecay();
}
const DecayTree& SF_DtoT1P0_T1toT2P1_T2toP2P3::getExampleDecay(){
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(-211, 100215)->addDgtr(211, 225)->addDgtr(-211, 211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoTP0_TtoVP1_VtoP2P3::exampleDecay(){
  return getExampleDecay();
}


// ----- D->TS ----------
const DecayTree& SF_DtoTS_TtoP0P1_StoP2P3::getExampleDecay(){
  if(0==_exampleDecay){
    _exampleDecay = new DecayTree(421);
    // remark: addDgtr always returns a pointer to the 
    // last daughter that was added, thus allowing these
    // chains:
    _exampleDecay->addDgtr( 225)->addDgtr( 221, -211);
    _exampleDecay->addDgtr( 9010221)->addDgtr(-321,  321);
  }
  return *_exampleDecay;
}
const DecayTree& SF_DtoTS_TtoP0P1_StoP2P3::exampleDecay(){
  return getExampleDecay();
}

// ----- D->PseudoTensor P0, various modes
const DecayTree& SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::getExampleDecay(){
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(-211, 10215)->addDgtr(211, 225)->addDgtr(211, -211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::exampleDecay(){
  return getExampleDecay();
}
const DecayTree& SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::getExampleDecay(){
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(-211, 10215)->addDgtr(211, 9000221)->addDgtr(211, -211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::exampleDecay(){
  return getExampleDecay();
}
const DecayTree& SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::getExampleDecay(){
  // (angular momenum of pi2->rho must be 1)
  _exampleDecay = new DecayTree(421);
  _exampleDecay->addDgtr(-211, 10215)->addDgtr(211, 113)->addDgtr(211, -211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::exampleDecay(){
  return getExampleDecay();
}

// ----- D->AP, A->TP, ---------
const DecayTree& SF_DtoAP0_AtoTP1_TtoP2P3::getExampleDecay(){
  _exampleDecay = new DecayTree(421);
  // example: D->a1(1260) pi, a1 -> f2(1270) pi 
  // (kinematically challenged, but in PDG - broad particles...)
  _exampleDecay->addDgtr(-211, 20213)->addDgtr(211, 225)->addDgtr(211, -211);
  return *_exampleDecay;
}
const DecayTree& SF_DtoAP0_AtoTP1_TtoP2P3::exampleDecay(){
  return getExampleDecay();
}

//=============================================================

bool SF_DtoVT_VtoP0P1_TtoP2P3_BASE::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoVT_VtoP0P1_TtoP2P3_S::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }

  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if(dgtr->getVal().SVPAT() == "V" && ! dgtr->isFinalState()) V = dgtr;
    else if(dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T = dgtr;
  }

  if(0==V || 0==T){
    cout << "ERROR in SF_DtoVT_VtoP0P1_TtoP2P3_BASE::parseTree"
	 << " Didn't find V or T " << V.get() << ", " << T.get() << endl;
    return false;
  }
  if(V->nDgtr() != 2){
    cout << "ERROR in SF_DtoVT_VtoP0P1_TtoP2P3_BASE::parseTree"
	 << " V should have 2 daughters, but it says it has "
	 << V->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[0] = V->getDgtrTreePtr(0);
  fsPS[1] = V->getDgtrTreePtr(1);
  //normalOrder(fsPS[0], fsPS[1]);

  if(T->nDgtr() != 2){
    cout << "ERROR in SF_DtoVT_VtoP0P1_TtoP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << T->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = T->getDgtrTreePtr(0);
  fsPS[3] = T->getDgtrTreePtr(1);
  //normalOrder(fsPS[2], fsPS[3]);

  // this->printYourself();
  return true;
}

double SF_DtoVT_VtoP0P1_TtoP2P3_P::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  const double units = GeV*GeV*GeV*GeV;
    
  TLorentzVector pV = p(0, evt) + p(1, evt);
  TLorentzVector qV = p(0, evt) - p(1, evt);
  TLorentzVector pT = p(2, evt) + p(3, evt);
  TLorentzVector qT = p(2, evt) - p(3, evt);
  TLorentzVector pD = pV + pT;
  TLorentzVector qD = pV - pT;
  
  double MV = mRes(V, evt);
  double MT = mRes(T, evt);
    
  ZTspin1 tV(qV, pV, MV);
  ZTspin2 tT(qT, pT, MT);  
    
  if (_useZemachTensors) {
      ZTspin1 tD(qD,pD,pD.M());
      return (tT.Contract(tV)).Dot(tD) / units;
  }
  
  double returnVal = (tT.Contract(tV)).Dot(pV) / units;

  if(dbThis){
    cout << " SF_DtoVT_VtoP0P1_TtoP2P3_P::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
  
}
void SF_DtoVT_VtoP0P1_TtoP2P3_P::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoVT_VtoP0P1_TtoP2P3_P"
     << "\n\t T1(V)_{\\mu}  T2(T)^{\\mu\\nu} pV_{\\nu}"
     << "\n\t implemented as: (tT.Contract(tV)).Dot(pV) / GeV^4"
     << "\n\t with T1(V) = ZTspin1 tV(qV, pV, MV), T2(T)= ZTspin2 tT(qT, pT, MT);"
     << "\n\t and pV = p(0) + p(1), qV = p(0) - p(1), pT = p(2) + p(3), qT = p(2) - p(3)"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}
// -------------------------------------

double SF_DtoVT_VtoP0P1_TtoP2P3_D::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pV = p(0, evt) + p(1, evt);
  TLorentzVector qV = p(0, evt) - p(1, evt);
  TLorentzVector pT = p(2, evt) + p(3, evt);
  TLorentzVector qT = p(2, evt) - p(3, evt);

  TLorentzVector pD = pT + pV;
  TLorentzVector qD = pT - qV;
  
  double MV = mRes(V, evt);
  double MT = mRes(T, evt);
  
  ZTspin1 tV(qV, pV, MV);
  ZTspin2 tT(qT, pT, MT);
  TLorentzVector vecT(tT.Contract(pV));

  const double units = GeV*GeV * GeV*GeV * GeV*GeV;
    
  if (_useZemachTensors) {
      ZTspin2 tD(qD,pD,pD.M());
      LorentzMatrix M = tD.Contract_1(tT);
      return LeviCivita(tV,pD,M) / units;
  }    

  double returnVal = LeviCivita(tV, vecT, qD, pD) /units;


  if(dbThis){
    cout << " SF_DtoVT_VtoP0P1_TtoP2P3_D::getVal "
	 << " returning " << returnVal
	 << endl;
    double checkVal = LeviCivita(p(0, evt), p(1, evt), p(2, evt), p(3, evt))*MV*MT/units;
    cout << "cross check: " << checkVal
	 << " ratio " << checkVal/returnVal << endl;
  }
  return returnVal;
  
}
void SF_DtoVT_VtoP0P1_TtoP2P3_D::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoVT_VtoP0P1_TtoP2P3_D"
     << "\n\t ZTspin1 tV(qV, pV, MV);"
     << "\n\t  ZTspin2 tT(qT, pT, MT);"
     << "\n\t   TLorentzVector vecT(tT.Contract(pV));"
     << "\n\t  const double units = GeV*GeV * GeV*GeV * GeV*GeV;"
     << "\n\t  double returnVal = LeviCivita(tV, vecT, qD, pD) /units;"
     << "\n\t with pV = p(0) + p(1), qV = p(0) - p(1), pT = p(2) + p(3), qT = p(2) - p(3)"
     << "\n\t and  pD = pT + pV, qD = pT - qV;"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

//------------------------------------------
// D->TT, various L states
//-------------------------------------------
// -----------------------------------------------
bool SF_DtoT1T2_T1toP0P1_T2toP2P3_BASE::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoT1T2_T1toP0P1_T2toP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }

  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if(dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T[i] = dgtr;
    else{
      cout << "ERROR in SF_DtoT1T2_T1toP0P1_T2toP2P3_BASE::parseTree"
	   << " expected exactly 2 Tensors, have "
	   << theDecay(pat).nDgtr();
      cout << "\n so daughter " << i << " is not a tensor." << endl;
      return false;
    }
  }

  if(0==T[0] || 0==T[1]){
    cout << "ERROR in SF_DtoT1T2_T1toP0P1_T2toP2P3_BASE::parseTree"
	 << " Didn't find T or S " << T[0].get() << ", " << T[1].get() << endl;
    return false;
  }
  if(T[0]->nDgtr() != 2){
    cout << "ERROR in SF_DtoT1T2_T1toP0P1_T2toP2P3_BASE::parseTree"
	 << " T[0] should have 2 daughters, but it says it has "
	 << T[0]->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[0] = T[0]->getDgtrTreePtr(0);
  fsPS[1] = T[0]->getDgtrTreePtr(1);
  // normalOrder(fsPS[0], fsPS[1]);

  if(T[1]->nDgtr() != 2){
    cout << "ERROR in SF_DtoT1T2_T1toP0P1_T2toP2P3_S::parseTree"
	 << " T[1] should have 2 daughters, but it says it has "
	 << T[1]->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = T[1]->getDgtrTreePtr(0);
  fsPS[3] = T[1]->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);

  // this->printYourself();
  return true;
}
// ------- same decay, s-wave -----
double SF_DtoT1T2_T1toP0P1_T2toP2P3_S::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pT1 = p(0, evt) + p(1, evt);
  TLorentzVector qT1 = p(0, evt) - p(1, evt);
  TLorentzVector pT2 = p(2, evt) + p(3, evt);
  TLorentzVector qT2 = p(2, evt) - p(3, evt);
  
  double MT1 = mRes(T[0], evt);
  double MT2 = mRes(T[1], evt);
  
  ZTspin2 tT1(qT1, pT1, MT1);
  ZTspin2 tT2(qT2, pT2, MT2);
  
  const double units = GeV*GeV*GeV*GeV;

  double returnVal = tT1.Contract_2(tT2) / units;

  if(dbThis){
    cout << " SF_DtoT1T2_T1toP0P1_T2toP2P3_S::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
  
}

void SF_DtoT1T2_T1toP0P1_T2toP2P3_S::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoT1T2_T1toP0P1_T2toP2P3_S"
     << "\n\t  ZTspin2 tT1(qT1, pT1, MT1);"
     << "\n\t  ZTspin2 tT2(qT2, pT2, MT2);"
     << "\n\t  return: tT1.Contract_2(tT2) / GeV^8"
     << "\n\t with pT1 = p(0) + p(1), qT = p(0) - p(1), pT2 = p(2) + p(3), qT2 = p(2)-p(3)"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}
// ------- same decay, P-wave -----
double SF_DtoT1T2_T1toP0P1_T2toP2P3_P::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pT1 = p(0, evt) + p(1, evt);
  TLorentzVector qT1 = p(0, evt) - p(1, evt);
  TLorentzVector pT2 = p(2, evt) + p(3, evt);
  TLorentzVector qT2 = p(2, evt) - p(3, evt);
  TLorentzVector pD  = pT1 + pT2;
  TLorentzVector qD  = pT1 - pT2;
  

  double MT1 = mRes(T[0], evt);
  double MT2 = mRes(T[1], evt);
  
  ZTspin2 tT1(qT1, pT1, MT1);
  ZTspin2 tT2(qT2, pT2, MT2);
  
  const double units = GeV*GeV * GeV*GeV * GeV*GeV ;
    
  if (_useZemachTensors) {
            ZTspin1 tD(qD,pD,pD.M());
            return LeviCivita(pD, tD, tT1.Contract_1(tT2)) / units;
  }  

  double returnVal = LeviCivita(pD, qD, tT1.Contract_1(tT2)) / units;

  if(dbThis){
    cout << " SF_DtoT1T2_T1toP0P1_T2toP2P3_P::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
  
}

void SF_DtoT1T2_T1toP0P1_T2toP2P3_P::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoT1T2_T1toP0P1_T2toP2P3_P"
     << "\n\t  ZTspin2 tT1(qT1, pT1, MT1);"
     << "\n\t  ZTspin2 tT2(qT2, pT2, MT2);"
     << "\n\t  pD = pT1 + pT2 "
     << "\n\t  qD = pT1 - pT2 "
     << "\n\t  return:  LeviCivita(pD, qD, tT1.Contract_1(tT2))/ GeV^4"
     << "\n\t  with pT1 = p(0) + p(1), qT = p(0) - p(1), pT2 = p(2) + p(3), qT2 = p(2)-p(3)"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}
// ------- same decay, D-wave -----
double SF_DtoT1T2_T1toP0P1_T2toP2P3_D::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pT1 = p(0, evt) + p(1, evt);
  TLorentzVector qT1 = p(0, evt) - p(1, evt);
  TLorentzVector pT2 = p(2, evt) + p(3, evt);
  TLorentzVector qT2 = p(2, evt) - p(3, evt);
  TLorentzVector qD  = pT1 - pT2;
  TLorentzVector pD  = pT1 + pT2;

  double MT1 = mRes(T[0], evt);
  double MT2 = mRes(T[1], evt);
  
  ZTspin2 tT1(qT1, pT1, MT1);
  ZTspin2 tT2(qT2, pT2, MT2);
  
  TLorentzVector X1(tT1.Contract(qD));
  TLorentzVector X2(tT2.Contract(qD));

  const double units = GeV*GeV*GeV * GeV*GeV*GeV;
    
  if (_useZemachTensors) {
            ZTspin1 tD(qD,pD,pD.M());
            double val = (tT1.Contract(tD)).Dot(tT2.Contract(tD));
            val -= 1./3. * tD.Dot(tD) * tT1.Contract_2(tT2);
            val += 1./3. * tD.Dot(tD) / (pD.Dot(pD)) * (tT1.Contract(pD)).Dot(tT2.Contract(pD));
            return val/units;
  }  

  double returnVal = X1.Dot(X2) / units;

  if(dbThis){
    cout << " SF_DtoT1T2_T1toP0P1_T2toP2P3_D::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
  
}

void SF_DtoT1T2_T1toP0P1_T2toP2P3_D::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoT1T2_T1toP0P1_T2toP2P3_D"
     << "\n\t  ZTspin2 tT1(qT1, pT1, MT1);"
     << "\n\t  ZTspin2 tT2(qT2, pT2, MT2);"
     << "\n\t  qD = pT1 - pT2"
     << "\n\t  return: (tT1.Contract(qD))* (tT2.Contract(qD))/ GeV^4"
     << "\n\t  with pT1 = p(0) + p(1), qT = p(0) - p(1), pT2 = p(2) + p(3), qT2 = p(2)-p(3)"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

// ------------------------------------------
// D->TP0, various subesquent decays of the T
//------------------------------------------------
bool SF_DtoTP0_TtoVP1_VtoP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoTP0_TtoVP1_VtoP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }
  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[0] = dgtr;
    else if(dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T = dgtr;
  }
  if(0==T || 0==fsPS[0]){
    cout << "ERROR in SF_DtoTP0_TtoVP1_VtoP2P3_BASE::parseTree"
	 << " Didn't find T or P0 " << T.get() << ", " << fsPS[0].get() << endl;
    return false;
  }
  
  if(T->nDgtr() != 2){
    cout << "ERROR in SF_DtoTP0_TtoVP1_VtoP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << T->nDgtr() << "."
	 << endl;
    return false;
  }
  for(int i=0; i< T->nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= T->getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "V" && ! dgtr->isFinalState()) V = dgtr;
    else if(dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[1] = dgtr;
  }
  if(0==V || 0==fsPS[1]){
    cout << "ERROR in SF_DtoTP0_TtoVP1_VtoP2P3_BASE::parseTree"
	 << " Didn't find V2 or P1 " << V.get() << ", " << fsPS[1].get() << endl;
    return false;
  }
  
  if(V->nDgtr() != 2){
    cout << "ERROR in SF_DtoTP0_TtoVP1_VtoP2P3_BASE::parseTree"
	 << " V should have 2 daughters, but it says it has "
	 << V->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = V->getDgtrTreePtr(0);
  fsPS[3] = V->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);
  
  // this->printYourself();
  return true;
}
// ------- note that T->VP1 has to be in D-wave for parity and ang mom ----
double SF_DtoTP0_TtoVP1_VtoP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pV = p(2, evt) + p(3, evt);
  TLorentzVector qV = p(2, evt) - p(3, evt);
  TLorentzVector pT = pV + p(1, evt);
  TLorentzVector qT = pV - p(1, evt);
  TLorentzVector pD = pT + p(0, evt);
  TLorentzVector qD = pT - p(0, evt);

  double MT = mRes(T, evt);
  double MV = mRes(V, evt);
  
  ZTspin2 t2T(qT, pT, MT);
  ZTspin1 tV(qV, pV, MV);
  
  TLorentzVector DT(t2T.Contract(qD));

  const double units = GeV*GeV * GeV*GeV * GeV*GeV;

  if (_useZemachTensors) {
            ZTspin1 tD(qD,pD,pD.M());
            ZTspin1 t1T(qT, pT, MT);
            SpinSumV P1T(pT,MT);
            SpinSumT P2T(pT,MT);
      
            TLorentzVector tmp = LeviCivita(t1T,pT,P1T.Dot(tV));
            
            return P2T.Sandwich(tD,tD,tmp,t1T)/ units + 1./3. * tD.Dot(tD) / (pT.Dot(pT)) * P2T.Sandwich(pD,pD,tmp,t1T)/ units;
  }  
    
  double returnVal = LeviCivita(pD, qD, DT, tV) / units;

  if(dbThis){
    cout << " SF_DtoTP0_TtoVP1_VtoP2P3::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
}

void SF_DtoTP0_TtoVP1_VtoP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoTP0_TtoVP1_VtoP2P3"
     << "\n\t  ZTspin2 tT1(qT, pT, MT);"
     << "\n\t  ZTspin1 tT2(qV, pV, MV);"
     << "\n\t  pV = p(2) + p(3); qV = p(2) - p(3);"
     << "\n\t  pT = pV + p(1); qT = pV - p(1);"
     << "\n\t  pD = pT + p(0); qD = pT - p(0);"
     << "\n\t  return: DT(tT.Contract(qD))/GeV^6"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}
// -----------------------------------------------
bool SF_DtoT1P0_T1toT2P1_T2toP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoT1P0_T1toT2P1_T2toP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }
  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[0] = dgtr;
    else if(dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T1 = dgtr;
  }
  if(0==T1 || 0==fsPS[0]){
    cout << "ERROR in SF_DtoT1P0_T1toT2P1_T2toP2P3_BASE::parseTree"
	 << " Didn't find T or P0 " << T1.get() << ", " << fsPS[0].get() << endl;
    return false;
  }
  
  if(T1->nDgtr() != 2){
    cout << "ERROR in SF_DtoT1P0_T1toT2P1_T2toP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << T1->nDgtr() << "."
	 << endl;
    return false;
  }
  for(int i=0; i< T1->nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= T1->getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T2 = dgtr;
    else if(dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[1] = dgtr;
  }
  if(0==T2 || 0==fsPS[1]){
    cout << "ERROR in SF_DtoT1P0_T1toT2P1_T2toP2P3_BASE::parseTree"
	 << " Didn't find T2 or P1 " << T2.get() << ", " << fsPS[1].get() << endl;
    return false;
  }
  
  if(T2->nDgtr() != 2){
    cout << "ERROR in SF_DtoT1P0_T1toT2P1_T2toP2P3_BASE::parseTree"
	 << " V should have 2 daughters, but it says it has "
	 << T2->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = T2->getDgtrTreePtr(0);
  fsPS[3] = T2->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);
  
  // this->printYourself();
  return true;
}
// ------- note that T->TP1 has to be in P-wave (or F, etc) for parity and ang mom ----
double SF_DtoT1P0_T1toT2P1_T2toP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pT2 = p(2, evt) + p(3, evt);
  TLorentzVector qT2 = p(2, evt) - p(3, evt);
  TLorentzVector pT1 = pT2 + p(1, evt);
  TLorentzVector qT1 = pT2 - p(1, evt);
  TLorentzVector pD = pT1 + p(0, evt);
  TLorentzVector qD = pT1 - p(0, evt);

  double MT1 = mRes(T1 ,evt);
  double MT2 = mRes(T2, evt);
  
  ZTspin2 tT1(qT1, pT1, MT1);
  ZTspin2 tT2(qT2, pT2, MT2);
  
  TLorentzVector DT(tT1.Contract(qT1));

  const double units = GeV*GeV * GeV*GeV * GeV*GeV;
    
  if (_useZemachTensors) {
      ZTspin1 LT1(qT1, pT1, MT1);
      ZTspin1 LT2(qT2, pT2, MT2);
      SpinSumT PT1(pT1,MT1);
      
      TLorentzVector tmp1 = LeviCivita(LT1,pT1,LT2);
      TLorentzVector tmp2 = LeviCivita(LT1,pT1,pT2);
      tmp2 *= LT2.Dot(LT2)*1./(3.*pT2.Dot(pT2));
            
      return (PT1.Sandwich(qD,qD,tmp1,LT2)+PT1.Sandwich(qD,qD,tmp2,pT2))/units;
  }  

  //returns always 0 since tT2 is a symmetric matrix ! 
  //should not be used !  
  double returnVal = LeviCivita(pD, DT, tT2) / units;   

  if(dbThis){
    cout << " SF_DtoT1P0_T1toT2P1_T2toP2P3::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
}

void SF_DtoT1P0_T1toT2P1_T2toP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoT1P0_T1toT2P1_T2toP2P3"
     << "\n\t  ZTspin2 tT1(qT, pT, MT);"
     << "\n\t  ZTspin1 tT2(qV, pV, MV);"
     << "\n\t    TLorentzVector DT(tT1.Contract(qT1));"
     << "\n\t  pV = p(2, evt) + p(3); qV = p(2) - p(3);"
     << "\n\t  pT = pV + p(1); qT = pV - p(1);"
     << "\n\t  pD = pT + p(0); qD = pT - p(0);"
     << "\n\t  return: LeviCivita(pD, DT, tT2)/GeV^6"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

// -------------------
//     D->TS
// ------------------
bool SF_DtoTS_TtoP0P1_StoP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoTS_TtoP0P1_StoP2P3::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }

  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if(dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T = dgtr;
    else if(dgtr->getVal().SVPAT() == "S" && ! dgtr->isFinalState()) S = dgtr;
  }

  if(0==S || 0==T){
    cout << "ERROR in SF_DtoTS_TtoP0P1_StoP2P3::parseTree"
	 << " Didn't find T or S " << T.get() << ", " << S.get() << endl;
    return false;
  }
  if(T->nDgtr() != 2){
    cout << "ERROR in SF_DtoTS_TtoP0P1_StoP2P3::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << T->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[0] = T->getDgtrTreePtr(0);
  fsPS[1] = T->getDgtrTreePtr(1);
  // normalOrder(fsPS[0], fsPS[1]);

  if(S->nDgtr() != 2){
    cout << "ERROR in SF_DtoTS_TtoP0P1_StoP2P3::parseTree"
	 << " S should have 2 daughters, but it says it has "
	 << S->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = S->getDgtrTreePtr(0);
  fsPS[3] = S->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);

  // this->printYourself();
  return true;
}


// -----------------------------------------------

double SF_DtoTS_TtoP0P1_StoP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pT = p(0, evt) + p(1, evt);
  TLorentzVector qT = p(0, evt) - p(1, evt);
  TLorentzVector pS = p(2, evt) + p(3, evt);
  TLorentzVector qD = pT - pS;
  TLorentzVector pD = pT + pS;
  
  double MT = mRes(T, evt);
  
  ZTspin2 tT(qT, pT, MT);
  
  const double units = GeV*GeV * GeV*GeV;

  if (_useZemachTensors) {
            ZTspin2 tD(qD,pD,pD.M());
            SpinSumT P2T(pT,MT);        
            SpinSumV P1T(pT,MT);        
            ZTspin1 t1T(qT, pT, MT);

            SymmLorentzMatrix qq(qD);
            //return tT.Contract_2(qq);
            if(dbThis){
                cout << "" << endl;
                cout << P2T.Sandwich(qD,qD,qT,qT)/units << endl;
                cout << tT.Contract_2(qq)/units << endl;
                cout << tT.Contract(qD).Dot(qD)/units << endl;
                cout << (qD.Dot(t1T)*qD.Dot(t1T) + 1./3. * t1T.Dot(t1T) * P1T.Sandwich(qD,qD))/units << endl;
                cout << "" << endl;
            }
            //return P2T.Sandwich(qD,qD,qT,qT)/units;
            return tT.Contract_2(tD)/units;
  }  
    
  double returnVal = (tT.Contract(qD)).Dot(qD) / units;

  if(dbThis){
      
      
    cout << " SF_DtoTS_TtoP0P1_StoP2P3::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
  
}

void SF_DtoTS_TtoP0P1_StoP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoTS_TtoP0P1_StoP2P3"
     << "\n\t  ZTspin2 tT(qT, pT, MT);"
     << "\n\t  return: (tT.Contract(qD)).Dot(qD) / GeV^4"
     << "\n\t with pT = p(0) + p(1), qT = p(0) - p(1), pS = p(2) + p(3), qD = pT - qV;"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

// ---------------
//    D->pseudoTensor P0, various modes
// --------------

// -----------------------------------------------
bool SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }
  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[0] = dgtr;
    else if(dgtr->getVal().SVPAT() == "PT" && ! dgtr->isFinalState()) PT = dgtr;
  }
  if(0==PT || 0==fsPS[0]){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3_BASE::parseTree"
	 << " Didn't find T or P0 " << PT.get() << ", " << fsPS[0].get() << endl;
    return false;
  }
  
  if(PT->nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << PT->nDgtr() << "."
	 << endl;
    return false;
  }
  for(int i=0; i< PT->nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= PT->getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "V" && ! dgtr->isFinalState()) V = dgtr;
    else if(dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[1] = dgtr;
  }
  if(0==V || 0==fsPS[1]){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3_BASE::parseTree"
	 << " Didn't find V2 or P1 " << V.get() << ", " << fsPS[1].get() << endl;
    return false;
  }
  
  if(V->nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3_BASE::parseTree"
	 << " V should have 2 daughters, but it says it has "
	 << V->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = V->getDgtrTreePtr(0);
  fsPS[3] = V->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);
  
  // this->printYourself();
  return true;
}
// ------- note that pseudoT->VP1 has to be in P-wave for parity and ang mom ----
double SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pV = p(2, evt) + p(3, evt);
  TLorentzVector qV = p(2, evt) - p(3, evt);
  TLorentzVector pT = pV + p(1, evt);
  TLorentzVector qT = pV - p(1, evt);
  TLorentzVector pD = pT + p(0, evt);
  TLorentzVector qD = pT - p(0, evt);

  double MT = mRes(PT, evt);
  double MV = mRes(V, evt);
  
  ZTspin1 tV(qV, pV, MV);
  ZTspin2 tT(qT, pT, MT);
  //SpinSumT ProT(pT, MT);
  
  const double units = GeV*GeV * GeV*GeV;
    
  if (_useZemachTensors) {
            SpinSumT PT(pT,MT);
            ZTspin1 LT(qT,pT,MT);
            ZTspin1 LD(qD,pD,pD.M());
            return PT.Sandwich(LD,LD,LT,tV)/units + 1./3. * LD.Dot(LD) / (pD.Dot(pD)) * PT.Sandwich(pD,pD,LT,tV)/units;
            //return PT.Sandwich(qD,qD,LT,tV)/units;
  }  

  double returnVal = (tT.Contract(qT)).Dot(tV) / units;

  if(dbThis){
    cout << " SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::getVal "
	 << " returning " << returnVal
	 << endl;
    SpinSumT ProT(pT, MT);
    cout << " compare to: "
	 << ProT.Sandwich(qD, qD, qT, tV) / units
	 << endl;
  }
  return returnVal;
}

void SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoPseudoTP0_PseudoTtoVP1_VtoP2P3"
     << "\n\t  ProT(pT, MT)"
     << "\n\t  ZTspin1 tT2(qV, pV, MV);"
     << "\n\t  pV = p(2) + p(3); qV = p(2) - p(3);"
     << "\n\t  pT = pV + p(1); qT = pV - p(1);"
     << "\n\t  pD = pT + p(0); qD = pT - p(0);"
     << "\n\t  return: ProT.Sandwich(qD, qD, qT, tV) / GeV^4"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}
// -----------------------------------------------
bool SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }
  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[0] = dgtr;
    else if(dgtr->getVal().SVPAT() == "PT" && ! dgtr->isFinalState()) PT = dgtr;
  }
  if(0==PT || 0==fsPS[0]){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3_BASE::parseTree"
	 << " Didn't find T or P0 " << PT.get() << ", " << fsPS[0].get() << endl;
    return false;
  }
  
  if(PT->nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << PT->nDgtr() << "."
	 << endl;
    return false;
  }
  for(int i=0; i< PT->nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= PT->getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "S" && ! dgtr->isFinalState()) S = dgtr;
    else if(dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[1] = dgtr;
  }
  if(0==S || 0==fsPS[1]){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3_BASE::parseTree"
	 << " Didn't find S or P1 " << S.get() << ", " << fsPS[1].get() << endl;
    return false;
  }
  
  if(S->nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3_BASE::parseTree"
	 << " V should have 2 daughters, but it says it has "
	 << S->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = S->getDgtrTreePtr(0);
  fsPS[3] = S->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);
  
  // this->printYourself();
  return true;
}
// ------- note that pseudoT->SP1 has to be in D-wave for parity and ang mom ----
double SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());
  
  TLorentzVector pS = p(2, evt) + p(3, evt);
  TLorentzVector qS = p(2, evt) - p(3, evt);
  TLorentzVector pT = pS + p(1, evt);
  TLorentzVector qT = pS - p(1, evt);
  TLorentzVector pD = pT + p(0, evt);
  TLorentzVector qD = pT - p(0, evt);

  double MT = mRes(PT, evt);
  //double MV = mRes(S);
  
  ZTspin2 tT(qT, pT, MT);
  
  const double units = GeV*GeV * GeV*GeV;
    
  if(_useZemachTensors){
      ZTspin2 tD(qD,pD,pD.M());
      return tD.Contract_2(tT)/units;
  }  

  double returnVal =  (tT.Contract(qD)).Dot(qD) / units;

  if(dbThis){
    cout << " SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
}

void SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoPseudoTP0_PseudoTtoSP1_StoP2P3"
     << "\n\t  ProT(pT, MT)"
     << "\n\t  ZTspin1 tT2(qV, pV, MV);"
     << "\n\t  pS = p(2) + p(3); qS = p(2) - p(3);"
     << "\n\t  pT = pV + p(1); qT = pV - p(1);"
     << "\n\t  pD = pT + p(0); qD = pT - p(0);"
     << "\n\t  (tT.Contract(qD)).Dot(qD) / GeV^4"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

// -------------------------------------------
bool SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }
  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[0] = dgtr;
    else if(dgtr->getVal().SVPAT() == "PT" && ! dgtr->isFinalState()) PT = dgtr;
  }
  if(0==PT || 0==fsPS[0]){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3_BASE::parseTree"
	 << " Didn't find T or P0 " << PT.get() << ", " << fsPS[0].get() << endl;
    return false;
  }
  
  if(PT->nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << PT->nDgtr() << "."
	 << endl;
    return false;
  }
  for(int i=0; i< PT->nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= PT->getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T = dgtr;
    else if(dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[1] = dgtr;
  }
  if(0==T || 0==fsPS[1]){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3_BASE::parseTree"
	 << " Didn't find T or P1 " << T.get() << ", " << fsPS[1].get() << endl;
    return false;
  }
  
  if(T->nDgtr() != 2){
    cout << "ERROR in SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3_BASE::parseTree"
	 << " V should have 2 daughters, but it says it has "
	 << T->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = T->getDgtrTreePtr(0);
  fsPS[3] = T->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);
  
  // this->printYourself();
  return true;
}
// ------- note that pseudoT->SP1 has to be in D-wave for parity and ang mom ----
double SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());

  double invGeV = 1;// 1.0/GeV;
  
  TLorentzVector pT = p(2, evt) + p(3, evt);
  TLorentzVector qT = p(2, evt) - p(3, evt);
  TLorentzVector pPT = pT + p(1, evt);
  TLorentzVector qPT = pT - p(1, evt);
  TLorentzVector pD = pPT + p(0, evt);
  TLorentzVector qD = pPT - p(0, evt);

  double MPT = mRes(PT, evt);
  double MT  = mRes(T , evt);
  //double MV = mRes(S);
  
  ZTspin2 tPT(qD*invGeV, pPT*invGeV, MPT*invGeV);
  ZTspin2 tT(qT*invGeV, pT*invGeV, MT*invGeV);

  const double units = GeV*GeV*GeV*GeV;
    
  if (_useZemachTensors) {
            SpinSumT PT1(pPT,MPT);
            ZTspin1 LT(qT,pT,MT);
      
            return (PT1.Sandwich(qD,qD,LT,LT) + 1./3. * LT.Dot(LT)/(pT.Dot(pT)) * PT1.Sandwich(qD,qD,pT,pT))/units;
  }

  double returnVal =  tPT.Contract_2(tT) / units;

  if(dbThis){
    cout << " SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::getVal "
	 << " returning " << returnVal
	 << " = " << tPT.Contract_2(tT) << " / " << units
	 << "\n MPT = " << MPT << ", MT = " << MT
	 << endl;
  }
  return returnVal;
}

void SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoPseudoTP0_PseudoTtoTP1_TtoP2P3"
     << "\n\t  ZTspin2 tPT(qD, pPT, MPT);"
     << "\n\t  ZTspin2 tT(qT, pT, MT);"
     << "\n\t  pT = p(2) + p(3); qT = p(2) - p(3);"
     << "\n\t  pPT = pT + p(1); qPT = pT - p(1);"
     << "\n\t  pD = pPT + p(0); qD = pPT - p(0);"
     << "\n\t  tPT.Contract_2(tT)/ GeV^8"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

// ==========================================
// -----------------------------------------------
// D->PA, A->TP (ang mom L=1 at each vtx)
// -----------------------------------------------

bool SF_DtoAP0_AtoTP1_TtoP2P3::parseTree(const DalitzEventPattern& pat){
  //  bool debugThis=false;
  if(fsPS.size() < 4) fsPS.reserve(4);
  if(theDecay(pat).nDgtr() != 2){
    cout << "ERROR in SF_DtoAP0_AtoTP1_TtoP2P3_BASE::parseTree"
	 << " expected exactly 2 daughers of D, have "
	 << theDecay(pat).nDgtr();
    return false;
  }
  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[0] = dgtr;
    else if(dgtr->getVal().SVPAT() == "A" && ! dgtr->isFinalState()) A = dgtr;
  }
  if(0==A || 0==fsPS[0]){
    cout << "ERROR in SF_DtoAP0_AtoTP1_TtoP2P3_BASE::parseTree"
	 << " Didn't find T or P0 " << A.get() << ", " << fsPS[0].get() << endl;
    return false;
  }
  
  if(A->nDgtr() != 2){
    cout << "ERROR in SF_DtoAP0_AtoTP1_TtoP2P3_BASE::parseTree"
	 << " T should have 2 daughters, but it says it has "
	 << A->nDgtr() << "."
	 << endl;
    return false;
  }
  for(int i=0; i< A->nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= A->getDgtrTreePtr(i);
    if     (dgtr->getVal().SVPAT() == "T" && ! dgtr->isFinalState()) T = dgtr;
    else if(dgtr->getVal().SVPAT() == "P" &&   dgtr->isFinalState()) fsPS[1] = dgtr;
  }
  if(0==T || 0==fsPS[1]){
    cout << "ERROR in SF_DtoAP0_AtoTP1_TtoP2P3_BASE::parseTree"
	 << " Didn't find T or P1 " << T.get() << ", " << fsPS[1].get() << endl;
    return false;
  }
  
  if(T->nDgtr() != 2){
    cout << "ERROR in SF_DtoAP0_AtoTP1_TtoP2P3_BASE::parseTree"
	 << " A should have 2 daughters, but it says it has "
	 << T->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[2] = T->getDgtrTreePtr(0);
  fsPS[3] = T->getDgtrTreePtr(1);
  // normalOrder(fsPS[2], fsPS[3]);
  
  // this->printYourself();
  return true;
}
// ---
double SF_DtoAP0_AtoTP1_TtoP2P3::getVal(IDalitzEvent& evt){
  bool dbThis=false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) parseTree(evt.eventPattern());

  TLorentzVector pT = p(2, evt) + p(3, evt);
  TLorentzVector qT = p(2, evt) - p(3, evt);
  TLorentzVector pA = pT + p(1, evt);
  TLorentzVector qA = pT - p(1, evt);
  TLorentzVector pD = pA + p(0, evt);
  TLorentzVector qD = pA - p(0, evt);

  double MA = mRes(A, evt);
  double MT = mRes(T, evt);
  //double MV = mRes(S);
  
  ZTspin1 tA(qA, pA, MA);
  ZTspin2 tT(qT, pT, MT);

  const double units = GeV*GeV*GeV*GeV;
    
  if (_useZemachTensors) {
      ZTspin1 LD(qD,pD,pD.M());
      SpinSumT PT(pT,MT);
      //return (tT.Contract(tA)).Dot(LD);
      return PT.Sandwich(qT,qT,LD,tA)/units;
  }

  double returnVal =  (tT.Contract(tA)).Dot(qD) / units;

  if(dbThis){
    cout << " SF_DtoAP0_AtoTP1_TtoP2P3::getVal "
	 << " returning " << returnVal
	 << endl;
  }
  return returnVal;
}

void SF_DtoAP0_AtoTP1_TtoP2P3::printYourself(ostream& os) const{
  //  bool debugThis = false;
  if(! ( fsPS[0] && fsPS[1] && fsPS[2] && fsPS[3]) ) return;
  os << "spin factor SF_DtoAP0_AtoTP1_TtoP2P3"
     << "\n\t  ZTspin1 tA(qA, pA, MA);"
     << "\n\t  ZTspin2 tT(qT, pT, MT);"
     << "\n\t  pT = p(2) + p(3); qT = p(2) - p(3);"
     << "\n\t  pA = pT + p(1); qA = pT - p(1);"
     << "\n\t  pD = pA + p(0); qD = pA - p(0);"
     << "\n\t  (tT.Contract(tA)).Dot(qD) / GeV^4"
     << "\n\t    parsed tree " << theDecay().oneLiner()
     << "\n      like this:" << endl;
  this->printParsing(os);
}

//====================================

//====================================
