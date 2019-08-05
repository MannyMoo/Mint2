#include <Mint/HadronicParameters.h>
#include <TMath.h>
#include <iostream>

using namespace std ;

int HadronicParameters::EventPair::binSign() const {
  return isPlus() ? 1 : -1 ;
}

bool HadronicParameters::EventPair::isPlus() const {
  return evt == evtPlus ;
}

HadronicParameters::BinSign::BinSign(HadronicParameters::ModelPtr model) :
  m_model(model)
{}

HadronicParameters::EventPair HadronicParameters::BinSign::classifyEvent(const IDalitzEvent& evt) const {
  HadronicParameters::EventPair evtPair ;
  evtPair.evt = MINT::counted_ptr<IDalitzEvent>(new DalitzEvent(evt)) ;
  DalitzEvent* cpevt = new DalitzEvent(evt) ;
  cpevt->CP_conjugateYourself() ;
  evtPair.cpevt = MINT::counted_ptr<IDalitzEvent>(cpevt) ;
  double f = m_model->RealVal(*evtPair.evt) ;
  double fcp = m_model->RealVal(*evtPair.cpevt) ;
  if(f > fcp){
    evtPair.evtPlus = evtPair.evt ;
    evtPair.evtMinus = evtPair.cpevt ;
    evtPair.Fplus = f ;
    evtPair.Fminus = fcp ;
  }
  else{
    evtPair.evtPlus = evtPair.cpevt ;
    evtPair.evtMinus = evtPair.evt ;
    evtPair.Fplus = fcp ;
    evtPair.Fminus = f ;
  }
  return evtPair ;
}

HadronicParameters::Bin::Bin(double _Fplus, double _Fminus, const complex<double>& X,
			     double sumw, double sumw2, double norm) :
  m_Fplus(_Fplus),
  m_Fminus(_Fminus),
  m_X(X),
  m_sumw(sumw),
  m_sumw2(sumw2),
  m_norm(norm),
  m_model(0),
  m_cpmodel(0)
{
  m_X *= sqrt(Fplus() * Fminus()) ;
}

HadronicParameters::Bin::Bin(ModelPtr model, ModelPtr cpmodel) :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(1.),
  m_model(model),
  m_cpmodel(cpmodel)
{}

void HadronicParameters::Bin::add(const HadronicParameters::EventPair& evtpair,
				  double weight) {
  m_Fplus += evtpair.Fplus ;
  m_Fminus += evtpair.Fminus ;
  complex<double> Aplus = m_model->ComplexVal(*evtpair.evtPlus) ;
  complex<double> Aminus = m_cpmodel->ComplexVal(*evtpair.evtPlus) ;
  m_X += conj(Aplus) * Aminus ;
  m_sumw += weight ;
  m_sumw2 += weight*weight ;
}

double HadronicParameters::Bin::Fplus() const {
  return m_Fplus / m_sumw / m_norm ;
}

double HadronicParameters::Bin::Fminus() const {
  return m_Fminus / m_sumw / m_norm ;
}

complex<double> HadronicParameters::Bin::Xplus() const {
  return m_X / m_sumw / sqrt(Fplus() * Fminus()) / m_norm ;
}

complex<double> HadronicParameters::Bin::Xminus() const {
  return conj(Xplus()) ;
}

double HadronicParameters::Bin::norm() const {
  return m_norm ;
}

void HadronicParameters::Bin::setNorm(double norm) {
  m_norm = norm ;
}

void HadronicParameters::Bin::Print() const {
  complex<double> Xp = Xplus() ;
  complex<double> Xm = Xminus() ;
  cout << "Fplus: " << Fplus() << ", Fminus: " << Fminus() << " Xplus: "
       << Xp.real() << " + " << Xp.imag() << "i, Xminus: " 
       << Xm.real() << " + " << Xm.imag() << "i." << endl ;
}

HadronicParameters::HadronicParameters(const DalitzEventPattern& pat, const DalitzEventPattern& cppat,
				       const HadronicParameters::Bins& bins, 
				       MINT::counted_ptr<IBinSign> binSign) :
  m_pat(pat),
  m_cppat(cppat),
  m_bins(bins),
  m_model(0),
  m_cpmodel(0),
  m_binSign(binSign)
{}

HadronicParameters::HadronicParameters(const DalitzEventPattern& pat, const DalitzEventPattern& cppat,
				       HadronicParameters::ModelPtr model,
				       HadronicParameters::ModelPtr cpmodel,
				       unsigned nbins,
				       MINT::counted_ptr<IBinSign> binSign) :
  m_pat(pat),
  m_cppat(cppat),
  m_bins(nbins, HadronicParameters::Bin(model, cpmodel)),
  m_model(model),
  m_cpmodel(cpmodel),
  m_binSign(binSign)
{}

unsigned HadronicParameters::nBins() const {
  return m_bins.size() ;
}

int HadronicParameters::binNumber(const HadronicParameters::EventPair& evtPair) const {
  complex<double> Amp = m_model->ComplexVal(*evtPair.evtPlus) ;
  complex<double> Ampbar = m_cpmodel->ComplexVal(*evtPair.evtPlus) ;
  double phasediff = arg(Amp/Ampbar) ;
  if(phasediff < 0.)
    phasediff += 2. * TMath::Pi() ;
  return (int(phasediff * nBins()/2./TMath::Pi() + 0.5) % nBins() + 1) * evtPair.binSign() ;
}

int HadronicParameters::binNumber(const IDalitzEvent& evt) const {
  return binNumber(m_binSign->classifyEvent(evt)) ;
}

void HadronicParameters::add(const HadronicParameters::EventPair& evtPair, double weight) {
  int binno = abs(binNumber(evtPair)) - 1;
  m_bins[binno].add(evtPair, weight) ;
}

void HadronicParameters::add(const IDalitzEvent& evt, double weight) {
  return add(m_binSign->classifyEvent(evt), weight) ;
}

void HadronicParameters::add(TRandom3& rndm, unsigned nevt) {
  for(unsigned i = 0 ; i < nevt ; ++i)
    add(rndm) ;
}

void HadronicParameters::add(TRandom3& rndm) {
  DalitzEvent evt(m_pat, (TRandom*)&rndm) ;
  add(evt) ;
}

double HadronicParameters::normalise() {
  for(auto& bin : m_bins)
    bin.setNorm(1.) ;
  double norm = 0. ;
  for(const auto& bin : m_bins)
    norm += bin.Fplus() + bin.Fminus() ;
  for(auto& bin : m_bins)
    bin.setNorm(norm) ;
  return norm ;
}

void HadronicParameters::Print() const {
  unsigned i = 1 ;
  for(const auto& bin : m_bins){
    cout << "Bin " << i << ":" << endl ;
    bin.Print() ;
    i += 1 ;
  }
}
  
