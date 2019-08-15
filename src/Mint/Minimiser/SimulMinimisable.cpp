#include <Mint/SimulMinimisable.h>

using MINT::Minimisable ;

namespace MINT {
  SimulMinimisable::SimulMinimisable(Minimisable* mini1, Minimisable* mini2) :
    Minimisable(new MinuitParameterSet()),
    m_minimisables(1, mini1)
  {
    m_minimisables.push_back(mini2) ;
    init() ;
  }
  
  SimulMinimisable::SimulMinimisable(const SimulMinimisable::Minimisables& minimisers) :
    Minimisable(new MinuitParameterSet()),
    m_minimisables(minimisers)
  {
    init() ;
  }

  SimulMinimisable::~SimulMinimisable() {
    delete getParSet() ;
  }

  void SimulMinimisable::init() {
    for(auto* mini : m_minimisables){
      MinuitParameterSet* parset = mini->getParSet() ;
      for(unsigned ipar = 0 ; ipar < parset->size() ; ++ipar)
	getParSet()->add(parset->getParPtr(ipar)) ;
    }
  }

  double SimulMinimisable::getVal() {
    double val = 0. ;
    for(auto* mini : m_minimisables)
      val += mini->getVal() ;
    return val ;
  }

  void SimulMinimisable::parametersChanged() {
    Minimisable::parametersChanged() ;
    for(auto* mini : m_minimisables)
      mini->parametersChanged() ;
  }

}
