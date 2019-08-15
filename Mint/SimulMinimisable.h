#ifndef __SIMULMINIMISABLE_H__
#define __SIMULMINIMISABLE_H__

#include <Mint/Minimisable.h>
#include <deque>

namespace MINT {
  class SimulMinimisable : public MINT::Minimisable {
    typedef std::deque<MINT::Minimisable*> Minimisables ;

    SimulMinimisable(MINT::Minimisable*, MINT::Minimisable*) ;
    SimulMinimisable(const Minimisables&) ;
    virtual ~SimulMinimisable() ;
    virtual double getVal() override ;
    virtual void parametersChanged() override ;

  protected :
    Minimisables m_minimisables ;
    void init() ;
  } ;
}

#endif
