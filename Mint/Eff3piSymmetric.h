#ifndef MINT_EFF_THREE_PI_SYMMETRIC_HH
#define MINT_EFF_THREE_PI_SYMMETRIC_HH

#include "Mint/IReturnRealForEvent.h"
#include "Mint/IDalitzEvent.h"
#include "math.h"

using namespace MINT;
using namespace std;

class Eff3piSymmetric : public IReturnRealForEvent<IDalitzEvent> {
    double m_p0;
    double m_p1;
    double m_p2;
    double m_p3;

  public:
    Eff3piSymmetric(double p0 = 0.599125, double p1 = 5.25471 * pow(10,-8), 
		                double p2 = -5.83156 * pow(10,-14), double p3 = -1.80603 * pow(10,-21));
    virtual double RealVal(IDalitzEvent& evt);

};

#endif


