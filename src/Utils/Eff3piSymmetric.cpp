#include "Mint/Eff3piSymmetric.h"

using namespace MINT;
using namespace std;

Eff3piSymmetric::Eff3piSymmetric(double p0, double p1, double p2, double p3):
  m_p0(p0),
  m_p1(p1),
  m_p2(p2),
  m_p3(p3)
{
}

double Eff3piSymmetric::RealVal(IDalitzEvent& evt){
    double s13 = evt.s(1,3);
    double s23 = evt.s(2,3);
    return m_p0 + m_p1*(s13+s23) + m_p2*(pow(s13,2) + pow(s23,2)) + m_p3*(pow(s13,3) + pow(s23,3));
}



