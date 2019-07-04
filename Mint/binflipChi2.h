#include "Mint/IMinimisable.h"
#include "Mint/MinuitParameterSet.h"
#include "Mint/Minimisable.h"

using namespace MINT;
using namespace std;


class binflipChi2 : public Minimisable{
    
   
  public:
    binflipChi2();
   
    double getVal(){
        float chi2 = 0 ;
        return chi2 ;
    }
};

