#include "Mint/IMinimisable.h"
#include "Mint/MinuitParameterSet.h"
#include "Mint/Minimisable.h"
#include "Mint/FitParameter.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "cmath"
#include "complex"

using namespace MINT;
using namespace std;


class binflipChi2 : public Minimisable{
    int m_nbinsPhase;
    int m_nbinsTime;
    vector<complex<double> > m_X;
    vector<double> m_r;
    vector<double> m_tAv;
    vector<double> m_tSqAv;
    TH2F m_pHistD0;
    TH2F m_pHistD0bar;
    TH2F m_nHistD0;
    TH2F m_nHistD0bar;
    FitParameter m_ReZcp;
    FitParameter m_ImZcp;
    FitParameter m_ReDz;
    FitParameter m_ImDz;
    int m_fakeData;
    vector<double> m_Fm;
    vector<double> m_Fp;
    int m_verbosity;

  public:
    binflipChi2(vector<complex<double> > X, vector<double> r, vector<double> tAv, vector<double> tSqAv, 
                       TH2F pHistD0, TH2F pHistD0bar, TH2F nHistD0, TH2F nHistD0bar, double ReZcp, double ImZcp, 
                       double ReDz, double ImDz, double stepSize, int fakeData = 0, vector<double> Fm = vector<double>(), 
		       vector<double> Fp = vector<double>(), int verbosity = 0);
    ~binflipChi2();
    double getVal();
    vector<vector<TGraph> > getFits();
    void genFakeData();
};

