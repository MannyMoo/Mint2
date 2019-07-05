#include "Mint/binflipChi2.h"
#include "typeinfo"

using namespace MINT;
using namespace std;

binflipChi2::binflipChi2(int nbinsPhase, int nbinsTime, float* Xreal, float* Ximag, float* r, float* tAv, float* tSqAv,
			 TH2F pHistD0, TH2F pHistD0bar, TH2F nHistD0, TH2F nHistD0bar, float ReZcp, float ImZcp, float ReDz, float ImDz):
  m_nbinsPhase(nbinsPhase),
  m_nbinsTime(nbinsTime),
  m_r(r),
  m_tAv(tAv),
  m_tSqAv(tSqAv),
  m_pHistD0(pHistD0),
  m_pHistD0bar(pHistD0bar),
  m_nHistD0(nHistD0),
  m_nHistD0bar(nHistD0bar),
  m_ReZcp("ReZcp", FitParameter::FIT, ReZcp, 0.0001, 0, 0, getParSet()),
  m_ImZcp("ImZcp", FitParameter::FIT, ImZcp, 0.0001, 0, 0, getParSet()),
  m_ReDz("ReDz", FitParameter::FIT, ReDz, 0.0001, 0, 0, getParSet()),
  m_ImDz("ImDz", FitParameter::FIT, ImDz, 0.0001, 0, 0, getParSet())
{
    for(int i = 0; i < m_nbinsPhase; i++){
        complex<float> Xb(Xreal[i], Ximag[i]);
        m_X.push_back(Xb);
    }
}

double binflipChi2::getVal(){    

    double chi2 = 0.;
    vector<vector<TGraph> > fits(m_nbinsPhase);
    complex<float> zcp(m_ReZcp.getCurrentFitVal(), m_ImZcp.getCurrentFitVal());
    complex<float> deltaz(m_ReDz.getCurrentFitVal(), m_ImDz.getCurrentFitVal());

    fits = getFits(zcp, deltaz);
    
    for(int b = 1; b <= m_nbinsPhase; b++){
    
        double* Rvals_pl = fits[0][b-1].GetY();
        double* Rvals_mi = fits[0][b-1].GetY();  
    
        for(int j = 0; j <= m_nbinsTime; j++){
	    double R_pl = Rvals_pl[j-1];
            double R_mi = Rvals_mi[j-1]; 
            float D0_term = 0, D0bar_term = 0;
            float D0_num = 0, D0_den = 0, D0bar_num = 0, D0bar_den = 0; 

            D0_num = pow(m_nHistD0.GetBinContent(j,b) - m_pHistD0.GetBinContent(j,b) * R_pl, 2);
            D0_den =   pow(m_nHistD0.GetBinError(j,b), 2) + pow(m_pHistD0.GetBinError(j,b) * R_pl, 2); 
            if(D0_den != 0){
	        D0_term = D0_num/D0_den;
	    }

            D0bar_num = pow(m_nHistD0bar.GetBinContent(j,b) - m_pHistD0bar.GetBinContent(j,b) * R_mi, 2);
            D0bar_den =   pow(m_nHistD0bar.GetBinError(j,b), 2) + pow(m_pHistD0bar.GetBinError(j,b) * R_mi, 2); 
            if(D0bar_den != 0){
	        D0bar_term = D0bar_num/D0bar_den;
	    }

            chi2 += D0_term + D0bar_term;
            
        }
    }
    
    return chi2;
}


vector<vector<TGraph> > binflipChi2::getFits(complex<float >zcp, complex<float> deltaz){

  vector<vector<TGraph> > fits(m_nbinsPhase);
    
    for(int i : {0,1}){

        for(int b = 1; b <= m_nbinsPhase; b++){
	    TGraph temp = TGraph(m_nbinsTime);
	    fits[i].push_back( temp );

            for(int j = 1; j <= m_nbinsTime; j++){
                float numerator = 0., denominator = 0.;
     
                numerator = m_r[b-1] * ( 1 + 0.25 * m_tSqAv[j-1] * ( pow(zcp,2) - pow(deltaz,2) ).real() );
                numerator += 0.25 * m_tSqAv[j-1] * abs(zcp + (float)pow(-1, i) * pow(deltaz, 2));
    	        numerator += sqrt(m_r[b-1]) * m_tAv[j-1] * ( conj(m_X[b-1]) * (zcp + (float)pow(-1, i) * deltaz) ).real();

                denominator = 1 + 0.25 * m_tSqAv[j-1] * ( pow(zcp,2) - pow(deltaz,2) ).real();
                denominator += m_r[b-1] * 0.25 * m_tSqAv[j-1] * abs(zcp + (float)pow(-1, i) * pow(deltaz, 2));
    	        denominator += sqrt(m_r[b-1]) * m_tAv[j-1] * ( m_X[b-1] * (zcp + (float)pow(-1, i) * deltaz) ).real();

                float Rval = numerator / denominator;
    	        fits[i][b-1].SetPoint(j-1, m_tAv[j-1], Rval);
            }
        }
    }
    
    return fits;
}

