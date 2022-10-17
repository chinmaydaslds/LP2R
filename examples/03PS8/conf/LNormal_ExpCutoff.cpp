/* Log-normal MWD with an exponential cut-off at the higher mass end.
   Parameters fixed in the main() */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>


double aaerfcc(const double x)
{
        double t,z,ans;

        z=fabs(x);
        t=1.0/(1.0+0.5*z);
        ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                t*(-0.82215223+t*0.17087277)))))))));
        return x >= 0.0 ? ans : 2.0-ans;
}

double LogNormalWt(const double Mw, const double PDI, const double M1, const double M2, double &MwBin)
{
double mu=log(Mw) - 1.50*log(PDI);
double sigma=sqrt(log(PDI));
double WtBin;
if(M1 < 0.0){ // (0, M2)
WtBin=0.5*aaerfcc( (mu + sigma*sigma - log(M2))/(sqrt(2.0)*sigma) );
double t2=aaerfcc( (mu + 2.0*sigma*sigma - log(M2))/(sqrt(2.0)*sigma) );
MwBin=0.50*Mw*t2/WtBin;  // t1=erfc(infinity)=0
            }
else{
 if(M2 < 0.0){ // (M1, infinity)
  WtBin=1.0 - 0.5*aaerfcc( (mu + sigma*sigma - log(M1))/(sqrt(2.0)*sigma) );
  double t1=aaerfcc( (mu + 2.0*sigma*sigma - log(M1))/(sqrt(2.0)*sigma) );
  MwBin=0.50*Mw*(2.0 - t1)/WtBin;
             }
 else{ // (M1, M2)
    double w1=0.5*aaerfcc( (mu + sigma*sigma - log(M1))/(sqrt(2.0)*sigma) );
    double w2=0.5*aaerfcc( (mu + sigma*sigma - log(M2))/(sqrt(2.0)*sigma) );
    WtBin=w2-w1;
    double t1=aaerfcc( (mu + 2.0*sigma*sigma - log(M1))/(sqrt(2.0)*sigma) );
    double t2=aaerfcc( (mu + 2.0*sigma*sigma - log(M2))/(sqrt(2.0)*sigma) );
    MwBin=0.50*Mw*(t2 - t1)/WtBin;
     }
    }
return WtBin;
}

int main()
{
double mw=390000.0;
double pdi=2.90;
double m_exp_cut=3.5e6;
int n=100; // 100 mw bins
std::vector<double> M_ar;
std::vector<double> wt_ar;
double mu=log(mw)-1.50*log(pdi);
double sigma=sqrt(log(pdi));
double lnMlow=mu + sigma*sigma-2.63*sqrt(2.0)*sigma;
double lnMhigh=mu + sigma*sigma+2.63*sqrt(2.0)*sigma;
double wtBin, MwBin, M1, M2;
M2=exp(lnMlow);
wtBin=LogNormalWt(mw, pdi, -1, M2, MwBin);
M_ar.push_back(MwBin);
wt_ar.push_back(wtBin);
double delta_lnM=(lnMhigh-lnMlow)/((double) (n-2));
for(int i=0; i< (n-2); i++){
 M1=exp(lnMlow + ((double) (i))*delta_lnM);
 M2=exp(lnMlow + ((double) (i+1))*delta_lnM);
wtBin=LogNormalWt(mw, pdi, M1, M2, MwBin);
M_ar.push_back(MwBin);
wt_ar.push_back(wtBin);
                           }
M1=exp(lnMhigh);
wtBin=LogNormalWt(mw, pdi, M1, -1, MwBin);
M_ar.push_back(MwBin);
wt_ar.push_back(wtBin);

double wttot=0.0;
for(int i=0; i<100; i++){
wt_ar[i] *= exp(-M_ar[i]/m_exp_cut);
wttot+=wt_ar[i];
                        }

std::fstream f1;
f1.open("PS390Kmcut.dat", std::ios::out);
for(int i=0; i<100; i++){
f1<<M_ar[i] << "  " << wt_ar[i]/wttot << std::endl;
                        }
f1.close();
return 0;
}
