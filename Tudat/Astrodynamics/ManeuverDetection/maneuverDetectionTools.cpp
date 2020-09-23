/*
 *      Developed by Maarten Schild, 2020
 */
#include <iostream>
#include "Tudat/Astrodynamics/ManeuverDetection/maneuverDetectionTools.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


namespace tudat
{
    namespace maneuver_detection
    {

        /*! Function to calculate semiMajorAxis of TLE
         *
         * \param tle Tle to calculate semi major axis from
         * \param gravitationalParameter Gravitational parameter of body
         * \return semiMajorAxis in km
         */
        double tleToSemiMajorAxis(tudat::ephemerides::Tle tle, double gravitationalParameter){
            double semiMajorAxis = std::pow(gravitationalParameter/std::pow(tle.getMeanMotion()/60, 2), 1.0/3.0);
            return semiMajorAxis;
        }
        std::map< double, double > tleToSemiMajorAxis(std::vector<tudat::ephemerides::Tle> tleSeries, double gravitationalParameter){
            std::map< double, double > semiMajorAxisSeries;
            for(auto tle : tleSeries){
                semiMajorAxisSeries[tle.getEpoch()] = tudat::maneuver_detection::tleToSemiMajorAxis(tle, gravitationalParameter);
            }
            return semiMajorAxisSeries;
        }

        double calcMedian(double arr[], int n){
            // First we sort the array
            std::sort(arr, arr+n);
            // check for even case
            if (n % 2 != 0)
                return (double)arr[n/2];
            return (double)(arr[(n-1)/2] + arr[n/2])/2.0;
        }

        int indexMedian(int l, int r)
        {
            int n = r - l + 1;
            n = (n + 1) / 2 - 1;
            return n + l;
        }

        // Function to calculate IQR
        std::tuple<double, double, double> IQR(double a[], int n)
        {            
            std::sort(a,a+n);            
            // Index of median of entire data
            int mid_index = indexMedian(0, n);
            // Median of first half
            double Q1 = a[indexMedian(0, mid_index)];
            // Median
            double Q2 = a[mid_index];
            // Median of second half
            double Q3 = a[indexMedian(mid_index + 1, n)];
            // IQR calculation
            //return (Q3 - Q1);
            return std::make_tuple(Q1, Q2, Q3);
        }

        std::tuple<double, double> theilSen(std::vector< double> x, std::vector< double> y) {
            //TODO: make private
            int len = x.size();
            double slopes[ len * len ];
            int count = 0;
            for (int i = 0; i < len; i++) {
                for (int j = i + 1; j < len; j++) {
                    if (x[i] != x[j]) {
                        slopes[count++] = (y[j] - y[i]) / (x[j] - x[i]);
                    }
                }
            }

            double medianSlope = calcMedian(slopes, count);

            double cuts[len];
            for (int k = 0; k < len; k++) {
                cuts[k] = y[k] - medianSlope * (double)x[k];
            }

            double slope = medianSlope;
            double intercept = calcMedian(cuts, len);
            return std::make_tuple(intercept, slope);
        }

        std::vector<double> slice(const std::vector<double> v, int start=0, int end=-1) {
            // TODO: make private
            int oldlen = v.size();
            int newlen;
            if (end == -1 or end >= oldlen){
                newlen = oldlen-start;
            } else {
                newlen = end-start;
            }

            std::vector<double> nv(newlen);

            for (int i=0; i<newlen; i++) {
                nv[i] = v[start+i];
            }
            return nv;
        }
        std::map< double, double > theilSenCorrection(std::map<double, double> semiMajorAxisSeries, int halfWindowSize){
            std::vector<double> time;
            std::vector<double> a;
            std::map< double, double > correctedSeries;
            for(auto &imap: semiMajorAxisSeries)
            {
                time.push_back(imap.first);
                a.push_back(imap.second);
            }
            double intercept;
            double slope;
            //for (int point = halfWindowSize; point< (semiMajorAxisSeries.size() - halfWindowSize); point++ )
            for (int point = halfWindowSize; point< (semiMajorAxisSeries.size() ); point++ )
            {
                //std::vector<double> timeWindow = slice(time, point - halfWindowSize, point + halfWindowSize + 1);
                std::vector<double> timeWindow = slice(time, point - halfWindowSize, point);
                std::vector<double> aWindow = slice(a, point - halfWindowSize, point + halfWindowSize + 1);
                std::tie(intercept, slope) = theilSen(timeWindow, aWindow);
                correctedSeries[time[point]] = a[point] - intercept - slope * time[point];
            }
            return correctedSeries;
        }
        std::map< double, Eigen::VectorXd > determineThreshold(std::map<double, double> correctedSeries, int halfWindowSize, double min , double max){
            // Lines below are for possible unit teste
            //int n = theilSenCorrected.size();
            //double testarr[] = { 1, 19, 7, 6, 5, 9, 12, 27, 18, 2, 15 };
            //double Q1, Q2, Q3;
            //std::tie(Q1, Q2, Q3) = tudat::maneuver_detection::IQR(correctedSemiMajorAxisArray, n);
            //double IQR = (Q3 - Q1);
            //double IQR = tudat::maneuver_detection::IQR(testarr, 11);
            //std::cout << Q1 << " " << Q2 << " " << Q3<< " IQR "<<IQR << std::endl;
            std::map< double, Eigen::VectorXd > thresholdMap;
            int lengthSeries = correctedSeries.size();
            double correctedArray[correctedSeries.size()];
            std::vector<double> correctedVector;
            std::vector<double> correctedTime;
            int i =0;
            for(const auto &map : correctedSeries){
                correctedArray[i++] = map.second;
                correctedVector.push_back(map.second);
                correctedTime.push_back(map.first);
            }

            double Q1, Q2, Q3;
            std::tie(Q1, Q2, Q3) = IQR(correctedArray, lengthSeries);
            //std::cout << " Q1: " << Q1 << " Q2: " << Q2 << " Q3: " << Q3 << " IQR " << Q3-Q1<< std::endl;
            double minThresholdUpper = std::max(Q3 + (Q3-Q1), max);
            double minThresholdLower = std::min(Q1 - (Q3-Q1),min);
            Eigen::VectorXd threshold(2);
            for(auto &imap: correctedSeries)
            {
                threshold[0] = minThresholdUpper;
                threshold[1] = minThresholdLower;
                thresholdMap[imap.first] = threshold;
            }
            for (int point = halfWindowSize; point< (correctedSeries.size() - halfWindowSize); point++ )
            //for (int point = backWindowSize; point< (correctedSeries.size() ); point++ )
            {
                std::vector<double> correctedWindow = tudat::maneuver_detection::slice(correctedVector, point - halfWindowSize, point + 1 +halfWindowSize);

                std::tie(Q1, Q2, Q3) = tudat::maneuver_detection::IQR(correctedWindow.data(), correctedWindow.size());
                double thresholdUpper = Q3 + (Q3-Q1);
                double thresholdLower = Q1 - (Q3-Q1);
                Eigen::VectorXd threshold(2);
                threshold[0] = std::max(thresholdUpper, minThresholdUpper);
                threshold[1] = std::min(thresholdLower, minThresholdLower);
                //std::cout << correctedTime[point] << std::endl;
                thresholdMap[correctedTime[point]] = threshold;
                //std::cout << "minThresholdLower: " << minThresholdLower << " thresholdLower: " << thresholdLower << " minThresholdUpper " << minThresholdUpper << " thresholdUpper " << thresholdUpper<< std::endl;

            }
            return thresholdMap;
        }
        Eigen::VectorXd detectManeuver(std::map<double, double> correctedSeries, std::map< double, Eigen::VectorXd > thresholdMap){
            std::vector<double> maneuvers;
            int n = 0;
            double t = 0;

            for(auto &imap: correctedSeries)
            {
                if (imap.second > thresholdMap[imap.first][0] || imap.second < thresholdMap[imap.first][1]){
                    //maneuvers.push_back(imap.first);
                    //std::cout << imap.first << " " <<imap.second << " " << thresholdMap[imap.first][0] << std::endl;
                    n++;

                } else if (n > 2){
                    maneuvers.push_back(imap.first);
                    //std::cout << imap.first << " " <<imap.second << " " << thresholdMap[imap.first][0] << std::endl;
                    n = 0;
                } else {
                    n = 0;

                }
            }
            Eigen::VectorXd maneuversEigen(maneuvers.size());
            for(int i = 0; i<maneuvers.size(); i++)
            {
                maneuversEigen[i] = maneuvers[i];
            }
            return maneuversEigen;
        }
        void fasper(std::vector<double> &x, std::vector<double> &y, double ofac, double hifac, std::vector<double> &px, std::vector<double> &py, int &nout, int &jmax, double &prob) {
            const int MACC = 4;
            int j,k,nwk,nfreq,nfreqt,n=x.size(),np=px.size();
            double ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt,hs2wt,
                    hypo,pmax,sterm,swt,var,xdif,xmax,xmin;
            nout=(int)(0.5*ofac*hifac*n);
            nfreqt=(int)(ofac*hifac*n*MACC);
            nfreq=64;
            while (nfreq < nfreqt) nfreq <<= 1;
            nwk=nfreq << 1;
            if (np < nout) {px.resize(nout); py.resize(nout);}
            avevar(y,ave,var);

            if (var == 0.0) throw("zero variance in fasper");
            xmin=x[0];
            xmax=xmin;
            for (j=1;j<n;j++) {
                if (x[j] < xmin) xmin=x[j];
                if (x[j] > xmax) xmax=x[j];
            }
            xdif=xmax-xmin;
            std::vector<double> wk1(nwk,0.);
            std::vector<double> wk2(nwk,0.);
            fac=nwk/(xdif*ofac);
            fndim=nwk;
            for (j=0;j<n;j++) {
                ck=fmod((x[j]-xmin)*fac,fndim);
                ckk=2.0*(ck++);
                ckk=fmod(ckk,fndim);
                ++ckk;
                spread(y[j]-ave,wk1,ck,MACC);
                spread(1.0,wk2,ckk,MACC);
            }
            realft(wk1,1);
            realft(wk2,1);
            df=1.0/(xdif*ofac);
            pmax = -1.0;
            for (k=2,j=0;j<nout;j++,k+=2) {
                hypo=sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]);
                hc2wt=0.5*wk2[k]/hypo;
                hs2wt=0.5*wk2[k+1]/hypo;
                cwt=sqrt(0.5+hc2wt);
                swt=SIGN(sqrt(0.5-hc2wt),hs2wt);
                den=0.5*n+hc2wt*wk2[k]+hs2wt*wk2[k+1];
                cterm=SQR(cwt*wk1[k]+swt*wk1[k+1])/den;
                sterm=SQR(cwt*wk1[k+1]-swt*wk1[k])/(n-den);
                px[j]=(j+1)*df;
                py[j]=(cterm+sterm)/(2.0*var);
                if (py[j] > pmax) pmax=py[jmax=j];
                }
                expy=exp(-pmax);
                effm=2.0*nout/ofac;
                prob=effm*expy;
                if (prob > 0.01) prob=1.0-pow(1.0-expy,effm);
            }

        void spread(const double y, std::vector<double> &yy, const double x, const int m) {
                    static int nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880};
                    int ihi,ilo,ix,j,nden,n=yy.size();
                    double fac;
                    if (m > 10) throw("factorial table too small in spread");
                    ix=(int)x;
                    if (x == (double)ix) yy[ix-1] += y;
                    else {
                        ilo=std::min(std::max((int)(x-0.5*m),0),(int)n-m);
                        ihi=ilo+m;
                        nden=nfac[m];
                        fac=x-ilo-1;
                        for (j=ilo+1;j<ihi;j++) fac *= (x-j-1);
                        yy[ihi-1] += y*fac/(nden*(x-ihi));
                        for (j=ihi-1;j>ilo;j--) {
                            nden=(nden/(j-ilo))*(j-ihi);
                            yy[j-1] += y*fac/(nden*(x-j));
                        }
                    }
                }
        void avevar(std::vector<double> &data, double &ave, double &var) {
                double s,ep;
                int j,n=data.size();
                ave=0.0;
                for (j=0;j<n;j++) ave += data[j];
                ave /= n;
                var=ep=0.0;
                for (j=0;j<n;j++) {
                    s=data[j]-ave;
                    ep += s;
                    var += s*s;
                }
                var=(var-ep*ep/n)/(n-1);
            }
        void realft(std::vector<double> &data, const int isign) {
            int i,i1,i2,i3,i4,n=data.size();
            double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp;
            double theta=3.141592653589793238/double(n>>1);
            if (isign == 1) {
            c2 = -0.5;
            four1(data,1);
            } else {
            c2=0.5; theta = -theta;
            }
            wtemp=sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=sin(theta);
            wr=1.0+wpr;
            wi=wpi;
            for (i=1;i<(n>>2);i++) {
            i2=1+(i1=i+i);
            i4=1+(i3=n-i1);
            h1r=c1*(data[i1]+data[i3]);
            h1i=c1*(data[i2]-data[i4]);
            h2r= -c2*(data[i2]+data[i4]);
            h2i=c2*(data[i1]-data[i3]);
            data[i1]=h1r+wr*h2r-wi*h2i;

            data[i2]=h1i+wr*h2i+wi*h2r;
            data[i3]=h1r-wr*h2r+wi*h2i;
            data[i4]= -h1i+wr*h2i+wi*h2r;
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
            }
            if (isign == 1) {
            data[0] = (h1r=data[0])+data[1];

            data[1] = h1r-data[1];
            } else {
            data[0]=c1*((h1r=data[0])+data[1]);
            data[1]=c1*(h1r-data[1]);
            four1(data,-1);
            }

        }
        void four1(double *data, const int n, const int isign) {
            int nn,mmax,m,j,istep,i;
            double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
            if (n<2 || n&(n-1)) throw("n must be power of 2 in four1");
            nn = n << 1;
            j = 1;
            for (i=1;i<nn;i+=2) {
            if (j > i) {
            SWAP(data[j-1],data[i-1]);
            SWAP(data[j],data[i]);
            }
            m=n;
            while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
            }
            j += m;
            }
            mmax=2;
            while (nn > mmax) {
                istep=mmax << 1;
                theta=isign*(6.28318530717959/mmax);
                wtemp=sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi=sin(theta);
                wr=1.0;
                wi=0.0;
                for (m=1;m<mmax;m+=2) {
                    for (i=m;i<=nn;i+=istep) {
                        j=i+mmax;
                        tempr=wr*data[j-1]-wi*data[j];
                        tempi=wr*data[j]+wi*data[j-1];
                        data[j-1]=data[i-1]-tempr;
                        data[j]=data[i]-tempi;
                        data[i-1] += tempr;
                        data[i] += tempi;
                    }
                    wr=(wtemp=wr)*wpr-wi*wpi+wr;
                    wi=wi*wpr+wtemp*wpi+wi;
                }
                mmax=istep;
            }
        }
        void four1(std::vector<double> &data, const int isign) {
            four1(&data[0],data.size()/2,isign);
        }

        double SQR(double x){
            return x*x;
            }
        double SIGN(const double &a, const double &b)
            {return (double)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

        void SWAP(double &a, double &b)
            {double dum=a; a=b; b=dum;}

        } // namespace maneuver_detection
    } // namespace tudat
