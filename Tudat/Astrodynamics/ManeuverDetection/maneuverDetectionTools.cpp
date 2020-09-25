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
            std::vector<double> correctedTimeDay;
            std::vector<double> correctedTime;
            int i =0;
            for(const auto &map : correctedSeries){
                correctedArray[i++] = map.second;
                correctedVector.push_back(map.second);
                correctedTimeDay.push_back(tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(map.first));
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
                std::vector<double> timeWindow = tudat::maneuver_detection::slice(correctedTimeDay, point - halfWindowSize, point + 1 +halfWindowSize);

                double amplitude = tudat::maneuver_detection::harmonicAnalysis(timeWindow, correctedWindow);

                std::tie(Q1, Q2, Q3) = tudat::maneuver_detection::IQR(correctedWindow.data(), correctedWindow.size());
                double thresholdUpper = Q3 + (Q3-Q1) + amplitude;
                double thresholdLower = Q1 - (Q3-Q1) - amplitude;
                Eigen::VectorXd threshold(2);
                threshold[0] = std::max(thresholdUpper, minThresholdUpper);
                threshold[1] = std::min(thresholdLower, minThresholdLower);
                //std::cout << correctedTime[point] << std::endl;
                thresholdMap[correctedTime[point]] = threshold;
                //std::cout << "minThresholdLower: " << minThresholdLower << " thresholdLower: " << thresholdLower << " minThresholdUpper " << minThresholdUpper << " thresholdUpper " << thresholdUpper<< std::endl;

            }
            return thresholdMap;
        }
        std::map< double, double > detectManeuver(std::map<double, double> correctedSeries, std::map< double, Eigen::VectorXd > thresholdMap){
            // TODO more useful map info
            std::map< double, double > maneuverMap;
            int n = 0;
            double t = 0;
            int count = 0;
            int allowcount = 0;
            bool maneuvering = false;
            int allowed = 2; //Allowed values within threshold in maneuvers
            for(auto &imap: correctedSeries)
            {
                if (imap.second > thresholdMap[imap.first][0] || imap.second < thresholdMap[imap.first][1]){
                    if (maneuvering == false){
                        maneuvering = true;
                        t = imap.first;
                    }
                    n++;
                    allowcount = 0;
                }
                else {
                    if (allowcount<allowed){
                        allowcount++;
                    }
                    else{
                        if (n > 2){
                            count++;
                            maneuverMap[t] = count;
                        }
                        allowcount = 0;
                        n = 0;
                        maneuvering=false;
                    }
                }
            }

            return maneuverMap;
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
        double harmonicAnalysis(std::vector<double> x, std::vector<double> y){

            // HARMONIC ANALYSIS
            /*
            double xar[] = {37.454011884736246, 95.07143064099162, 73.1993941811405, 59.86584841970366, 15.601864044243651, 15.599452033620265, 5.8083612168199465, 86.61761457749351, 60.11150117432088, 70.80725777960456, 2.0584494295802447, 96.99098521619943, 83.24426408004217, 21.233911067827616, 18.182496720710063, 18.34045098534338, 30.42422429595377, 52.475643163223786, 43.194501864211574, 29.122914019804192, 61.18528947223795, 13.949386065204184, 29.214464853521815, 36.63618432936917, 45.606998421703594, 78.51759613930136, 19.967378215835975, 51.42344384136116, 59.24145688620425, 4.645041271999773, 60.75448519014384, 17.052412368729154, 6.505159298527952, 94.88855372533332, 96.56320330745594, 80.83973481164611, 30.46137691733707, 9.767211400638388, 68.42330265121569, 44.01524937396013, 12.203823484477883, 49.51769101112702, 3.4388521115218396, 90.9320402078782, 25.87799816000169, 66.2522284353982, 31.171107608941096, 52.00680211778108, 54.67102793432797, 18.485445552552704, 96.95846277645586, 77.51328233611146, 93.9498941564189, 89.48273504276489, 59.78999788110851, 92.18742350231169, 8.84925020519195, 19.59828624191452, 4.522728891053807, 32.53303307632643, 38.8677289689482, 27.134903177389592, 82.87375091519293, 35.67533266935893, 28.093450968738075, 54.26960831582485, 14.092422497476264, 80.21969807540397, 7.455064367977082, 98.68869366005173, 77.22447692966574, 19.87156815341724, 0.5522117123602399, 81.54614284548342, 70.68573438476172, 72.90071680409874, 77.12703466859458,
                            7.4044651734090365, 35.84657285442726, 11.586905952512971, 86.31034258755935, 62.329812682755794, 33.08980248526492, 6.355835028602364, 31.09823217156622, 32.518332202674706, 72.96061783380641, 63.75574713552131, 88.72127425763266, 47.22149251619493, 11.959424593830171, 71.3244787222995, 76.07850486168974, 56.127719756949624, 77.0967179954561, 49.379559636439076, 52.27328293819941, 42.75410183585496, 2.541912674409519, 10.789142699330444};
            double yar[] ={3.8201923923728454, 4.844782225244905, -2.8833670573663337, -3.0749463272387807, -4.720037151286264, -4.735276443806648, 2.415400063180868, -4.04375725158461, 4.23077948675483, 2.309624161558726, 4.551529640704266, -0.8126654182150839, -5.023779889016205, -4.720499571975827, -1.4601984654898927, 0.7656767237708686, 4.879144763221514, 2.183010915189709, -2.5438960072511136, 3.5254303820480963, -1.6897214514609042, -4.05271825403947, -3.9191504350845587, -2.7408639177566405, -4.65227104853875, -1.670221358875789, -2.9186807584441277, 4.87902965592326, -4.951438753491215, -1.9417638902008152, 5.1707600779058955, 4.191954480594582, -0.459731285767004, -4.321691090924678, -4.836234354804626, -0.6045134964453002, 3.333017135592395, 4.985491472574076, 4.941699009263856, 1.447664497031278, -3.2260489570498856, -1.7534761433074701, 4.68341829650406, -4.716334184033022, -3.65048703322235, -5.086528340250541, -0.27778087051711997, 0.4991443085461853, 0.4692626811590124, 1.573627393247171, -3.625950513296567, -1.295423261175356, -4.040975245877911, 1.5482657711861838, 3.4899134957105797, -1.8998979487563088, -1.583955436448821, -4.755345529086185, -2.1691562278880236, -2.7609918600892702, -3.097655488418317, 2.785753424541316, -3.3699668903611366, 0.6900376493225708, 4.932092701309594, -4.5316307781806895, 4.766048896766128, -4.187860617800865, 3.7725289937681747, 2.0950426393771533, -4.556130591377167, -3.431634768595316, -4.111863813699266, -3.7912114880928245, 1.7836967802211103, -4.7410793934279685, 3.32936782926948, 4.891565632183207, -1.2055764285637787, -5.0601019673768395, -1.9130596599211067, -0.2841867096934201, 4.8437418906218905, 2.1233727726085263, 4.706651507529654, -1.614880878269942, -3.2641973737380154, 4.88862146136544, 4.381103963587335, -4.254053822217937, -3.379859935177536, -0.6409788911213476, 4.954766360962028, 3.2743411779518934, 4.752759027856264, 3.744397575949538, -4.533865685643313, 5.019177541764259, -3.5242215090044366, 3.7823936689291484};

            std::vector<double> x (xar, xar + sizeof(xar) / sizeof(double) );
            std::vector<double> y (yar, yar + sizeof(yar) / sizeof(double) );
            std::vector<double> xs, ys;
            for(auto &imap: theilSenCorrected){
                xs.push_back(imap.first/(physical_constants::JULIAN_DAY));
                ys.push_back(imap.second);
            }
            */
            double ofac = 4;
            double hifac = 1;
            std::vector<double> px(100);
            std::vector<double> py(100);
            int nout = 100;
            int jmax = 0;
            double prob = 0;
            tudat::maneuver_detection::fasper(x, y, ofac, hifac, px, py, nout, jmax, prob);

            std::map< double, double > harmonicsTestMap;
            for (int l = 0; l<px.size(); l++){
                harmonicsTestMap[px[l]] = py[l];
            }
            double sig = 0;
            double n = y.size();
            for (int ind = 0; ind<n; ind++){
                sig+=y[ind]*y[ind];
            }
            sig = sig/(n-1);

            //std::cout <<"size: "<< px.size()<< " jmax: "<< jmax <<" px: " << px[jmax] <<" py: "<< py[jmax]<< " prob: "<< prob << " amplitude: "<< amp<< std::endl;
            //std::cout <<"min: "<< px[0]<< " max: "<< px[nout-1] << " size:"<< px.size()<< std::endl;


            //const std::string filePath( __FILE__ );
            /*
            const std::string folder = "C:/tudatBundle/tudatApplications/MScThesis/output/harmonics/";
            boost::gregorian::date date = tudat::basic_astrodynamics::convertJulianDayToCalendarDate(x[0]);
            std::string x0 = std::to_string(date.year())+"-"+std::to_string(date.month())+"-"+std::to_string(date.day());

            boost::gregorian::date datef = tudat::basic_astrodynamics::convertJulianDayToCalendarDate(x[x.size()-1]);
            std::string xf = std::to_string(datef.year())+"-"+std::to_string(datef.month())+"-"+std::to_string(datef.day());

            tudat::input_output::writeDataMapToTextFile(harmonicsTestMap, x0+"_"+xf+"harmonics.dat",
                                                        folder,
                                                        "",
                                                        std::numeric_limits< double >::digits10,
                                                        std::numeric_limits< double >::digits10,
                                                        "," );
            //  END HARMONICS
            */
            double amp;
            if (prob<0.001){
                amp = std::sqrt((2/n)*2*sig*py[jmax]);
            } else{
                amp = 0;
            }
            //std::cout <<"amp: "<< amp<< " date: "<< tudat::basic_astrodynamics::convertJulianDayToCalendarDate(x[0])<< " - "<< tudat::basic_astrodynamics::convertJulianDayToCalendarDate(x[x.size()-1])<<std::endl;
            return amp;
        }        } // namespace maneuver_detection
    } // namespace tudat
