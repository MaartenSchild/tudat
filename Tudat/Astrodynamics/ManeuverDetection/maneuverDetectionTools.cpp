/*
 *      Developed by Maarten Schild, 2020
 */
#include <iostream>
#include "Tudat/Astrodynamics/ManeuverDetection/maneuverDetectionTools.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

/*      Maneuver detection algorithm components
 *      Author: Maarten Schild
 *      Date: 21-03-2021
 *
 *      Used for MSc Thesis "Sun-synchronous Spacecraft compliance with International Space Debris Guidelines"
 *      http://resolver.tudelft.nl/uuid:f105b7fc-b9d6-484e-9c70-c76ba994d0a4
 *
 *
 *      Converts TLE series to semimajor axis series
 *      Theil-Sen-Siegel Slope correction
 *      Harmonic detection
 *      Treshold Generation
 *      Event Detection
 *
 *
 *
 */

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
            int n = x.size();
            //int len = n*(n-1)/2;
            int len = (n)*(n);
            double slopes[ len  ];
            double intercepts[ len  ];
            int count = 0;
            for (int i = 0; i < n; i++) {
                //for (int j = i + 1; j < n; j++) {
                for (int j = 0; j < n; j++) {
                    if (x[i] != x[j]) {
                        slopes[count++] = (y[j] - y[i]) / (x[j] - x[i]);
                        intercepts[count] = (x[j]*y[i] - x[i]*y[j]) / (x[j] - x[i]);
                        //std::cout << intercepts[count] << std::endl;
                    }
                }
            }

            double medianSlope = calcMedian(slopes, count);
            double cuts[n];
            for (int k = 0; k < n; k++) {
                cuts[k] = y[k] - medianSlope * (double)x[k];
            }
            double intercept = calcMedian(cuts, n);
            double slope = medianSlope;
            return std::make_tuple(intercept, slope);

        }

        std::vector<double> slice(const std::vector<double> v, int start=0, int end=-1) {

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
            for (int point = halfWindowSize; point< (semiMajorAxisSeries.size() - halfWindowSize); point++ )
            //for (int point = halfWindowSize; point< (semiMajorAxisSeries.size() ); point++ )
            {
                //std::vector<double> timeWindow = slice(time, point - halfWindowSize, point + halfWindowSize + 1);
                std::vector<double> timeWindow = slice(time, point - halfWindowSize, point);
                //std::vector<double> aWindow = slice(a, point - halfWindowSize, point + halfWindowSize + 1);
                std::vector<double> aWindow = slice(a, point - halfWindowSize, point);
                std::tie(intercept, slope) = theilSen(timeWindow, aWindow);
                correctedSeries[time[point]] = a[point] - intercept - slope * time[point];
            }
            return correctedSeries;
        }
        std::map< double, Eigen::VectorXd > determineThreshold(std::map<double, double> correctedSeries, int halfWindowSize, double min , double max){

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
            double amplitude, freq;
            int halfHarmonicWindowSize = 100;
            std::tie(amplitude, freq) = tudat::maneuver_detection::harmonicAnalysis(correctedTimeDay, correctedVector);
            //std::cout << "Amplitude: "<< amplitude << ">>" << freq <<std::endl;

            //double amplitudeArray[correctedSeries.size()-halfWindowSize*2] = {amplitude};
            std::vector<double> amplitudeVector(correctedSeries.size(), amplitude);
            for (int point = halfHarmonicWindowSize; point< (correctedSeries.size() - halfHarmonicWindowSize); point+=50 ){
                std::vector<double> correctedWindow = tudat::maneuver_detection::slice(correctedVector, point - halfHarmonicWindowSize, point + 1 +halfHarmonicWindowSize);
                std::vector<double> timeWindow = tudat::maneuver_detection::slice(correctedTimeDay, point - halfHarmonicWindowSize, point + 1 +halfHarmonicWindowSize);
                std::tie(amplitude, freq) = tudat::maneuver_detection::harmonicAnalysis(timeWindow, correctedWindow);
                //std::cout <<point << std::endl;
                for (int p = point-25; p < point+26; p++){
                    amplitudeVector[p] = amplitude;
                    //std::cout <<p << std::endl;
                }
            }


            for (int point = halfWindowSize; point< (correctedSeries.size() - halfWindowSize); point++ )
            //for (int point = halfWindowSize; point< (correctedSeries.size() ); point++ )
            {
                std::vector<double> correctedWindow = tudat::maneuver_detection::slice(correctedVector, point - halfWindowSize, point + 1 +halfWindowSize);
                std::vector<double> timeWindow = tudat::maneuver_detection::slice(correctedTimeDay, point - halfWindowSize, point + 1 +halfWindowSize);
                //std::vector<double> correctedWindow = tudat::maneuver_detection::slice(correctedVector, point - halfWindowSize, point );
                //std::vector<double> timeWindow = tudat::maneuver_detection::slice(correctedTimeDay, point - halfWindowSize, point );                

                double amplitude = amplitudeVector[point];
                std::tie(Q1, Q2, Q3) = tudat::maneuver_detection::IQR(correctedWindow.data(), correctedWindow.size());
                double thresholdUpper = Q3 + (Q3-Q1) + amplitude;
                double thresholdLower = Q1 - (Q3-Q1) - amplitude;
                Eigen::VectorXd threshold(2);
                threshold[0] = std::max(thresholdUpper, minThresholdUpper);
                threshold[1] = std::min(thresholdLower, minThresholdLower);
                //std::cout << correctedTime[point] << std::endl;
                thresholdMap[correctedTime[point]] = threshold;
                boost::gregorian::date date = tudat::basic_astrodynamics::convertJulianDayToCalendarDate(correctedTimeDay[point]);
                //std::cout << date.year() <<"-" << date.month()<< " Amplitude: "<< amplitude << "<<" <<Q1 - (Q3-Q1) << "<<" << Q3 + (Q3-Q1) << std::endl;
            }
            //std::cout << "test2" << std::endl;
            return thresholdMap;
        }
        std::map< double, double > detectManeuver(std::map<double, double> correctedSeries, std::map< double, Eigen::VectorXd > thresholdMap, int skipN){            
            std::map< double, double > maneuverMap;
            int n = 0;
            double t = 0;
            int count = 0;
            int allowcount = 0;
            bool maneuvering = false;
            int allowed = 2; //Allowed values within threshold in maneuvers
            int i = 0;
            for(auto &imap: correctedSeries)
            {
                if (i > skipN){
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
                i++;
            }
            return maneuverMap;
        }

        std::vector<double> linspace(double start, double end, int num_in)
        {
          std::vector<double> linspaced;
          double num = static_cast<double>(num_in);

          if (num == 0) { return linspaced; }
          if (num == 1)
            {
              linspaced.push_back(start);
              return linspaced;
            }

          double delta = (end - start) / (num - 1);

          for(int i=0; i < num-1; ++i)
            {
              linspaced.push_back(start + delta * i);
            }
          linspaced.push_back(end);
          return linspaced;
        }
        std::vector<double> lombScargle(std::vector<double> x, std::vector<double> y, std::vector<double> freqs)
        {
            // Follows the Scipy implementation of lombScargle
            std::vector<double> pgram;
            double c, s, xc, xs, cc, ss, cs;
            double tau, c_tau, s_tau, c_tau2, s_tau2, cs_tau;

            for(int i = 0; i < freqs.size(); i++){
                    xc = 0.;
                    xs = 0.;
                    cc = 0.;
                    ss = 0.;
                    cs = 0.;

                    for (int j = 0; j < x.size(); j++){
                        c = cos(freqs[i] * x[j]);
                        s = sin(freqs[i] * x[j]);

                        xc += y[j] * c;
                        xs += y[j] * s;
                        cc += c * c;
                        ss += s * s;
                        cs += c * s;
                    }

                    tau = atan2(2 * cs, cc - ss) / (2 * freqs[i]);
                    c_tau = cos(freqs[i] * tau);
                    s_tau = sin(freqs[i] * tau);
                    c_tau2 = c_tau * c_tau;
                    s_tau2 = s_tau * s_tau;
                    cs_tau = 2 * c_tau * s_tau;

                    pgram.push_back( 0.5 * ((pow(c_tau * xc + s_tau * xs, 2) /
                        (c_tau2 * cc + cs_tau * cs + s_tau2 * ss)) +
                        (pow(c_tau * xs - s_tau * xc, 2) /
                        (c_tau2 * ss - cs_tau * cs + s_tau2 * cc))));
            }
            return pgram;
        }
        std::tuple<double, double> harmonicAnalysis(std::vector<double> x, std::vector<double> y){
            double minFreq = 0.01;
            double maxFreq = 1.0;
            int nFreqs = 500;

            int n = x.size();

            std::vector<double> freqs = linspace(minFreq, maxFreq, nFreqs);
            std::vector<double> pgram = lombScargle(x, y, freqs);
            auto it = std::max_element(std::begin(pgram), std::end(pgram));
            double pmax = *it;
            double fmax = freqs[std::distance(pgram.begin(), it)];
            double amp = std::sqrt(4*pmax/n);

            return std::make_tuple(amp, fmax);
        }        } // namespace maneuver_detection
    } // namespace tudat
