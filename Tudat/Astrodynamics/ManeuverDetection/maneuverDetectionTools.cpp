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
    } // namespace maneuver_detection
} // namespace tudat
