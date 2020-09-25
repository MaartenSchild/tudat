/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *     Developed by Maarten Schild, 2020
 *
 *
*/
#ifndef MANEUVERDETECTIONTOOLS_H
#define MANEUVERDETECTIONTOOLS_H

#include <iostream>
#include "Tudat/Astrodynamics/Ephemerides/tleEphemeris.h"

namespace tudat
{
    namespace maneuver_detection
    {

    /*! Function to log message to console
     *
     * \param message string to log
     * \return void
     */
    double  tleToSemiMajorAxis(tudat::ephemerides::Tle tle, double gravitationalParameter);
    std::map< double, double > tleToSemiMajorAxis(std::vector<tudat::ephemerides::Tle> tleSeries, double gravitationalParameter);



    double calcMedian(double arr[], int n);
    std::tuple<double, double, double> IQR(double a[], int n);

    std::tuple<double, double> theilSen(std::vector< double> x, std::vector< double> y);

    std::map< double, double > theilSenCorrection(std::map<double, double> semiMajorAxisSeries, int windowSize);
    std::map< double, Eigen::VectorXd > determineThreshold(std::map<double, double> correctedSeries, int halfWindowSize, double min = 0, double max = 0);

    std::vector<double> slice(const std::vector<double> v, int start, int end);
    std::map< double, double > detectManeuver(std::map<double, double> correctedSeries, std::map< double, Eigen::VectorXd > thresholdMap);

    void fasper(std::vector<double> &x, std::vector<double> &y, const double ofac, const double hifac,
                std::vector<double> &px, std::vector<double> &py, int &nout, int &jmax, double &prob);

    void spread(const double y, std::vector<double> &yy, const double x, const int m);
    void avevar(std::vector<double> &data, double &ave, double &var);
    void realft(std::vector<double> &data, const int isign);
    void four1(double *data, const int n, const int isign);
    void four1(std::vector<double> &data, const int isign);
    double SQR(double x);
    double SIGN(const double &a, const double &b);
    void SWAP(double &a, double &b);
    double harmonicAnalysis(std::vector<double> x, std::vector<double> y);


    } // namespace maneuver_detection
} // namespace tudat
#endif // MANEUVERDETECTIONTOOLS_H
