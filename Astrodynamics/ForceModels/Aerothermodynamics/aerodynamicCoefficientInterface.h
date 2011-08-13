/*! \file aerodynamicCoefficientInterface.cpp
 *    This file contains the definition of the AerodynamicCoefficientInterface
 *    base class included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author           :  F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@tudelft.nl
 *
 *    Date created      : 08 June 2011
 *    Last modified     : 10 August 2011
 *
 *    References
 *
 *    Notes
 *
 *    The computeCoefficients() function is not yet implented in any derived classes and is
 *    therefor not pure virtual in this base class.
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110608    F.M. Engelen      First creation of code.
 *      110714    D. Dirkx          Class name change and other minor changes during code check.
 *      110715    F.M. Engelen      Added the virtual compute function.
 *      110810    J. Leloux         Corrected doxygen documentation (function variable name).
 */

#ifndef AERODYNAMICCOEFFICIENTINTERFACE_H
#define AERODYNAMICCOEFFICIENTINTERFACE_H

//includestatements
#include "linearAlgebra.h"

//! Base class to hold an Aerodynamic Coefficient Interface. This interface can, for instance,
//! be a database of coefficients or an aerodynamic analysis code which generates coefficients.
class AerodynamicCoefficientInterface
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     */
    AerodynamicCoefficientInterface( );

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~AerodynamicCoefficientInterface( );

    //! Sets reference area.
    /*!
     *  Sets reference area used to non-dimensionalize aerodynamic
     *  forces and moments.
     *  \param referenceArea Aerodynamic reference area.
     */
    void setReferenceArea( const double& referenceArea );

    //! Gets reference area.
    /*!
     *  Gets reference area used to non-dimensionalize aerodynamic
     *  forces and moments.
     *  \return Aerodynamic reference area.
     */
    double getReferenceArea( );

    //! Sets reference length.
    /*!
     *  Sets reference length used to non-dimensionalize aerodynamic
     *  moments.
     *  \param referenceLength Aerodynamic reference length.
     */
    void setReferenceLength( const double& referenceLength );

    //! Gets reference length.
    /*!
     *  Gets reference length used to non-dimensionalize aerodynamic
     *  moments.
     *  \return Aerodynamic reference length.
     */
    double getReferenceLength( );

    //! Sets lateral reference length.
    /*!
     *  Sets lateral reference length used to non-dimensionalize aerodynamic
     *  moments.
     *  \param lateralReferenceLength Aerodynamic reference length.
     */
    void setLateralReferenceLength( const double& lateralReferenceLength );

    //! Gets lateral reference length.
    /*!
     *  Gets lateral reference length used to non-dimensionalize aerodynamic
     *  moments.
     *  \return Aerodynamic  lateral reference length.
     */
    double getLateralReferenceLength( );

    //! Sets moment reference point.
    /*!
     *  Sets the point w.r.t. which the arm of the aerodynamic
     *  moment on a vehicle panel is determined.
     *  \param referencePoint Aerodynamic reference point.
     */
    void setMomentReferencePoint( const Vector3d& referencePoint );

    //! Gets moment reference point.
    /*!
     *  Gets the point w.r.t. which the arm of the aerodynamic
     *  moment on a vehicle panel is determined.
     *  \return Aerodynamic reference point.
     */
    VectorXd getMomentReferencePoint( );

    //! Retreives the current force coefficients.
    /*!
     *  Retreives the force coefficients that have been set as the 'current' force coefficients,
     *  i.e. at the current flight condition.
     * \return current force coefficients.
     */
    Vector3d getCurrentForceCoefficients( );

    //! Get the moment coefficients
    /*!
     *  Retreives the moment coefficients that have been set as the 'current' moment coefficients,
     *  i.e. at the current flight condition.
     * \return current moment coefficients.
     */
    Vector3d getCurrentMomentCoefficients( );

    //! Function to set the force Coefficients.
    /*!
     *  Sets the current force coefficients, i.e. at the current flight condition.
     * \param currentForceCoefficients the current force coefficients.
     */
    void setCurrentForceCoefficients( const Vector3d& currentForceCoefficients );

    //! Function to set the moment Coefficients.
    /*!
     *  Sets the current moment coefficients, i.e. at the current flight condition.
     * \param currentMomentCoefficients the current force coefficients.
     */
    void setCurrentMomentCoefficients( const Vector3d& currentMomentCoefficients );

    //! Function to compute the aerodynamic coefficients at current flight condition.
    /*!
     *  This function calculates the current force and moment coefficients and is to be
     *  implemented in derived classes. Such a calculation would be performed by, for instance,
     *  including a pointer to a Vehicle object, from which the current flight condition can
     *  be retrieved. The function here is not made pure virtual pending the inclusion of
     *  such functionality in deived classes.
     */
    virtual void computeCurrentCoefficients( ){}

protected:

    //! The current force coefficients
    /*!
     * The force coefficients at the current flight condition.
     */
    Vector3d currentForceCoefficients_;

    //! The current moment coefficients
    /*!
     * The moment coefficients at the current flight condition.
     */
    Vector3d currentMomentCoefficients_;

    //! Aerodynamic reference length.
    /*!
     *  Reference length with which aerodynamic moments are non-
     *  dimensionalized.
     */
    double referenceLength_ ;

    //! Aerodynamic reference area.
    /*!
     *  Reference area with which aerodynamic forces and moments are non-
     *  dimensionalized.
     */
    double referenceArea_ ;

    //! Lateral Aerodynamic reference length
    /*!
     *  Lateral reference length with which aerodynamic moments are non-
     *  dimensionalized.
     */
    double lateralReferenceLength_;

    //! Aerodynamic moment reference point.
    /*!
     *  Point w.r.t. which the arm of the moment on a vehicle panel is
     *  determined.
     */
    Vector3d momentReferencePoint_ ;

private:

};

#endif // AERODYNAMICCOEFFICIENTINTERFACE_H

// End of file.