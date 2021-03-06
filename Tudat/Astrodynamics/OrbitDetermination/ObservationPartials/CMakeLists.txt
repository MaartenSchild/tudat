 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(OBSERVATION_PARTIALS_SOURCES
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/angularPositionPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/rotationMatrixPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/positionPartials.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayRangePartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/twoWayDopplerPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/nWayRangePartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/differencedOneWayRangeRatePartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayDopplerPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/observationPartial.cpp"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/lightTimeCorrectionPartial.cpp"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/firstOrderRelativisticPartial.cpp"
)

# Set the header files.
set(OBSERVATION_PARTIALS_HEADERS
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/angularPositionPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/rotationMatrixPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/positionPartials.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayRangePartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/nWayRangePartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/differencedOneWayRangeRatePartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/oneWayDopplerPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/twoWayDopplerPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/observationPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/eulerAngleObservablePartials.h"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/lightTimeCorrectionPartial.h"
  "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/firstOrderRelativisticPartial.h"
)

# Set the header files.
set(OBSERVATION_PARTIALS_TEST_HEADERS
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/observationPartialTestFunctions.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/numericalObservationPartial.h"
)

set(OBSERVATION_PARTIALS_TEST_SOURCES
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/observationPartialTestFunctions.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/numericalObservationPartial.cpp"
)
# Add static libraries.
add_library(tudat_observation_partials STATIC ${OBSERVATION_PARTIALS_SOURCES} ${OBSERVATION_PARTIALS_HEADERS})
setup_tudat_library_target(tudat_observation_partials "${SRCROOT}{OBSERVATIONPARTIALSDIR}")

add_library(tudat_observation_partials_test_functions STATIC ${OBSERVATION_PARTIALS_TEST_SOURCES} ${OBSERVATION_PARTIALS_TEST_HEADERS})
setup_tudat_library_target(tudat_observation_partials "${SRCROOT}{OBSERVATIONPARTIALSDIR}/UnitTests")

add_executable(test_OneWayRangePartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestOneWayRangePartials.cpp")
setup_custom_test_program(test_OneWayRangePartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_OneWayRangePartials tudat_observation_partials_test_functions  ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_NWayRangePartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestNWayRangePartials.cpp")
setup_custom_test_program(test_NWayRangePartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_NWayRangePartials tudat_observation_partials_test_functions  ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_OneWayDifferencedRangePartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestDifferencedOneWayRangePartials.cpp")
setup_custom_test_program(test_OneWayDifferencedRangePartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_OneWayDifferencedRangePartials tudat_observation_partials_test_functions  ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_OneWayDopplerPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestOneWayDopplerPartials.cpp")
setup_custom_test_program(test_OneWayDopplerPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_OneWayDopplerPartials tudat_observation_partials_test_functions  ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_TwoWayDopplerPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestTwoWayDopplerPartials.cpp")
setup_custom_test_program(test_TwoWayDopplerPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_TwoWayDopplerPartials tudat_observation_partials_test_functions  ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_AngularPositionPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestAngularPositionPartials.cpp")
setup_custom_test_program(test_AngularPositionPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_AngularPositionPartials tudat_observation_partials_test_functions  ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_PositionPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestPositionPartials.cpp")
setup_custom_test_program(test_PositionPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_PositionPartials tudat_observation_partials_test_functions ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EulerAngleObservationPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestEulerAngleObservationPartials.cpp")
setup_custom_test_program(test_EulerAngleObservationPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_EulerAngleObservationPartials tudat_observation_partials_test_functions ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RotationMatrixPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestRotationMatrixPartials.cpp")
setup_custom_test_program(test_RotationMatrixPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_RotationMatrixPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_FullRotationPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestFullRotationPartials.cpp")
setup_custom_test_program(test_FullRotationPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_FullRotationPartials ${TUDAT_PROPAGATION_LIBRARIES} ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_LightTimePartials "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}/UnitTests/unitTestLightTimePartials.cpp")
setup_custom_test_program(test_LightTimePartials "${SRCROOT}${LIGHTTIMECORRECTIONPARTIALSDIR}")
target_link_libraries(test_LightTimePartials tudat_observation_partials_test_functions ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
