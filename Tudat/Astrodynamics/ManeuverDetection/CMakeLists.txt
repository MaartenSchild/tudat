#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(MANEUVERDETECTION_SOURCES
  "${SRCROOT}${MANEUVERDETECTIONDIR}/maneuverDetectionTools.cpp"
)

# Set the header files.
set(MANEUVERDETECTION_HEADERS
  "${SRCROOT}${MANEUVERDETECTIONDIR}/maneuverDetectionTools.h"
)

# Add static libraries.
# Add static libraries.
add_library(tudat_maneuver_detection STATIC ${MANEUVERDETECTION_SOURCES} ${MANEUVERDETECTION_HEADERS})
setup_tudat_library_target(tudat_maneuver_detection "${SRCROOT}${MANEUVERDETECTIONDIR}")
# Add unit tests.
include(tudatLinkLibraries)
