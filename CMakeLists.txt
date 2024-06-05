#cmake -DAspect_DIR=/opt/aspect/src/build .; make
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(Aspect 2.4.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

IF (NOT Aspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find a valid ASPECT build/installation directory. "
    "....\n")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "melt_petrol")
PROJECT(${TARGET})
ADD_LIBRARY(${TARGET} SHARED melt_petrol.cc)
ASPECT_SETUP_PLUGIN(${TARGET})


SET(TARGET "c_heating")
PROJECT(${TARGET})
ADD_LIBRARY(${TARGET} SHARED c_heating.cc)
ASPECT_SETUP_PLUGIN(${TARGET})


SET(TARGET "initial_composition")
PROJECT(${TARGET})
ADD_LIBRARY(${TARGET} SHARED initial_composition.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
