CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(Aspect 2.4.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

IF (NOT Aspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find a valid ASPECT build/installation directory. "
    "....\n")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "melt_global_mod")
PROJECT(${TARGET})
ADD_LIBRARY(${TARGET} SHARED melt_global_loop.cc)
ASPECT_SETUP_PLUGIN(${TARGET})