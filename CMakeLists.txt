# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(VESPA)

# Check if we have a fortran compiler
find_program(FORT NAMES $ENV{FC}
                        gfortran
                        ifort)
if(NOT FORT)
  project(VESPA C CXX)
  message(STATUS
          "No FORTRAN compiler")
else()
  project(VESPA C CXX Fortran)
  message(STATUS
          "Found FORTRAN compiler")
endif()

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${FISSION_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(VESPA VESPA.cc ${sources} ${headers})
target_link_libraries(VESPA ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(VESPA_SCRIPTS
  icons.mac
  gui.mac
  init.mac
  init_vis.mac
  vis.mac
  example-Efficiency-LS301-4inch.mac
  )

foreach(_script ${VESPA_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX

install(TARGETS VESPA DESTINATION bin)

# Check if we override built-in version with new version (libFission.*) of the library
if(USE_FISSION_NEW)
  find_library(LFISSION Fission
               PATHS ${CMAKE_SOURCE_DIR}/../lib ${FISSION_DIR}/lib
               DOC "fission library"
               NO_DEFAULT_PATH
              )
  if(NOT LFISSION)
    message(FATAL_ERROR
            "fission library not found")
  else()
    target_link_libraries(VESPA 
                          ${LFISSION} 
                         )
    configure_file(${LFISSION} .
                  COPYONLY
                  )
    message(STATUS ${LFISSION})
  endif()
  message(STATUS
          "Overriding fission library libFission.* with version in ${FISSION_DIR}")
  include_directories(${CMAKE_SOURCE_DIR}/../include ${FISSION_DIR}/include /Users/alfgook/geant4/fission_v2.0.5/include)
  add_definitions(-DFISSION_NEW)

  # Check if we use FREYA 
  if(USE_FREYA)
    if(NOT FORT)
      message(FATAL_ERROR
              "No Fortran compiler found, so FREYA is not included")
    endif()
    message(STATUS
            "Using FREYA")
    add_definitions(-DUSEFREYA)
  endif()
endif()