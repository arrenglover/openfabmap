# This Cmake file written by Michael Warren,
# Queensland University of Technology, Australia
# https://wiki.qut.edu.au/display/cyphy/Michael+Warren
# Last updated 2014-11-03

project(openFABMAP)

cmake_minimum_required(VERSION 3.0)

## Cmake setup #################################################################

# make a lib directory in the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# tell cmake that the library goes in the library directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# make a binary directory in the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# tell cmake that the binaries goes in the binary directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Compiler warning level /W3 (msvc) or -Wall (gcc)
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W3" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()


## Required Packages ###########################################################

# OpenMP speedups
message(STATUS "")
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message("Found OpenMP")
endif(OPENMP_FOUND)

# Find OpenCV.
# If it's not found, set OpenCV_DIR to the directory with OpenCVConfig.cmake
if(WIN32)
    set(OpenCV_PATHS
        $ENV{OPENCV_HOME}
        $ENV{OPENCV_DIR}/../../
        C:/opencv/
        C:/OpenCV2.2/
        C:/OpenCV2.3/
        C:/OpenCV2.4/
        )
else() # Linux
    set(OpenCV_PATHS
        $ENV{OPENCV_HOME}/build
        /usr/local/share/OpenCV/
        /usr/share/OpenCV
        )
endif()
find_package(OpenCV REQUIRED NO_MODULE PATHS ${OpenCV_PATHS})

message("OpenCV version: ${OpenCV_VERSION}")

if(OpenCV_VERSION VERSION_LESS "2.4.0")
    add_definitions(-DUSENONFREE)
elseif(OpenCV_VERSION VERSION_LESS "3.0.0")
    if("${OpenCV_LIBRARIES}" MATCHES opencv_nonfree)
        message("Non-free packages found. Using e.g. SIFT/SURF")
        add_definitions(-DUSENONFREE)
    endif()
else()
    if("${OpenCV_LIBRARIES}" MATCHES opencv_xfeatures2d)
        message("Non-free packages found. Using e.g. SIFT/SURF")
        add_definitions(-DUSENONFREE)
    endif()
endif()


## openFABMAP library ##########################################################

# List sources
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src OPENFABMAP_FILES)

# Include the headers
file(GLOB OPENFABMAP_IMPL_INCS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
file(GLOB OPENFABMAP_INCS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

# tell cmake about the library
add_library(openFABMAP ${OPENFABMAP_FILES}
    ${OPENFABMAP_IMPL_INCS} ${OPENFABMAP_INCS})

# Tell CMake where the headers are
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${OpenCV_INCLUDE_DIRS})

# Link against the required libraries in OpenCV >=2.2
target_link_libraries(openFABMAP ${OpenCV_LIBRARIES})

## openFABMAPcli executable ####################################################

# Tell the project where settings / doxygen / readme files are (for Qt Creator)
file(GLOB SETTINGS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/samples/settings.yml")
file(GLOB DOXY_FILES "${CMAKE_CURRENT_SOURCE_DIR}"
    ".travis.yml" "doxygen/*.dox" "doxygen/*.dox.*" "doxygen/*.xml")
set(README_FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

# Copy the settings file across when building (not used for now)
#FILE(COPY ${SETTINGS_FILE} DESTINATION ${CMAKE_BINARY_DIR}/bin)

# Tell cmake about the binary
add_executable(openFABMAPcli ${CMAKE_SOURCE_DIR}/samples/openFABMAPcli.cpp
    ${SETTINGS_FILE} ${DOXY_FILES} ${README_FILES})

# Tell openFABMAPcli to link against its required libs
target_link_libraries(openFABMAPcli openFABMAP ${OpenCV_LIBRARIES} )


## Doxygen API documentation ###################################################

# Add the 'doc' target to generate API documentation with Doxygen
set(DOXYGEN_CREATE_DOCS false CACHE BOOL "Create documentation with Doxygen")
if(DOXYGEN_CREATE_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        # Process and copy configuration file
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        # Uncomment ALL to make every time. Otherwise use "make doc"
        add_custom_target(doc #ALL
            ${DOXYGEN_EXECUTABLE}
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
    else(DOXYGEN_FOUND)
        message(WARNING "Doxygen package not found, no documentation target")
    endif(DOXYGEN_FOUND)
endif(DOXYGEN_CREATE_DOCS)
