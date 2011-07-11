# Locate OpenCv-0.9.9 install directory

# This module defines
# ANN_HOME where to find include, lib, bin, etc.

# Find all the ann stuff with pkg-config

FIND_PATH( ANN_PATH ANN.h
# in the openfabmap directory
${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/ANN
# installation selected by user
$ENV{ANN_HOME}/include
# system placed in /usr/local/include
/usr/local/include/ann
# system placed in /usr/include
/usr/include/ann

)

if( ANN_PATH )
    MESSAGE( STATUS "Looking for ANN - found")
    MESSAGE( STATUS "ANN include path: "${ANN_PATH} )
    INCLUDE_DIRECTORIES( ${ANN_PATH} )
    if(${ANN_PATH} MATCHES ${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/ANN)
        LINK_DIRECTORIES(${ANN_PATH})
    else(${ANN_PATH} MATCHES ${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/ANN)
        LINK_DIRECTORIES(${ANN_PATH}/lib)
    endif(${ANN_PATH} MATCHES ${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/ANN)
    SET ( ANN_FOUND 1 )
else( ANN_PATH )
    message( STATUS "Looking for ANN  - not found" )
    SET ( ANN_FOUND 0 )
endif( ANN_PATH )
