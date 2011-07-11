# Locate OpenCv-0.9.9 install directory

# This module defines
# OPENCV_HOME where to find include, lib, bin, etc.

# Find all the opencv stuff with pkg-config

FIND_PATH( OPENCV_PATH cv.h
# installation selected by user
$ENV{OPENCV_HOME}/include
# system placed in /usr/local/include
/usr/local/include/opencv
# system placed in /usr/include
/usr/include/opencv
# system placed in windows
C:/OpenCV2.1/
)

if( OPENCV_PATH )
    MESSAGE( STATUS "Looking for OpenCV - found")
	MESSAGE( STATUS "OpenCV path: "${OPENCV_PATH} )
    INCLUDE_DIRECTORIES( ${OPENCV_PATH}/include/opencv )
    INCLUDE_DIRECTORIES( ${OPENCV_PATH})
	IF(WIN32)
		LINK_DIRECTORIES( ${OPENCV_PATH}/lib/$(ConfigurationName))
		LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${OPENCV_PATH}/lib)
	ELSE(WIN32)
		LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${OPENCV_PATH}/lib)
	ENDIF(WIN32)
    SET ( OPENCV_FOUND 1 )
else( OPENCV_PATH )
    message( STATUS "Looking for OpenCV  - not found" )
    SET ( OPENCV_FOUND 0 )
endif( OPENCV_PATH )
