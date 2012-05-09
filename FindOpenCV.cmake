# Locate OpenCv-0.9.9 install directory

# This module defines
# OPENCV_FOUND whether OpenCV 2.1 was found
# OPENCV2_FOUND whether the OpenCV 2.2 was found

# USE_OPENCV_1 true if the application should use OpenCV <=2.1,
# false if the application should use OpenCV >=2.2

# OPENCV_PATH where the OpenCV 2.1 or lesser files are (WIN32 only)
# OPENCV2_PATH where the OpenCV 2.2 or greater files are (WIN32 only)

# OPENCV_INCLUDE_PATH where the OpenCV 2.1 or lesser header files are
# OPENCV2_INCLUDE_PATH where the OpenCV 2.2 or greater header files are

# OPENCV_LIB_PATH where the OpenCV 2.1 or lesser library files are
# OPENCV2_LIB_PATH where the OpenCV 2.2 or greater library files are

# OPENCV_RELEASE_LIBS the list of OpenCV 2.1 or lesser release version libs (WIN32 only)
# OPENCV2_RELEASE_LIBS the list of OpenCV 2.2 or greater release version libs (WIN32 only)
# OPENCV_DEBUG_LIBS the list of OpenCV 2.1 or lesser debug version libs (WIN32 only)
# OPENCV2_DEBUG_LIBS the list of OpenCV 2.2 or greater debug version libs (WIN32 only)

# depending on conditions, sets USE_OPENCV_1
# Find all the opencv stuff with pkg-config

IF(WIN32)
	FIND_PATH( OPENCV_PATH include/opencv/cv.h 
	$ENV{OPENCV_HOME}
	C:/OpenCV2.1/
	)

	FIND_PATH( OPENCV2_PATH include/opencv.hpp
	$ENV{OPENCV_HOME}
	C:/OpenCV2.2/
	C:/OpenCV2.3/
	)
	
	if( OPENCV_PATH )
		MESSAGE( STATUS "Looking for OpenCV - found")
		MESSAGE( STATUS "OpenCV path: "${OPENCV_PATH} )
		SET ( OPENCV_FOUND 1 )
		SET(OPENCV_INCLUDE_PATH ${OPENCV_PATH}/include/opencv)
		SET(OPENCV_LIB_PATH ${OPENCV_PATH}/lib/)
		file(GLOB OPENCV_RELEASE_LIBS "${OPENCV_LIB_PATH}*[0-9][0-9][0-9].lib")
		file(GLOB OPENCV_DEBUG_LIBS "${OPENCV_LIB_PATH}*[0-9][0-9][0-9]d.lib")
	else( OPENCV_PATH )
		message( STATUS "Looking for OpenCV  - not found" )
		SET ( OPENCV_FOUND 0 )
	endif( OPENCV_PATH )
	
	if( OPENCV2_PATH )
		MESSAGE( STATUS "Looking for OpenCV2.2 or greater - found")
		MESSAGE( STATUS "OpenCV2.2 path: "${OPENCV2_PATH} )
		SET ( OPENCV2_FOUND 1 )
		SET(OPENCV2_INCLUDE_PATH ${OPENCV2_PATH}/include/opencv2 ${OPENCV2_PATH}/include)
		SET(OPENCV2_LIB_PATH ${OPENCV2_PATH}/lib/)
		file(GLOB OPENCV2_RELEASE_LIBS "${OPENCV2_LIB_PATH}*[0-9][0-9][0-9].lib")
		file(GLOB OPENCV2_DEBUG_LIBS "${OPENCV2_LIB_PATH}*[0-9][0-9][0-9]d.lib")
	else( OPENCV2_PATH )
		message( STATUS "Looking for OpenCV2.2 or greater  - not found" )
		SET ( OPENCV2_FOUND 0 )
	endif( OPENCV2_PATH )

ELSE(WIN32)
	FIND_PATH( OPENCV_INCLUDE_PATH cv.h
	# installation selected by user
	$ENV{OPENCV_HOME}/include
	# system placed in /usr/local/include
	/usr/local/include/opencv
	# system placed in /usr/include
	/usr/include/opencv
	)

	FIND_PATH( OPENCV2_INCLUDE_PATH opencv.hpp
	# installation selected by user
	$ENV{OPENCV_HOME}/include
	# system placed in /usr/local/include
	/usr/local/include/opencv2
	# system placed in /usr/include
	/usr/include/opencv2
	)
	
	if( OPENCV_INCLUDE_PATH )
		MESSAGE( STATUS "Looking for OpenCV - found")
		MESSAGE( STATUS "OpenCV include path: "${OPENCV_PATH} )
		SET ( OPENCV_FOUND 1 )
	else( OPENCV_INCLUDE_PATH )
		message( STATUS "Looking for OpenCV  - not found" )
		SET ( OPENCV_FOUND 0 )
	endif( OPENCV_INCLUDE_PATH )
	
	if( OPENCV2_INCLUDE_PATH )
		MESSAGE( STATUS "Looking for OpenCV2.2 or greater - found")
		MESSAGE( STATUS "OpenCV2.2 include path: "${OPENCV2_INCLUDE_PATH} )
		SET ( OPENCV2_FOUND 1 )
	else( OPENCV2_INCLUDE_PATH )
		message( STATUS "Looking for OpenCV2.2 or greater  - not found" )
		SET ( OPENCV2_FOUND 0 )
	endif( OPENCV2_INCLUDE_PATH )

	
ENDIF(WIN32)

# Figure out which OpenCV version we should be using. OpenCV2.2 or greater is defined as better
IF(OPENCV_FOUND AND NOT OPENCV2_FOUND)
    SET(USE_OPENCV_1  TRUE)
ELSEIF (OPENCV_FOUND AND OPENCV2_FOUND)
    SET(USE_OPENCV_1  FALSE)
ELSEIF (NOT OPENCV_FOUND AND OPENCV2_FOUND)
    SET(USE_OPENCV_1  FALSE)
ELSE(OPENCV_FOUND AND NOT OPENCV2_FOUND)
    SET(USE_OPENCV_1  FALSE)
ENDIF(OPENCV_FOUND AND NOT OPENCV2_FOUND)

# Include the right directories depending on which OpenCV we want
IF(OPENCV_FOUND AND USE_OPENCV_1)
		INCLUDE_DIRECTORIES( ${OPENCV_INCLUDE_PATH} )
ENDIF(OPENCV_FOUND AND USE_OPENCV_1)

IF(OPENCV2_FOUND AND NOT USE_OPENCV_1)
		INCLUDE_DIRECTORIES( ${OPENCV2_INCLUDE_PATH})
ENDIF(OPENCV2_FOUND AND NOT USE_OPENCV_1)
