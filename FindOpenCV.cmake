# Locate OpenCV >=2.2 install directory

# This module defines
# OPENCV2_FOUND whether the OpenCV 2.2 was found

# OPENCV2_PATH where the OpenCV 2.2 or greater files are (WIN32 only)

# OPENCV2_INCLUDE_PATH where the OpenCV 2.2 or greater header files are

# OPENCV2_LIB_PATH where the OpenCV 2.2 or greater library files are

# OPENCV2_RELEASE_LIBS the list of OpenCV 2.2 or greater release version libs (WIN32 only)
# OPENCV2_DEBUG_LIBS the list of OpenCV 2.2 or greater debug version libs (WIN32 only)

IF(WIN32)

	FIND_PATH( OPENCV2_PATH include/opencv.hpp
	$ENV{OPENCV_HOME}
	C:/OpenCV2.2/
	C:/OpenCV2.3/
	)
	
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
	FIND_PATH( OPENCV2_INCLUDE_PATH opencv.hpp
	# installation selected by user
	$ENV{OPENCV_HOME}/include
	# system placed in /usr/local/include
	/usr/local/include/opencv2
	# system placed in /usr/include
	/usr/include/opencv2
	)
	
	if( OPENCV2_INCLUDE_PATH )
		MESSAGE( STATUS "Looking for OpenCV2.2 or greater - found")
		MESSAGE( STATUS "OpenCV2.2 include path: "${OPENCV2_INCLUDE_PATH} )
		SET ( OPENCV2_FOUND 1 )
	else( OPENCV2_INCLUDE_PATH )
		message( STATUS "Looking for OpenCV2.2 or greater  - not found" )
		SET ( OPENCV2_FOUND 0 )
	endif( OPENCV2_INCLUDE_PATH )

	
ENDIF(WIN32)
IF(OPENCV2_FOUND)
		INCLUDE_DIRECTORIES( ${OPENCV2_INCLUDE_PATH})
ENDIF(OPENCV2_FOUND)
