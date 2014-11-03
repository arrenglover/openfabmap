# This Cmake file written by Kirk MacTavish
# Institute for Aerospace Studies, University of Toronto
# http://http://asrl.utias.utoronto.ca/~kam/

# FindOpenCV.cmake: Locate OpenCV >=2.2 headers and libs (for Windows and Linux)

# This module defines
# OpenCV_DIR the location of the OpenCV 2.2 cmake config file.

IF(WIN32)

	FIND_PATH( OpenCV_DIR OpenCVConfig.cmake
		$ENV{OPENCV_HOME}
		$ENV{OPENCV_DIR}/../../
		C:/opencv/
		C:/OpenCV2.2/
		C:/OpenCV2.3/
		C:/OpenCV2.4/
	)
	
	if( OpenCV_DIR )
		include(${OpenCV_DIR}/OpenCVConfig.cmake)
		MESSAGE( STATUS "Found OpenCV2.2+: " ${OpenCV_DIR} )
	else( OpenCV_DIR )
		message( STATUS "Could not find OpenCV2.2 or greater - set OpenCV_DIR manually" )
		set(OPENCV_FOUND 0)
	endif( OpenCV_DIR )

ELSE(WIN32) # Linux

	FIND_PATH( OpenCV_DIR OpenCVConfig.cmake
	# installation selected by user
	$ENV{OPENCV_HOME}/build
	)
	
	if( OpenCV_DIR )
		include(${OpenCV_DIR}/OpenCVConfig.cmake)
		MESSAGE( STATUS "Found OpenCV2.2+: " ${OpenCV_DIR} )
	else( OpenCV_DIR )
		message( STATUS "Could not find OpenCV2.2 or greater - set OpenCV_DIR manually" )
		set(OPENCV_FOUND 0)
	endif( OpenCV_DIR )

ENDIF(WIN32)
