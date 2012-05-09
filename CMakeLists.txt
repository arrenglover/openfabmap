PROJECT(openFABMAP)

cmake_minimum_required(VERSION 2.6)


############ Cmake setup ##################

# tell cmake where the FindXXX.cmake functions are
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

# make a lib directory in the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# tell cmake that the library goes in the library directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# make a binary directory in the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# tell cmake that the binaries goes in the binary directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

# required packages
FIND_PACKAGE(OpenCV REQUIRED)

############ end CMake setup ##################


############ openFABMAP library ##################
# tell cmake about the library
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src OPENFABMAP_FILES) 

# include the headers
file(GLOB OPENFABMAP_INCS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")

# tell cmake about the library
ADD_LIBRARY(openFABMAP ${OPENFABMAP_FILES} ${OPENFABMAP_INCS})
############ end openFABMAP library ##################

# Tell CMake where the headers are
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link against the required libraries in OpenCV >=2.2
IF(OPENCV2_FOUND)
    SET( OPENCV2_INSTALLED 1)
    message( STATUS "Using OpenCV >= 2.2" )
	IF(WIN32)
		FOREACH (LIB ${OPENCV2_DEBUG_LIBS})
			TARGET_LINK_LIBRARIES(openFABMAP
			debug ${LIB})
		ENDFOREACH(LIB)
		FOREACH (LIB ${OPENCV2_RELEASE_LIBS})
			TARGET_LINK_LIBRARIES(openFABMAP
			optimized ${LIB})
		ENDFOREACH(LIB)
	ELSE(WIN32)
		TARGET_LINK_LIBRARIES(openFABMAP opencv_highgui opencv_core opencv_features2d)
	ENDIF(WIN32)
ENDIF(OPENCV2_FOUND)


############ openFABMAPcli executable ##################
if(OPENCV2_FOUND)

	# tell the project where the settings file is
	file(GLOB SETTINGS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/samples/settings.yml")
	FILE(COPY ${SETTINGS_FILE} DESTINATION ${CMAKE_BINARY_DIR}/bin)

    # tell cmake about the binary
    ADD_EXECUTABLE(openFABMAPcli ${CMAKE_SOURCE_DIR}/samples/openFABMAPcli.cpp ${SETTINGS_FILE})
	
	# Tell openFABMAPcli to link against its required libs
    TARGET_LINK_LIBRARIES(openFABMAPcli openFABMAP)
	
endif(OPENCV2_FOUND)

############ end openFABMAP executable ##################


