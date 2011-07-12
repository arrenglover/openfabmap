PROJECT(openFABMAP)

cmake_minimum_required(VERSION 2.6)


############ Cmake setup ##################

# tell cmake where the FindXXX.cmake functions are
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

# make a lib directory in the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# tell cmake that the library goes in the library directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# required packages
FIND_PACKAGE(OpenCV REQUIRED)

############ end CMake setup ##################


############ ANN library ##################
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/ANN ANN_FILES)

# include the headers
file(GLOB ANN_INCS "${CMAKE_CURRENT_SOURCE_DIR}/ANN/*.h")

# tell cmake about the library
ADD_LIBRARY(ANN STATIC ${ANN_FILES} ${ANN_INCS})
############ end ANN library ##################


############ OpenSURF library ##################
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/OpenSURF OPENSURF_FILES)

# include the headers
file(GLOB OPENSURF_INCS "${CMAKE_CURRENT_SOURCE_DIR}/OpenSURF/*.h")

# tell cmake about the library
ADD_LIBRARY(OpenSURF STATIC ${OPENSURF_FILES} ${OPENSURF_INCS})
############ end OpenSURF library ##################


############ openFABMAP library ##################
# tell cmake about the library
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP OPENFABMAP_FILES) 

# include the headers
file(GLOB OPENFABMAP_INCS "${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/*.h")

# tell cmake about the library
ADD_LIBRARY(openFABMAP STATIC ${OPENFABMAP_FILES} ${OPENFABMAP_INCS})

# Tell openFABMAP that it needs ANN
TARGET_LINK_LIBRARIES(openFABMAP ANN)
############ end openFABMAP library ##################


# Tell CMake where the headers are
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Link the libs against OpenCV if it was found
IF(OPENCV_FOUND)
    IF(WIN32)
        TARGET_LINK_LIBRARIES(openFABMAP debug highgui210d debug cv210d debug cxcore210d optimized highgui210 optimized cv210 optimized cxcore210)
	TARGET_LINK_LIBRARIES(OpenSURF highgui210 cv210 cxcore210)
    ELSE(WIN32)
        TARGET_LINK_LIBRARIES(openFABMAP highgui)
	TARGET_LINK_LIBRARIES(OpenSURF highgui)
    ENDIF(WIN32)
ENDIF(OPENCV_FOUND)

############ openFABMAP executable ##################
if(OPENCV_FOUND)
    # make a binary directory in the build directory
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

    # tell cmake that the binary goes in the binary directory
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

	# tell the project where the settings file is
	file(GLOB SETTINGS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/settings.txt")

    # tell cmake about the binary
    ADD_EXECUTABLE(openFABMAPexe ${CMAKE_SOURCE_DIR}/openFABMAPexe/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/ConfigFile/ConfigFile.cpp ${CMAKE_CURRENT_SOURCE_DIR}/ConfigFile/ConfigFile.h ${SETTINGS_FILE})

    # tell the executable where the fabmap headers are. Need to put these somewhere better like /usr/local/include/openfabmap
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP)
	
	# Tell openFABMAPexe to link against it's required libs
    TARGET_LINK_LIBRARIES(openFABMAPexe openFABMAP ANN OpenSURF)
	
endif(OPENCV_FOUND)

############ end openFABMAP executable ##################


