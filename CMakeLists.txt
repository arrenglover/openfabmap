PROJECT(openFABMAP)

cmake_minimum_required(VERSION 2.6)

############ ANN library ##################
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/ANN ANN_FILES)
# tell cmake about the library
ADD_LIBRARY(ANN STATIC ${ANN_FILES})




############ ANN library ##################
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/OpenSURF OPENSURF_FILES)
# tell cmake about the library
ADD_LIBRARY(OpenSURF STATIC ${OPENSURF_FILES})



############ openFABMAP library ##################
# tell cmake about the library
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP OPENFABMAP_FILES) 

# add the headers to the project
#file(GLOB_RECURSE INCS "*.h")

# tell cmake about the library
ADD_LIBRARY(openFABMAP STATIC ${OPENFABMAP_FILES})

# create symbolic lib target for calling target lib_XXX
#ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )

# tell cmake where the FindXXX.cmake functions are
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# make a lib directory in the build directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# tell cmake that the library goes in the library directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# required packages
FIND_PACKAGE(OpenCV REQUIRED)
#FIND_PACKAGE(ANN REQUIRED)

IF(OPENCV_FOUND)
    IF(WIN32)
        TARGET_LINK_LIBRARIES(openFABMAP highgui210 cv210 cxcore210)
	TARGET_LINK_LIBRARIES(OpenSURF highgui210 cv210 cxcore210)
    ELSE(WIN32)
        TARGET_LINK_LIBRARIES(openFABMAP highgui)
	TARGET_LINK_LIBRARIES(OpenSURF highgui)
    ENDIF(WIN32)
ENDIF(OPENCV_FOUND)

#IF(ANN_FOUND)
    TARGET_LINK_LIBRARIES(openFABMAP ANN)
#ENDIF(ANN_FOUND)


# make an include directory in the build directory
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
# get all the headers
#file(GLOB CMAKE_HEADER_LIST ${CMAKE_SOURCE_DIR}/openFABMAP/ *.h)
# copy all the headers
#foreach(_file ${CMAKE_HEADER_LIST})
#    configure_file(_file ${CMAKE_CURRENT_BINARY_DIR}/include/tst.h COPYONLY)
#endforeach
#copy(${CMAKE_SOURCE_DIR}/include/global.h ${CMAKE_CURRENT_BINARY_DIR}/include/global.h)

#MESSAGE(STATUS ${CMAKE_INCLUDE_PATH})
#set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_CURRENT_BINARY_DIR}/include/)
#MESSAGE(STATUS ${CMAKE_INCLUDE_PATH})

############ end openFABMAP library ##################







############ openFABMAP executable ##################
if(OPENCV_FOUND)
    # make a binary directory in the build directory
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

    # tell cmake that the binary goes in the binary directory
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

    # tell cmake about the binary
    ADD_EXECUTABLE(exampleopenFABMAP ${CMAKE_SOURCE_DIR}/exampleopenFABMAP/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP/ConfigFile/ConfigFile.cpp)

    # tell the executable where the fabmap headers are. Need to put these somewhere better like /usr/local/include/openfabmap
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP)

    #AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/openFABMAP ALL_FILES) 

    TARGET_LINK_LIBRARIES(exampleopenFABMAP openFABMAP ANN OpenSURF)
endif(OPENCV_FOUND)

############ end openFABMAP executable ##################

