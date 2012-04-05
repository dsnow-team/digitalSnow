PROJECT(Deformations)


#DGtal
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(DGtal REQUIRED)
MESSAGE(STATUS "Found DGtal: ${DGTAL_LIBRARIES}")
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


#Inclusion de Boost
include(FindBoost)
find_package(Boost 1.36.0 REQUIRED program_options)
message(STATUS "Found Boost: ${Boost_LIBRARIES} ")
link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#fftw
INCLUDE_DIRECTORIES(/usr/include/)
LINK_DIRECTORIES(/usr/lib/)

##########################################
ADD_EXECUTABLE(3dVolViewer 3dVolViewer)
TARGET_LINK_LIBRARIES(3dVolViewer ${DGTAL_LIBRARIES} ${BOOST_LIBRARIES})

ADD_EXECUTABLE(testLocalDeformation testLocalDeformation)
TARGET_LINK_LIBRARIES(testLocalDeformation ${DGTAL_LIBRARIES} ${BOOST_LIBRARIES})

ADD_EXECUTABLE(testLocalDeformation2d testLocalDeformation2d)
TARGET_LINK_LIBRARIES(testLocalDeformation2d ${DGTAL_LIBRARIES} ${BOOST_LIBRARIES})

ADD_EXECUTABLE(testDisk testDisk)
TARGET_LINK_LIBRARIES(testDisk ${DGTAL_LIBRARIES} ${BOOST_LIBRARIES} fftw3)

ADD_EXECUTABLE(deformation2d deformation2d)
TARGET_LINK_LIBRARIES(deformation2d ${DGTAL_LIBRARIES} ${BOOST_LIBRARIES} fftw3)

ADD_EXECUTABLE(deformation3d deformation3d)
TARGET_LINK_LIBRARIES(deformation3d ${DGTAL_LIBRARIES} ${BOOST_LIBRARIES} fftw3)

ADD_EXECUTABLE(imageBlurring imageBlurring)
TARGET_LINK_LIBRARIES(imageBlurring ${DGTAL_LIBRARIES} fftw3)
