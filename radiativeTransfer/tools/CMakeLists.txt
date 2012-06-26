PROJECT(digitalSnow-tools)

cmake_minimum_required(VERSION 2.4)


FIND_PACKAGE(DGtal 0.5 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


SET(SRCS_Tools
  resizeDCRF
  volSubSample
  Noff2Pbrt)


FOREACH(FILE ${SRCS_Tools})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO)
  add_test(${FILE} ${FILE})
ENDFOREACH(FILE)
