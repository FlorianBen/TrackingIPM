find_path(BLOSC_INCLUDE_DIR blosc.h)

find_library(BLOSC_LIBRARIES NAMES blosc)

if(BLOSC_INCLUDE_DIR AND BLOSC_LIBRARIES)
    set(BLOSC_FOUND TRUE)
    message(STATUS "Found BLOSC library: ${BLOSC_LIBRARIES}")
endif()