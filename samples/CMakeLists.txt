# ##############################################################################
# Samples CMake file
# ##############################################################################
# Author: Florian Benedetti Date: 2021/03/15

add_subdirectory(sample_001)
add_subdirectory(sample_002)
add_subdirectory(sample_003)
add_subdirectory(sample_004)
add_subdirectory(sample_005)
if(USE_DEALII_FEM)
  add_subdirectory(sample_007)
endif()
add_subdirectory(sample_008)
add_subdirectory(sample_009)

