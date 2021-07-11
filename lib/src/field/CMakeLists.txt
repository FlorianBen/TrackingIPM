# ##############################################################################
# CMake file for the io lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_FIELD
    ${CMAKE_CURRENT_SOURCE_DIR}/bunch.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/field_map.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/field_fdm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fields.cpp)

set(HEADER_FIELD
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/bunch.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_map.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_comsol.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_map.hpp)

if(USE_DEALII_FEM)
  list(APPEND LIB_FIELD ${CMAKE_CURRENT_SOURCE_DIR}/field_fem.cpp)
  list(APPEND HEADER_FIELD
     ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_fem.hpp)
endif()

install(FILES ${HEADER_FIELD}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/field/)

# Dependencies
find_package(TBB REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
if(USE_DEALII_FEM)
  find_package(MPI REQUIRED)
  find_package(deal.II REQUIRED)
endif()

add_library(field SHARED ${HEADER_FIELD} ${LIB_FIELD})
add_library(SpaceCharge::field ALIAS field)

target_link_libraries(
  field
  PUBLIC core
  PUBLIC TBB::tbb
  PUBLIC OpenMP::OpenMP_CXX
  PUBLIC GSL::gsl)

if(USE_DEALII_FEM)
  target_link_libraries(field PUBLIC MPI::MPI_CXX)
  DEAL_II_SETUP_TARGET(field)
  target_include_directories(field PUBLIC ${DEAL_II_INCLUDE_DIRS})
endif()

target_include_directories(
  field PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
               $<INSTALL_INTERFACE:include>)

target_compile_features(io PUBLIC cxx_std_17)
