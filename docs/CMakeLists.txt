# ##############################################################################
# Cmake file for the documentation
# ##############################################################################
# Author: Florian Benedetti 
# Date: 2021/02/02

find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "README.md")
set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/docs/")
set(DOXYGEN_BUILTIN_STL_SUPPORT TRUE)
set(DOXYGEN_GENERATE_XML TRUE)
set(DOXYGEN_EXCLUDE_PATTERNS "nanoflann.hpp" "csv.hpp")

doxygen_add_docs(libSpaceChargeDoc ${CMAKE_SOURCE_DIR}/include
                 "${CMAKE_SOURCE_DIR}/README.md")

set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

set(SPHINX_SOURCE ${CMAKE_SOURCE_DIR}/docs)

add_custom_command(
  OUTPUT ${SPHINX_INDEX_FILE}
  COMMAND ${SPHINX_EXECUTABLE} -b html
          -Dbreathe_projects.SpaceCharge=${CMAKE_CURRENT_BINARY_DIR}/xml
          ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${SPHINX_SOURCE}/index.rst libSpaceChargeDoc
  MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
  COMMENT "Generating documentation with Sphinx")

add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})