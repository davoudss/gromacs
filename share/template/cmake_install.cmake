# Install script for directory: /NOBACKUP/davoudss/gromacs-4.6.1/share/template

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/NOBACKUP/davoudss/gromacs")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/gromacs/template/CMakeLists.txt")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/gromacs/template" TYPE FILE RENAME "CMakeLists.txt" FILES "/NOBACKUP/davoudss/gromacs-4.6.1/share/template/CMakeLists.txt.template")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/gromacs/template/README;/NOBACKUP/davoudss/gromacs/share/gromacs/template/template.c;/NOBACKUP/davoudss/gromacs/share/gromacs/template/Makefile.pkg")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/gromacs/template" TYPE FILE FILES
    "/NOBACKUP/davoudss/gromacs-4.6.1/share/template/README"
    "/NOBACKUP/davoudss/gromacs-4.6.1/share/template/template.c"
    "/NOBACKUP/davoudss/gromacs-4.6.1/share/template/Makefile.pkg"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/gromacs/template/cmake/FindGROMACS.cmake")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/gromacs/template/cmake" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/share/template/cmake/FindGROMACS.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")

