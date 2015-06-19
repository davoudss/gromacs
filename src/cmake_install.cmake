# Install script for directory: /NOBACKUP/davoudss/gromacs-4.6.1/src

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man7/gromacs.7")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man7" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man7/gromacs.7")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib/cmake_install.cmake")
  INCLUDE("/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib/cmake_install.cmake")
  INCLUDE("/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/cmake_install.cmake")
  INCLUDE("/NOBACKUP/davoudss/gromacs-4.6.1/src/tools/cmake_install.cmake")
  INCLUDE("/NOBACKUP/davoudss/gromacs-4.6.1/src/ngmx/cmake_install.cmake")
  INCLUDE("/NOBACKUP/davoudss/gromacs-4.6.1/src/contrib/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

