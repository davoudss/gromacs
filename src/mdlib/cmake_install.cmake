# Install script for directory: /NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib

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

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "libraries")
  FOREACH(file
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libmd_d.so.8"
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libmd_d.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHECK
           FILE "${file}"
           RPATH "\$ORIGIN/../lib")
    ENDIF()
  ENDFOREACH()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/lib/libmd_d.so.8;/NOBACKUP/davoudss/gromacs/lib/libmd_d.so")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/lib" TYPE SHARED_LIBRARY FILES
    "/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib/libmd_d.so.8"
    "/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib/libmd_d.so"
    )
  FOREACH(file
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libmd_d.so.8"
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libmd_d.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
           NEW_RPATH "\$ORIGIN/../lib")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "libraries")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/lib/pkgconfig/libmd_d.pc")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/lib/pkgconfig" TYPE FILE RENAME "libmd_d.pc" FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib/libmd.pc")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")

