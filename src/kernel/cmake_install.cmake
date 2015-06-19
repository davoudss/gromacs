# Install script for directory: /NOBACKUP/davoudss/gromacs-4.6.1/src/kernel

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
   "/NOBACKUP/davoudss/gromacs/share/man/man1/grompp.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/grompp.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/tpbconv.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/tpbconv.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/pdb2gmx.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/pdb2gmx.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/g_protonate.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/g_protonate.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/gmxdump.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/gmxdump.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/g_x2top.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/g_x2top.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/gmxcheck.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/gmxcheck.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/share/man/man1/mdrun.1")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/share/man/man1" TYPE FILE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/man/man1/mdrun.1")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "grompp")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/grompp_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/grompp_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/grompp_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "grompp")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "tpbconv")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/tpbconv_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/tpbconv_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/tpbconv_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "tpbconv")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "pdb2gmx")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/pdb2gmx_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/pdb2gmx_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "pdb2gmx")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "g_protonate")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/g_protonate_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/g_protonate_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_protonate_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "g_protonate")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "gmxdump")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/gmxdump_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/gmxdump_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxdump_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "gmxdump")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "g_x2top")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/g_x2top_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/g_x2top_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_x2top_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "g_x2top")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "gmxcheck")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/gmxcheck_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/gmxcheck_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "gmxcheck")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "g_luck")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/g_luck_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/g_luck_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/g_luck_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "g_luck")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "mdrun")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d"
         RPATH "\$ORIGIN/../lib")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/bin/mdrun_d")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/bin" TYPE EXECUTABLE FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/mdrun_d")
  IF(EXISTS "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d"
         OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel:/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
         NEW_RPATH "\$ORIGIN/../lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/bin/mdrun_d")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "mdrun")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "libraries-gmxpreprocess")
  FOREACH(file
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libgmxpreprocess_d.so.8"
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libgmxpreprocess_d.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHECK
           FILE "${file}"
           RPATH "\$ORIGIN/../lib")
    ENDIF()
  ENDFOREACH()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/lib/libgmxpreprocess_d.so.8;/NOBACKUP/davoudss/gromacs/lib/libgmxpreprocess_d.so")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/lib" TYPE SHARED_LIBRARY FILES
    "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/libgmxpreprocess_d.so.8"
    "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/libgmxpreprocess_d.so"
    )
  FOREACH(file
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libgmxpreprocess_d.so.8"
      "$ENV{DESTDIR}/NOBACKUP/davoudss/gromacs/lib/libgmxpreprocess_d.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/NOBACKUP/davoudss/gromacs-4.6.1/src/mdlib:/NOBACKUP/davoudss/gromacs-4.6.1/src/gmxlib:"
           NEW_RPATH "\$ORIGIN/../lib")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "libraries-gmxpreprocess")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/NOBACKUP/davoudss/gromacs/lib/pkgconfig/libgmxpreprocess_d.pc")
FILE(INSTALL DESTINATION "/NOBACKUP/davoudss/gromacs/lib/pkgconfig" TYPE FILE RENAME "libgmxpreprocess_d.pc" FILES "/NOBACKUP/davoudss/gromacs-4.6.1/src/kernel/libgmxpreprocess.pc")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "development")

