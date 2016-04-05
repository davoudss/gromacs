# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. Example variables are:
#   CPACK_GENERATOR                     - Generator used to create package
#   CPACK_INSTALL_CMAKE_PROJECTS        - For each project (path, name, component)
#   CPACK_CMAKE_GENERATOR               - CMake Generator used for the projects
#   CPACK_INSTALL_COMMANDS              - Extra commands to install components
#   CPACK_INSTALLED_DIRECTORIES           - Extra directories to install
#   CPACK_PACKAGE_DESCRIPTION_FILE      - Description file for the package
#   CPACK_PACKAGE_DESCRIPTION_SUMMARY   - Summary of the package
#   CPACK_PACKAGE_EXECUTABLES           - List of pairs of executables and labels
#   CPACK_PACKAGE_FILE_NAME             - Name of the package generated
#   CPACK_PACKAGE_ICON                  - Icon used for the package
#   CPACK_PACKAGE_INSTALL_DIRECTORY     - Name of directory for the installer
#   CPACK_PACKAGE_NAME                  - Package project name
#   CPACK_PACKAGE_VENDOR                - Package project vendor
#   CPACK_PACKAGE_VERSION               - Package project version
#   CPACK_PACKAGE_VERSION_MAJOR         - Package project version (major)
#   CPACK_PACKAGE_VERSION_MINOR         - Package project version (minor)
#   CPACK_PACKAGE_VERSION_PATCH         - Package project version (patch)

# There are certain generator specific ones

# NSIS Generator:
#   CPACK_PACKAGE_INSTALL_REGISTRY_KEY  - Name of the registry key for the installer
#   CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS - Extra commands used during uninstall
#   CPACK_NSIS_EXTRA_INSTALL_COMMANDS   - Extra commands used during install


SET(CPACK_BINARY_BUNDLE "")
SET(CPACK_BINARY_CYGWIN "")
SET(CPACK_BINARY_DEB "OFF")
SET(CPACK_BINARY_DRAGNDROP "")
SET(CPACK_BINARY_NSIS "OFF")
SET(CPACK_BINARY_OSXX11 "")
SET(CPACK_BINARY_PACKAGEMAKER "")
SET(CPACK_BINARY_RPM "OFF")
SET(CPACK_BINARY_STGZ "ON")
SET(CPACK_BINARY_TBZ2 "OFF")
SET(CPACK_BINARY_TGZ "ON")
SET(CPACK_BINARY_TZ "ON")
SET(CPACK_BINARY_ZIP "")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_COMPONENTS_ALL "")
SET(CPACK_COMPONENT_GROUP_MDRUN_DESCRIPTION "GROMACS executable for running simulations")
SET(CPACK_COMPONENT_GROUP_TOOLS_DESCRIPTION "All GROMACS executable tools")
SET(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
SET(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
SET(CPACK_GENERATOR "STGZ;TGZ;TZ")
SET(CPACK_INSTALL_CMAKE_PROJECTS "/NOBACKUP/davoudss/gromacs-4.6.1;Gromacs;ALL;/")
SET(CPACK_INSTALL_PREFIX "/NOBACKUP/davoudss/gromacs")
SET(CPACK_MODULE_PATH "/NOBACKUP/davoudss/gromacs-4.6.1/cmake;/NOBACKUP/davoudss/gromacs-4.6.1/cmake/Platform")
SET(CPACK_NSIS_DISPLAY_NAME "Gromacs 4.6.1")
SET(CPACK_NSIS_INSTALLER_ICON_CODE "")
SET(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
SET(CPACK_NSIS_PACKAGE_NAME "Gromacs 4.6.1")
SET(CPACK_OUTPUT_CONFIG_FILE "/NOBACKUP/davoudss/gromacs-4.6.1/CPackConfig.cmake")
SET(CPACK_PACKAGE_CONTACT "gmx-users@gromacs.org")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake-2.8/Templates/CPack.GenericDescription.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gromacs - a toolkit for high-performance molecular simulation")
SET(CPACK_PACKAGE_FILE_NAME "Gromacs-4.6.1-Linux")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Gromacs 4.6.1")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "Gromacs 4.6.1")
SET(CPACK_PACKAGE_NAME "Gromacs")
SET(CPACK_PACKAGE_RELOCATABLE "true")
SET(CPACK_PACKAGE_VENDOR "gromacs.org")
SET(CPACK_PACKAGE_VERSION "4.6.1")
SET(CPACK_PACKAGE_VERSION_MAJOR "4")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PROJECT_CONFIG_FILE "/NOBACKUP/davoudss/gromacs-4.6.1/CPackInit.cmake")
SET(CPACK_RESOURCE_FILE_LICENSE "/NOBACKUP/davoudss/gromacs-4.6.1/COPYING")
SET(CPACK_RESOURCE_FILE_README "/NOBACKUP/davoudss/gromacs-4.6.1/admin/InstallInfo.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "/NOBACKUP/davoudss/gromacs-4.6.1/admin/InstallWelcome.txt")
SET(CPACK_SET_DESTDIR "ON")
SET(CPACK_SOURCE_CYGWIN "")
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2;TZ")
SET(CPACK_SOURCE_IGNORE_FILES "\\.isreposource$;\\.git/;\\.gitignore$")
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "/NOBACKUP/davoudss/gromacs-4.6.1;/;/NOBACKUP/davoudss/gromacs-4.6.1/man;man")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/NOBACKUP/davoudss/gromacs-4.6.1/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_TBZ2 "ON")
SET(CPACK_SOURCE_TGZ "ON")
SET(CPACK_SOURCE_TZ "ON")
SET(CPACK_SOURCE_ZIP "OFF")
SET(CPACK_SYSTEM_NAME "Linux")
SET(CPACK_TOPLEVEL_TAG "Linux")
