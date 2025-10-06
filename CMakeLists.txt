cmake_minimum_required(VERSION 3.19...3.30)

project(ROBOTICSLAB_SCREW_THEORY_SOLVERS
        VERSION 0.1.0
        DESCRIPTION "A collection of solvers based on Screw Theory"
        HOMEPAGE_URL "https://robots.uc3m.es/screw-theory-solvers/"
        LANGUAGES CXX)

# Let the user specify a configuration (only single-config generators).
if(NOT CMAKE_CONFIGURATION_TYPES)
    # Possible values.
    set(_configurations Debug Release MinSizeRel RelWithDebInfo)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${_configurations})

    foreach(_conf ${_configurations})
        set(_conf_string "${_conf_string} ${_conf}")
    endforeach()

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING
                 "Choose the type of build, options are:${_conf_string}")

    if(NOT CMAKE_BUILD_TYPE)
        # Encourage the user to specify build type.
        message(STATUS "Setting build type to 'Release' as none was specified.")
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE Release)
    endif()
endif()

# Pick up our cmake modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Hard dependencies.
find_package(orocos_kdl 1.4 REQUIRED)

# Soft dependencies.
find_package(GTestSources 1.8 QUIET)
find_package(pybind11 QUIET)
find_package(Doxygen QUIET)

# Standard installation directories.
include(GNUInstallDirs)

# Control where libraries and executables are placed during the build.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

# Create targets if specific requirements are satisfied.
include(CMakeDependentOption)

# Acknowledge this is a CTest-friendly project.
enable_testing()

# Add main contents.
add_subdirectory(libraries)
add_subdirectory(tests)
add_subdirectory(bindings)
add_subdirectory(doc)

# Collect public exported dependencies.
get_property(_exported_dependencies GLOBAL PROPERTY _exported_dependencies)

if(DEFINED _exported_dependencies)
    list(REMOVE_DUPLICATES _exported_dependencies)
endif()

# Store the package in the user registry.
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)

include(CMakePackageConfigHelpers)
set(_package_name "roboticslab-screw-theory-solvers")

# Generate local-dir <package>-config.cmake.
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in"
                                "${CMAKE_CURRENT_BINARY_DIR}/${_package_name}-config.cmake"
                                INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
                                NO_SET_AND_CHECK_MACRO
                                NO_CHECK_REQUIRED_COMPONENTS_MACRO
                                INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}")

# Generate install-dir <package>-config.cmake.
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in"
                                "${CMAKE_CURRENT_BINARY_DIR}/${_package_name}-config-install.cmake"
                                INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
                                NO_SET_AND_CHECK_MACRO
                                NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Generate <package>-config-version.cmake.
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${_package_name}-config-version.cmake"
                                    VERSION "${PROJECT_VERSION}"
                                    COMPATIBILITY SameMajorVersion)

# Generate <package>-targets.cmake.
export(EXPORT RL_STS_LIBS
        NAMESPACE ROBOTICSLAB::
        FILE ${_package_name}-targets.cmake)

# Install <package>-config.cmake.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${_package_name}-config-install.cmake"
        RENAME ${_package_name}-config.cmake
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# Install <package>-config-version.cmake.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${_package_name}-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# Install <package>-targets.cmake.
install(EXPORT RL_STS_LIBS
        NAMESPACE ROBOTICSLAB::
        FILE ${_package_name}-targets.cmake
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# Add uninstall target.
include(AddUninstallTarget)
