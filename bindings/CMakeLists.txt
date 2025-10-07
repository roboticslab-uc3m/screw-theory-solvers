if((NOT Python3_FOUND OR NOT pybind11_FOUND) AND (NOT DEFINED CREATE_BINDINGS_PYTHON OR CREATE_BINDINGS_PYTHON))
    message(WARNING "Python3/pybind11 package not found, disabling bindings")
endif()

cmake_dependent_option(CREATE_BINDINGS_PYTHON "Enable/disable Python bindings" OFF
                       "Python3_FOUND;pybind11_FOUND" OFF)

if(CREATE_BINDINGS_PYTHON)
    set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version used by PyBind11")

    set(PYTHON_PACKAGES_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/dist-packages"
        CACHE STRING "Install location of the python package")

    pybind11_add_module(PySTS PySTS/PySTS.h
                              PySTS/PySTS.cpp
                              PySTS/solvers.cpp)

    target_link_libraries(PySTS PRIVATE ROBOTICSLAB::ScrewTheoryLib)

    install(TARGETS PySTS
            DESTINATION "${PYTHON_PACKAGES_INSTALL_DIR}")
else()
    set(CREATE_BINDINGS_PYTHON OFF CACHE BOOL "Enable/disable Python bindings" FORCE)
endif()
