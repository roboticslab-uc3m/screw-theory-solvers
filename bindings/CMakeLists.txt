if(NOT pybind11_FOUND AND (NOT DEFINED CREATE_BINDINGS_PYTHON OR CREATE_BINDINGS_PYTHON))
    message(WARNING "pybind11 package not found, disabling tests")
endif()

cmake_dependent_option(CREATE_BINDINGS_PYTHON "Enable/disable Python bindings" ON
                       pybind11_FOUND OFF)

if(CREATE_BINDINGS_PYTHON)
    # The default for Python3_FIND_UNVERSIONED_NAMES is LAST, which
    # means find_package(Python3) will prefer to find /usr/bin/python3.11
    # over /usr/bin/python3.  We actually prefer to find the system version,
    # so set this to FIRST.
    set(Python3_FIND_UNVERSIONED_NAMES FIRST)

    set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version used by PyBind11")

    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

    # Installation path is determined reliably on most platforms using sysconfig.
    if(Python3_VERSION VERSION_GREATER_EQUAL 3.12)
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "import os;import sysconfig;relative_site_packages = sysconfig.get_path('purelib').replace(sysconfig.get_path('data'), '').lstrip(os.path.sep);print(relative_site_packages)"
                        OUTPUT_VARIABLE Python_INSTDIR
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(1,0,prefix=''))"
                        OUTPUT_VARIABLE Python_INSTDIR
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()

    pybind11_add_module(PySTS PySTS/PySTS.h
                              PySTS/PySTS.cpp
                              PySTS/solvers.cpp)

    target_link_libraries(PySTS PRIVATE ROBOTICSLAB::ScrewTheoryLib)

    install(TARGETS PySTS
            DESTINATION "${Python_INSTDIR}")
else()
    set(CREATE_BINDINGS_PYTHON OFF CACHE BOOL "Enable/disable Python bindings" FORCE)
endif()
