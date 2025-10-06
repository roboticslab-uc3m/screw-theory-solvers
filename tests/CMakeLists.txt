if(NOT GTestSources_FOUND AND (NOT DEFINED ENABLE_tests OR ENABLE_tests))
    message(WARNING "GTestSources package not found, disabling tests")
endif()

cmake_dependent_option(ENABLE_tests "Enable/disable unit tests" ON
                       GTestSources_FOUND OFF)

if(ENABLE_tests)

    add_subdirectory(${GTestSources_SOURCE_DIR} ${CMAKE_BINARY_DIR}/googletest)

    include_directories(${GTestSources_INCLUDE_DIR})

    # testScrewTheory

    if(ENABLE_ScrewTheoryLib)
        add_executable(testScrewTheory testScrewTheory.cpp)

        target_link_libraries(testScrewTheory ROBOTICSLAB::ScrewTheoryLib
                                              gtest_main)

        gtest_discover_tests(testScrewTheory)
    endif()

else()

    set(ENABLE_tests OFF CACHE BOOL "Enable/disable unit tests" FORCE)

endif()
