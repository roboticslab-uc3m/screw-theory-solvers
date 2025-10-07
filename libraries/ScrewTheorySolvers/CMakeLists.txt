option(ENABLE_ScrewTheorySolvers "Enable/disable ScrewTheorySolvers library" ON)

if(ENABLE_ScrewTheorySolvers)

    add_library(ScrewTheorySolvers SHARED ScrewTheoryTools.hpp
                                          ScrewTheoryTools.cpp
                                          MatrixExponential.hpp
                                          MatrixExponential.cpp
                                          ProductOfExponentials.hpp
                                          ProductOfExponentials.cpp
                                          ScrewTheoryIkProblem.hpp
                                          ScrewTheoryIkProblem.cpp
                                          ScrewTheoryIkProblemBuilder.cpp
                                          ScrewTheoryIkSubproblems.hpp
                                          PadenKahanSubproblems.cpp
                                          PardosGotorSubproblems.cpp
                                          ConfigurationSelector.hpp
                                          ConfigurationSelector.cpp
                                          ConfigurationSelectorLeastOverallAngularDisplacement.cpp
                                          ConfigurationSelectorHumanoidGait.cpp)

    set_property(TARGET ScrewTheorySolvers PROPERTY PUBLIC_HEADER MatrixExponential.hpp
                                                                  ProductOfExponentials.hpp
                                                                  ScrewTheoryIkProblem.hpp
                                                                  ConfigurationSelector.hpp)

    target_link_libraries(ScrewTheorySolvers PUBLIC ${orocos_kdl_LIBRARIES})

    target_include_directories(ScrewTheorySolvers PUBLIC ${orocos_kdl_INCLUDE_DIRS}
                                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

    target_compile_features(ScrewTheorySolvers PUBLIC cxx_std_17) # nested namespace definition

    install(TARGETS ScrewTheorySolvers
            EXPORT RL_STS)

    set_property(GLOBAL APPEND PROPERTY _exported_dependencies orocos_kdl)

    add_library(ROBOTICSLAB::ScrewTheorySolvers ALIAS ScrewTheorySolvers)

endif()
