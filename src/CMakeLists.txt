add_library(flow_solver 
    ExactRiemannFlux_t.cpp
    ExactSolver_t.cpp
    GodonovSolver_t.cpp
    Mesh1D_t.cpp
    RiemannProblem.cpp
    Solver_t.cpp)

target_include_directories(flow_solver PUBLIC 
    $<INSTALL_INTERFACE:include/flow_solver/>    
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>)

target_compile_features(flow_solver PRIVATE cxx_std_20)
if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set_property(TARGET flow_solver PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

install(TARGETS flow_solver 
    EXPORT flow_solver-export
    DESTINATION lib)

if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_INSTALL_PREFIX)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ 
        DESTINATION include/flow_solver
        FILES_MATCHING PATTERN "*.h")
endif()

install(EXPORT flow_solver-export
  FILE
    flow_solverConfig.cmake
  NAMESPACE
    flow_solver::
  DESTINATION
    lib/cmake/flow_solver
)

add_executable(assignment_3 
    main.cpp)

target_link_libraries(assignment_3 PRIVATE flow_solver)

target_compile_features(assignment_3 PRIVATE cxx_std_20)
if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set_property(TARGET assignment_3 PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

install(TARGETS assignment_3 DESTINATION bin)