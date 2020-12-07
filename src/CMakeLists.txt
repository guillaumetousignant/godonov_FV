add_library(FVM "")
add_subdirectory(entities)
add_subdirectory(fluxes)
add_subdirectory(functions)
add_subdirectory(limiters)
add_subdirectory(solvers)

target_include_directories(FVM PUBLIC 
    $<INSTALL_INTERFACE:include/FVM/>    
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>)

target_compile_features(FVM PRIVATE cxx_std_20)
if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set_property(TARGET FVM PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(FVM PUBLIC optimized OpenMP::OpenMP_CXX)
endif()

install(TARGETS FVM 
    EXPORT FVM-export
    DESTINATION lib)

if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_INSTALL_PREFIX)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ 
        DESTINATION include/FVM
        FILES_MATCHING PATTERN "*.h")
endif()

install(EXPORT FVM-export
  FILE
    FVMConfig.cmake
  NAMESPACE
    FVM::
  DESTINATION
    lib/cmake/FVM
)

add_executable(flow_solver
    main.cpp)

target_link_libraries(flow_solver PRIVATE FVM)

target_compile_features(flow_solver PRIVATE cxx_std_20)
if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set_property(TARGET flow_solver PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()

install(TARGETS flow_solver DESTINATION bin)