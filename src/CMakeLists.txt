add_executable(assignment_3 
    main.cpp)

target_include_directories(assignment_3 PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(assignment_3 PRIVATE cxx_std_20)

install(TARGETS assignment_3 DESTINATION bin)