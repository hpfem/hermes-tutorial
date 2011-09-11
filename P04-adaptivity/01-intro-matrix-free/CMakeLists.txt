project(P04-01-intro-matrix-free)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)
