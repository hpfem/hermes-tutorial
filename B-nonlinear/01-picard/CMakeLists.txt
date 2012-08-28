project(B-01-picard)
add_executable(${PROJECT_NAME} definitions.cpp main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")
IF(EXISTS tests)
  add_subdirectory(tests)
  enable_testing()
ENDIF(EXISTS tests)

