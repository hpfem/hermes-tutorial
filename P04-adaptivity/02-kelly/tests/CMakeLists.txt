project(test-P04-02-kelly)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

add_test(${PROJECT_NAME} ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
