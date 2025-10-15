cmake_minimum_required(VERSION 3.20)

project(genesis_physics VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 20)

add_subdirectory(include)

file(GLOB_RECURSE SRC_FILES src/*.cpp src/*.hpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
    src/
    include/
)

target_link_libraries(${PROJECT_NAME} PRIVATE
)