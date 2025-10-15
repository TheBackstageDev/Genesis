cmake_minimum_required(VERSION 3.20)

project(genesis_physics VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.2
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
    SYSTEM)
FetchContent_MakeAvailable(SFML)

add_subdirectory(include)

file(GLOB_RECURSE SRC_FILES src/*.cpp src/*.hpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
    src/
    include/
    include/sfml/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    SFML::Graphics
)