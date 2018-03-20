add_definitions(-std=c++17)

add_subdirectory(math_h)
include_directories(${MATH_H_INC})

add_executable(he3gg he3gg.cpp)
