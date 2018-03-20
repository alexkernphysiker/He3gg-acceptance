add_definitions(-std=c++17)

add_subdirectory(math_h)
include_directories(${MATH_H_INC})
add_subdirectory(FitGen)
include_directories(${FITGEN_INC})
add_subdirectory(WASA-libs)
include_directories(${WASA_LIBS_INC})

add_executable(he3gg he3gg.cpp common.cpp)
target_link_libraries(he3gg kinematics)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/distributions DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
