# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.14)

# Name of the executable
add_executable(c_tester main.c)

# Link the Fortran library to the C test program
target_link_libraries(c_tester PRIVATE xhcff)

# Include directories
target_include_directories(c_tester PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)

