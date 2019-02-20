cmake_minimum_required(VERSION 3.5)

project(freel)

set(VERSION 2.0)

enable_language(Fortran)


set(CMAKE_C_COMPILER gcc)
set (CMAKE_C_FLAGS "-g")
set(CMAKE_BUILD_TYPE Debug)

set(TARGET freel)

file(GLOB SOURCE_FILES src/*.f90)



add_executable(${TARGET} ${SOURCE_FILES})
#set(EXECUTABLE_OUTPUT_PATH ../.)
