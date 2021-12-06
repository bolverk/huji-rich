cmake_minimum_required(VERSION 3.20.2)
project(rich)
enable_language(Fortran CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CXX Compiler: " ${CMAKE_CXX_COMPILER})
message("C Compiler: " ${CMAKE_C_COMPILER})
message("Fortran Compiler: " ${CMAKE_Fortran_COMPILER})
message("CXX Standard: " ${CMAKE_CXX_STANDARD})
message("Build Type: " ${CMAKE_BUILD_TYPE})

message("CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
message("CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
message("CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
message("CMAKE_Fortran_FLAGS = " ${CMAKE_Fortran_FLAGS})
message("CMAKE_Fortran_FLAGS_RELEASE = " ${CMAKE_Fortran_FLAGS_RELEASE})
message("CMAKE_Fortran_FLAGS_DEBUG = " ${CMAKE_Fortran_FLAGS_DEBUG})
message("Root_dir = " ${PROJECT_ROOT_DIR})

set(proj_dir "${PROJECT_SOURCE_DIR}")

include_directories("${PROJECT_BINARY_DIR}/generated/cppversion")
find_package(pybind11 REQUIRED PATHS ${PYBIND11} NO_DEFAULT_PATH)
message("pybind11 directory: " ${pybind11_DIR})
message("python exe path: " ${PYTHON_EXECUTABLE})

# files to compile
file(GLOB_RECURSE SRC_CPP FOLLOW_SYMLINKS "${PROJECT_SOURCE_DIR}/*.cpp")
file(GLOB R3D "${PROJECT_ROOT_DIR}/r3d/src/*.c")
file(GLOB TEST_DIR "${TEST_DIR}/*.cpp")

find_library(HDF5 hdf5 PATHS ${HDF5_LIB_DIRECTORY} NO_DEFAULT_PATH)
find_library(HDF5_HL hdf5_hl PATHS ${HDF5_LIB_DIRECTORY} NO_DEFAULT_PATH)
find_library(HDF5_CPP hdf5_cpp PATHS ${HDF5_LIB_DIRECTORY} NO_DEFAULT_PATH)
set(HDF5_LIBS ${HDF5_CPP} ${HDF5_HL} ${HDF5})
if (NOT HDF5_HL)
    message(FATAL_ERROR "NO HDF5_HL")
endif()
if (NOT HDF5)
    message(FATAL_ERROR "NO HDF5")
endif()

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_ROOT_DIR})
include_directories(${BOOST_INCLUDE})
include_directories(${HDF5_INCLUDE})

message("proj_dir ${proj_dir}")
message("HDF5 ${HDF5_LIBS}")

add_executable(${EXE_NAME} ${TEST_DIR} ${R3D} ${SRC_CPP})
target_link_libraries(${EXE_NAME} PRIVATE pybind11::embed pybind11::module ${HDF5_LIBS})
