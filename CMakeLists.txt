cmake_minimum_required(VERSION 3.20)
project(Stoner-SCF CXX)

set(CMAKE_CXX_STANDARD 20)

# CPM bootstrap
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.0/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# Eigen3
CPMAddPackage(
  NAME Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  OPTIONS "EIGEN_BUILD_DOC OFF" "BUILD_TESTING OFF"
)

# OpenMP
find_package(OpenMP REQUIRED)

add_executable(Stoner-SCF main.cpp hamiltonian.cpp scf.cpp sweep.cpp)

target_link_libraries(Stoner-SCF
  PRIVATE
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
)

target_compile_options(Stoner-SCF PRIVATE
  -O3
  -march=native
)

# GoogleTest
CPMAddPackage(
  NAME googletest
  GITHUB_REPOSITORY google/googletest
  GIT_TAG v1.15.2
  OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt ON"
)

enable_testing()

add_executable(Stoner-SCF-tests tests.cpp hamiltonian.cpp scf.cpp sweep.cpp)

target_link_libraries(Stoner-SCF-tests
  PRIVATE
    GTest::gtest
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
)

target_compile_options(Stoner-SCF-tests PRIVATE -O2 -march=native -UNDEBUG)

include(GoogleTest)
gtest_discover_tests(Stoner-SCF-tests)