##===------------------------------------------------------------------------------*- CMake -*-===##
##                          _
##                         | |
##                       __| | __ ___      ___ ___
##                      / _` |/ _` \ \ /\ / / '_  |
##                     | (_| | (_| |\ V  V /| | | |
##                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
##
##
##  This file is distributed under the MIT License (MIT).
##  See LICENSE.txt for details.
##
##===------------------------------------------------------------------------------------------===##
include(GoogleTest)

set(executable ${PROJECT_NAME}UnittestCodeGen)
add_executable(${executable} TestCodeGen.cpp TestMain.cpp)
target_add_dawn_standard_props(${executable})
target_link_libraries(${executable} DawnOptimizer DawnCompiler DawnSupport DawnUnittest gtest)
file(COPY input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
gtest_discover_tests(${executable} TEST_PREFIX "Dawn::Unit::CodeGen::"
  EXTRA_ARGS ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CXX_COMPILER} ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${CMAKE_CUDA_COMPILER}
  DISCOVERY_TIMEOUT 30)
