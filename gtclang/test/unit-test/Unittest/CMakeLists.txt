##===------------------------------------------------------------------------------*- CMake -*-===##
##                         _       _
##                        | |     | |
##                    __ _| |_ ___| | __ _ _ __   __ _
##                   / _` | __/ __| |/ _` | '_ \ / _` |
##                  | (_| | || (__| | (_| | | | | (_| |
##                   \__, |\__\___|_|\__,_|_| |_|\__, | - GridTools Clang DSL
##                    __/ |                       __/ |
##                   |___/                       |___/
##
##
##  This file is distributed under the MIT License (MIT).
##  See LICENSE.txt for details.
##
##===------------------------------------------------------------------------------------------===##
include(GoogleTest)

set(test_name ${PROJECT_NAME}UnittestUnittest)
add_executable(${test_name}
  TestUnittestEnvironment.cpp
  TestMain.cpp
)

target_add_gtclang_standard_props(${test_name})
target_link_libraries(${test_name} ${PROJECT_NAME} ${PROJECT_NAME}Unittest gtest)

set_target_properties(${test_name} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/unittest
)

gtest_discover_tests(${test_name} TEST_PREFIX "GTClang::Unit::Unittest::" DISCOVERY_TIMEOUT 30)
