//===--------------------------------------------------------------------------------*- C++ -*-===//
//                          _
//                         | |
//                       __| | __ ___      ___ ___
//                      / _` |/ _` \ \ /\ / / '_  |
//                     | (_| | (_| |\ V  V /| | | |
//                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "dawn/Compiler/DawnCompiler.h"
#include "dawn/Compiler/Options.h"
#include "dawn/IIR/IIR.h"
#include "dawn/IIR/StencilInstantiation.h"
#include "dawn/Optimizer/PassSetCaches.h"
#include "dawn/Serialization/IIRSerializer.h"
#include "dawn/Unittest/CompilerUtil.h"
#include "test/unit-test/dawn/Optimizer/TestEnvironment.h"

#include <fstream>
#include <gtest/gtest.h>

using namespace dawn;

namespace {

class TestPassPreprocessor : public ::testing::Test {
protected:
  dawn::OptimizerContext::OptimizerContextOptions options_;
  std::unique_ptr<OptimizerContext> context_;

  void runTest(const std::string& filename, const bool expected = true) {
    std::shared_ptr<iir::StencilInstantiation> instantiation =
        CompilerUtil::load(filename, options_, context_, TestEnvironment::path_);
    ASSERT_EQ(instantiation != nullptr, expected);
  }
};

TEST_F(TestPassPreprocessor, BracketAccessTest) {
  runTest("input/BracketAccessTest.sir");
}

TEST_F(TestPassPreprocessor, DoMethodTest) {
  runTest("input/DoMethodTest.sir");
}

TEST_F(TestPassPreprocessor, DISABLED_GlobalsTest) {
  runTest("input/GlobalsTest.sir", false);
}

TEST_F(TestPassPreprocessor, StructTest) {
  runTest("input/StructTest.sir");
}

TEST_F(TestPassPreprocessor, VarDeclarationsTest) {
  runTest("input/VarDeclarationsTest.sir");
}

TEST_F(TestPassPreprocessor, VerticalRegionTest) {
  runTest("input/VerticalRegionTest.sir");
}

} // anonymous namespace
